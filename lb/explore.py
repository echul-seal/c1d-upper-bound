import numpy as np,json,sys,time
from scipy.optimize import linprog
from scipy import ndimage
from multiprocessing import Pool
N=320;xe=np.linspace(-.5,.5,N+1);w=1/N
VP=np.array([s*t for t in (0.25,0.5,1,1.5,2,3,4,6,8,12,16,24,32) for s in (1,-1)])
S=np.sinh(np.pi*VP)/(np.pi*VP)
LAP=np.array([(np.exp(2*np.pi*a*xe[1:])-np.exp(2*np.pi*a*xe[:-1]))/(2*np.pi*a)/s for a,s in zip(VP,S)])
AUB=np.vstack([np.full(N,w),LAP])
def sinc(z):return np.sin(np.pi*z)/(np.pi*z)
def cells(z):
    k=2*np.pi*(z.imag-1j*z.real);return (np.exp(k*xe[1:])-np.exp(k*xe[:-1]))/k
def topsum(G,e):
    g=np.sort(G)[::-1].clip(0);m=e/(2*w);k=int(m)
    return 2*g[:k].sum()+2*(m-k)*g[k] if k<N else 2*g.sum()
def lpmax(G,e):
    r=linprog(-G,A_ub=AUB,b_ub=np.r_[e,np.ones(len(VP))],bounds=(0,2),method='highs')
    return -r.fun if r.status==0 else 9.
def thr(z):
    s=sinc(z);C=cells(z);t0=np.angle(s);out=[];v=abs(z.imag)
    for crit in ('D','P0','P','L0'):
        lo,hi=0.0,0.4
        for _ in range(13):
            e=(lo+hi)/2
            if crit=='D':mg=abs(s)-(-np.expm1(-np.pi*v*e)/(np.pi*v) if v>0 else e)*np.exp(np.pi*v)
            elif crit=='P0':mg=abs(s)-topsum((np.exp(-1j*t0)*C).real,e)
            elif crit=='P':mg=max((np.exp(-1j*(t0+d))*s).real-topsum((np.exp(-1j*(t0+d))*C).real,e) for d in np.linspace(-.4,.4,17))
            else:mg=abs(s)-lpmax((np.exp(-1j*t0)*C).real,e)
            lo,hi=(e,hi) if mg>0 else (lo,e)
        out.append(lo)
    return out
if __name__=='__main__':
    t=time.time();us=np.arange(0.75,1.551,0.016);vs=np.arange(0,1.001,0.016)
    Z=[complex(u,v) for v in vs for u in us]
    with Pool(3) as p:T=np.array(p.map(thr,Z,chunksize=20)).reshape(len(vs),len(us),4)
    np.savez('explore.npz',T=T,us=us,vs=vs);res={}
    for j,c in enumerate(('D','P0','P','L0')):
        lo,hi=0.,0.4
        for _ in range(20):
            e=(lo+hi)/2;lab,_=ndimage.label(T[:,:,j]>e);b=lab[0]
            if set(b[(us<1)&(b>0)])&set(b[(us>1)&(b>0)]):lo=e
            else:hi=e
        res[c]={'eta':lo,'bound':1/(1-lo)}
    res['time']=time.time()-t;json.dump(res,open('explore.json','w'));print(res)
