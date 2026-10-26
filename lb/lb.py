import numpy as np, json
from scipy import ndimage
from fractions import Fraction as F
pi=np.pi
def phi(u,v):
    z=np.hypot(u,v)
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.sqrt(np.sin(pi*u)**2+np.sinh(pi*v)**2)*np.exp(-pi*np.abs(v))/(pi*z)
def psi(v,eta,mode):
    v=np.abs(np.asarray(v,float))
    if mode=='crude': return np.full(v.shape,eta)
    with np.errstate(divide='ignore',invalid='ignore'):
        r=-np.expm1(-pi*v*eta)/(pi*v)
    return np.where(v>0,r,eta)
def grid_eta(mode,du=0.002):
    u=np.arange(-1+du/2,2,du); v=np.arange(0,2.5,du); U,V=np.meshgrid(u,v); P=phi(U,V)
    def conn(eta):
        lab,_=ndimage.label(P>psi(V,eta,mode)); b=lab[0]
        return len(set(b[(u<1)&(b>0)])&set(b[(u>1)&(b>0)]))>0
    lo,hi=0.,0.6
    for _ in range(30):
        mid=(lo+hi)/2
        if conn(mid): lo=mid
        else: hi=mid
    return lo
def star(eta,mode,K=720,dr=0.0025,w=6):
    th=np.linspace(0,pi,K+1); rr=np.arange(dr,1.99,dr); R=len(rr)
    Z=1+rr[None,:]*np.exp(1j*th[:,None]); M=phi(Z.real,Z.imag)-psi(Z.imag,eta,mode)
    M[np.abs(Z)<0.05]=-9; M[0,rr>=1]=-9
    val=M[0].copy(); arg=np.zeros((K+1,R),int); idx=np.arange(R)
    for k in range(1,K+1):
        best=np.full(R,-9.); ba=idx.copy()
        for d in range(-w,w+1):
            src=np.full(R,-9.)
            if d>=0: src[d:]=val[:R-d]
            else: src[:d]=val[-d:]
            b=src>best; best=np.where(b,src,best); ba=np.where(b,idx-d,ba)
        val=np.minimum(best,M[k]); arg[k]=ba
    i=int(np.argmax(val)); p=[i]
    for k in range(K,0,-1): i=arg[k][i]; p.append(i)
    return Z[np.arange(K+1),p[::-1]]
def dense_min(zs,eta,mode,n=40):
    t=np.linspace(0,1,n)[None,:]; P=zs[:-1,None]+(zs[1:,None]-zs[:-1,None])*t
    return (phi(P.real,P.imag)-psi(P.imag,eta,mode)).min()
out={}
for mode in ('crude','refined'):
    g=grid_eta(mode); lo,hi=0.,g+0.02
    for _ in range(22):
        mid=(lo+hi)/2
        if dense_min(star(mid,mode),mid,mode)>1e-6: lo=mid
        else: hi=mid
    out[mode]={'grid':g,'star':lo}
    print(mode,'grid eta* %.6f (bound %.6f)  star eta %.6f (bound %.6f)'%(g,1/(1-g),lo,1/(1-lo)),flush=True)
e=float(np.floor(out['refined']['star']*1e4)/1e4-1e-4); zs=star(e,'refined'); V=[]
for k,z in enumerate(zs):
    x=round(float(z.real),6); y=0.0 if k in (0,len(zs)-1) else max(round(float(z.imag),6),1e-6); V.append([x,y])
zz=np.array([complex(a,b) for a,b in V])
print('cert eta',e,'bound',1/(1-e),'dense margin',dense_min(zz,e,'refined',400),'ends',V[0],V[-1])
json.dump({'eta':str(F(str(round(e,6)))),'vertices':[[str(F(str(a))),str(F(str(b)))] for a,b in V]},open('cert_lb.json','w'))
json.dump(out,open('lb_num.json','w'))
