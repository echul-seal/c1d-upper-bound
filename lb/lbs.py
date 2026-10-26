import numpy as np,json,sys,time
from multiprocessing import Pool
pi=np.pi;N=2048;h=1.0/N;xs=-0.5+(np.arange(N)+0.5)*h
D=np.linspace(-1.3,1.3,53);ETAS=[float(e) for e in sys.argv[1].split(',')]
def marg(z):
    s=np.sin(pi*z)/(pi*z);ph=np.angle(s)+D
    K=np.exp(2*pi*xs*z.imag)[None,:]*np.cos(2*pi*xs[None,:]*z.real+ph[:,None])
    K=-np.sort(-K,axis=1);K=K.clip(0);C=np.cumsum(K,1);A=(np.exp(-1j*ph)*s).real;out=[]
    for e in ETAS:
        m=e/(2*h);j=int(m);T=2*h*(C[:,j-1]+(m-j)*K[:,j]);out.append((A-T).max()/np.exp(pi*z.imag))
    return out
TH=np.linspace(0,pi,181);RR=np.arange(0.04,1.5,0.01)
if __name__=='__main__':
    t=time.time();Z=[1+r*np.exp(1j*a) for a in TH for r in RR]
    with Pool(4) as p:M=np.array(p.map(marg,Z,chunksize=64)).reshape(len(TH),len(RR),len(ETAS))
    np.save('polarM.npy',M);res=[]
    for q,e in enumerate(ETAS):
        G=M[:,:,q].copy();G[0,RR>=0.97]=-9;G[-1,RR>=0.97]=-9;B=G[0].copy();P=[]
        for a in range(1,len(TH)):
            W=np.full((7,len(RR)),-9.0)
            for s in range(-3,4):
                W[s+3]=np.roll(B,s)
                if s>0:W[s+3,:s]=-9
                if s<0:W[s+3,s:]=-9
            P.append(W.argmax(0)-3);B=np.minimum(W.max(0),G[a])
        r=int(B.argmax());val=float(B[r]);path=[r]
        for a in range(len(TH)-1,0,-1):
            r=int(r-P[a-1][r]);path.append(r)
        path=path[::-1];res.append({'eta':e,'minmarg':val,'path':[float(RR[k]) for k in path]});print(e,val,flush=True)
    json.dump({'th':TH.tolist(),'res':res,'time':time.time()-t},open('lbs.json','w'))
