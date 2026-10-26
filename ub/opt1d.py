import sys,json,time,numpy as np,torch
from scipy.optimize import linprog
from scipy.linalg import toeplitz
torch.set_default_dtype(torch.float64)
a=sys.argv;n=int(a[1]);init=a[2];out=a[3];T=float(a[4]);th=int(a[5]);seed=int(a[6]);pmax=float(a[7]);slp=int(a[8])
torch.set_num_threads(th);t0=time.time()
v0=np.array(json.load(open(init))['v'],float);m=len(v0)
v=v0[np.minimum(((np.arange(n)+0.5)*m/n).astype(int),m-1)].copy()
if seed:v=v*(1+0.003*np.random.default_rng(seed).standard_normal(n))
v*=n/v.sum()
def R(v):
    c=np.convolve(v,v);return 2*len(v)*np.abs(c).max()/v.sum()**2
best=R(v);bv=v.copy()
def save(tag):
    json.dump({'n':n,'ratio':best,'tag':tag,'t':time.time()-t0,'v':bv.tolist()},open(out,'w'))
print('init',best,flush=True);save('init')
L=1<<int(np.ceil(np.log2(2*n)))
def J(x,p):
    X=torch.fft.rfft(x,L);c=torch.fft.irfft(X*X,L)[:2*n-1];q=c.abs();mx=q.max().detach()
    return 2*n*mx*((q/mx)**p).sum()**(1/p)/x.sum()**2
p=1024.0;w=bv.copy()
while p<=pmax and time.time()-t0<T*(0.45 if slp else 1.0):
    x=torch.tensor(w,requires_grad=True)
    opt=torch.optim.LBFGS([x],lr=1,max_iter=1500,history_size=100,line_search_fn='strong_wolfe',tolerance_grad=1e-16,tolerance_change=1e-20)
    def cl():
        opt.zero_grad();f=J(x,p);f.backward();return f
    for _ in range(3):
        opt.step(cl);w=x.detach().numpy().copy();w*=n/w.sum();r=R(w)
        if r<best:best=r;bv=w.copy();save('p%g'%p)
    print('p',p,'R',r,'best',best,'t',round(time.time()-t0),flush=True);p*=2
v=bv.copy();rt=3e-4*np.abs(v).max();it=0
while slp and time.time()-t0<T and rt>1e-14:
    c=np.convolve(v,v);M=np.abs(c).max();B=2*rt*np.abs(v).sum()+n*rt*rt
    A=np.nonzero(np.abs(c)>=M-2.5*B)[0];k=len(A)
    Tm=toeplitz(np.r_[v,np.zeros(n-1)],np.r_[v[0],np.zeros(n-1)])[A]
    Aub=np.zeros((2*k,n+1));Aub[:k,:n]=2*Tm;Aub[k:,:n]=-2*Tm;Aub[:,n]=-1
    t1=time.time()
    res=linprog(np.r_[np.zeros(n),1.0],A_ub=Aub,b_ub=np.r_[-c[A],c[A]],A_eq=np.r_[np.ones(n),0.0][None,:],b_eq=[0.0],bounds=[(-rt,rt)]*n+[(0,None)],method='highs')
    it+=1
    if res.status!=0:
        rt*=0.5;print('lp fail',res.status,flush=True);continue
    d=res.x[:n];r0=R(v);ok=0
    for al in (1.0,0.5,0.25):
        r=R(v+al*d)
        if r<r0:v=v+al*d;ok=1;break
    rt=rt*1.6 if (ok and al==1.0) else rt*0.5
    r=R(v)
    if r<best:best=r;bv=v.copy();save('slp%d'%it)
    if it%5==0:print('slp',it,'k',k,'R',r,'rt',rt,'lp',round(time.time()-t1,1),'t',round(time.time()-t0),flush=True)
save('final');print('final',best,flush=True)
