import json,sys,time,numpy as np
from fractions import Fraction as Fr
from flint import arb,acb,ctx
from multiprocessing import get_context
ctx.prec=80
N=16384;h=1.0/N;xs=-0.5+(np.arange(N)+0.5)*h;PI=arb.pi();X=[arb(float(x),h/2) for x in xs];D=np.linspace(-1.4,1.4,141)
G={}
def setup(fn):
    d=json.load(open(fn));eta=Fr(d['eta']);V=[(Fr(a),Fr(b)) for a,b in d['vertices']];K=len(V)-1
    assert 0<eta<1 and V[0][1]==0==V[K][1] and all(y>0 for _,y in V[1:K]) and 1<V[0][0]<2 and 0<V[K][0]<1
    assert all((V[k][0]-1)*V[k+1][1]-V[k][1]*(V[k+1][0]-1)>0 for k in range(K))
    m=eta*N/2;G.update(eta=eta,V=V,K=K,jm=int(m),fr=float(m-int(m)))
def phi(z):
    jm,fr=G['jm'],G['fr'];s=np.sin(np.pi*z)/(np.pi*z);ph=np.angle(s)+D
    Kx=np.exp(2*np.pi*xs*z.imag)[None,:]*np.cos(2*np.pi*xs[None,:]*z.real+ph[:,None])
    Ks=-np.sort(-Kx,axis=1);Ks=Ks.clip(0);T=2*h*(Ks[:,:jm].sum(1)+fr*Ks[:,jm])
    return float(ph[((np.exp(-1j*ph)*s).real-T).argmax()])
def A(q):return arb(q.numerator)/q.denominator
def box(p,q):
    jm,fr=G['jm'],G['fr'];mx=(p[0]+q[0])/2;my=(p[1]+q[1])/2;r=float(max(abs(q[0]-p[0]),abs(q[1]-p[1])))/2*(1+1e-9)+1e-300
    f=phi(complex(float(mx),float(my)));rot=acb(0,-arb(f)).exp()
    Z=acb(A(mx)+arb(0,r),A(my)+arb(0,r));W=acb(0,-2)*PI*Z
    a=(rot*(PI*Z).sin()/(PI*Z)).real;lo=float(a.mid())-float(a.rad())
    up=sorted((max(float(k.mid())+float(k.rad()),0.0) for k in ((rot*(W*x).exp()).real for x in X)),reverse=True)
    return lo-2*h*(sum(up[:jm])+fr*up[jm])-1e-9
def seg(k):
    V=G['V'];st=[(V[k],V[k+1],0)];n=0;mn=9.0
    while st:
        p,q,dep=st.pop();n+=1;g=box(p,q)
        if g>0:
            mn=min(mn,g);continue
        if dep>=10:return (k,False,n,g)
        c=((p[0]+q[0])/2,(p[1]+q[1])/2);st+=[(p,c,dep+1),(c,q,dep+1)]
    return (k,True,n,mn)
if __name__=='__main__':
    for fn in sys.argv[1:]:
        setup(fn);t=time.time()
        with get_context("fork").Pool(4) as P:R=P.map(seg,range(G['K']))
        ok=all(r[1] for r in R);b=1/(1-G['eta'])
        print(fn,'P-criterion arb verifier:','PASS' if ok else 'FAIL','segments',G['K'],'boxes',sum(r[2] for r in R),'minmargin %.3g'%min(r[3] for r in R),'time %.0fs'%(time.time()-t),flush=True)
        if ok:
            print('Certified: C_1d >= 1/(1-eta) = %s = %.9f'%(b,float(b)),flush=True);json.dump(json.load(open(fn)),open('cert_P.json','w'));break
        print(' fail segs',[r[0] for r in R if not r[1]][:12],flush=True)
