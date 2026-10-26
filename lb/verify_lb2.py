# Independent check (mpmath interval arithmetic) of a C_1d lower-bound certificate. Lemma: if max|f*f|=1 and (int f)^2>=1-eta,
# e=1_[-1/2,1/2]-f*f has 0<=e<=2, int e<=eta, and (f*f)^=sinc-e^=(f^)^2. If on the certified polygon P around z=1 (plus its mirror,
# covered by z->conj z, phi->-phi, e(x)->e(-x)) Re(e^{-i phi}sinc z) > sup_e Re(e^{-i phi}e^(z)) (bathtub over N cells), then
# sinc-t e^ never vanishes there (0<=t<=1), so (f^)^2 has exactly one zero inside: impossible. Hence C_1d >= 1/(1-eta).
import json,sys,time,math
import numpy as np
from fractions import Fraction as Fr
from mpmath import iv
from multiprocessing import get_context
iv.prec=64
fn=sys.argv[1];part=int(sys.argv[2]);npart=int(sys.argv[3]);NS=(4096,32768)
d=json.load(open(fn));eta=Fr(d['eta']);V=[(Fr(a),Fr(b)) for a,b in d['vertices']];K=len(V)-1
assert 0<eta<1 and V[0][1]==0==V[K][1] and all(y>0 for _,y in V[1:K]) and 1<V[0][0]<2 and 0<V[K][0]<1
assert all((V[k][0]-1)*V[k+1][1]-V[k][1]*(V[k+1][0]-1)>0 for k in range(K))
PI=iv.pi
def I(q):return iv.mpf(q.numerator)/q.denominator
def R(a,b):
    a,b=min(a,b),max(a,b);return iv.mpf([I(a).a,I(b).b])
TAB={}
for N in NS:
    m=eta*N/2;jm=int(m);TAB[N]=([2*PI*R(Fr(k,N)-Fr(1,2),Fr(k+1,N)-Fr(1,2)) for k in range(N)],jm,I(m-jm))
xm=-0.5+(np.arange(4096)+0.5)/4096;DD=np.linspace(-1.4,1.4,57);ef=float(eta)*2048;ej=int(ef);ex=ef-ej
def phase(zc):
    s=np.sin(np.pi*zc)/(np.pi*zc);ph=np.angle(s)+DD
    Kx=np.exp(2*np.pi*xm*zc.imag)[None,:]*np.cos(2*np.pi*xm[None,:]*zc.real+ph[:,None])
    Ks=-np.sort(-Kx,axis=1);Ks=Ks.clip(0);g=(np.exp(-1j*ph)*s).real-(Ks[:,:ej].sum(1)+ex*Ks[:,ej])/2048
    return float(ph[g.argmax()])
def lower(U,W,P):
    cP,sP=iv.cos(P),iv.sin(P);eW,emW=iv.exp(PI*W),iv.exp(-PI*W);sh=(eW-emW)/2;ch=(eW+emW)/2
    su,cu=iv.sin(PI*U),iv.cos(PI*U);den=PI*(U**2+W**2)
    return ((cP*(U*su*ch+W*cu*sh)+sP*(U*cu*sh-W*su*ch))/den).a
def upper(U,W,P,N):
    X2,jm,FR=TAB[N];ups=[]
    for x in X2:
        b=(iv.exp(x*W)*iv.cos(x*U+P)).b
        if b>0:ups.append(math.nextafter(float(b),math.inf))
    ups.sort(reverse=True);ups+=[0.0]
    S=iv.mpf(math.nextafter(math.fsum(ups[:jm]),math.inf))+FR*ups[jm]
    return (2*S/N).b
def seg(k):
    st=[(V[k],V[k+1],0)];n=[0,0];mn=9.0;dm=0
    while st:
        p,q,dep=st.pop();U=R(p[0],q[0]);W=R(p[1],q[1])
        P=iv.mpf(phase(complex(float(p[0]+q[0])/2,float(p[1]+q[1])/2)));lo=lower(U,W,P);ok=0
        for i,N in enumerate(NS):
            if i and dep<3:break
            n[i]+=1;up=upper(U,W,P,N)
            if lo>up:ok=1;mn=min(mn,float(lo)-float(up));break
        if ok:dm=max(dm,dep);continue
        if dep>=14:return (k,False,n,float(lo)-float(up),dep)
        c=((p[0]+q[0])/2,(p[1]+q[1])/2);st+=[(p,c,dep+1),(c,q,dep+1)]
    return (k,True,n,mn,dm)
if __name__=='__main__':
    out=open('lb2_%d.log'%part,'w');t=time.time()
    upper(R(V[0][0],V[1][0]),R(V[0][1],V[1][1]),iv.mpf(0.5),4096);print('t4096 %.2fs'%(time.time()-t),file=out,flush=True)
    ks=list(range(part,K,npart));res=[];t=time.time()
    with get_context('fork').Pool(4) as P:
        for r in P.imap_unordered(seg,ks):
            res.append(r);print(len(res),len(ks),r,'t %.0f'%(time.time()-t),file=out,flush=True)
    ok=all(r[1] for r in res);b=1/(1-eta)
    print(fn,'mpmath-iv verifier part %d/%d:'%(part,npart),'PASS' if ok else 'FAIL','segments',len(ks),'boxes',[sum(r[2][i] for r in res) for i in (0,1)],'minmargin %.3g'%min(r[3] for r in res),'time %.0fs'%(time.time()-t),file=out,flush=True)
    if ok:print('=> C_1d >= 1/(1-eta) = %s = %.9f'%(b,float(b)),file=out,flush=True)
