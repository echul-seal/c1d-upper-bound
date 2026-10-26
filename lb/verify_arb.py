import json,sys,time
from fractions import Fraction as Fr
from flint import arb,acb,ctx
ctx.prec=96
d=json.load(open(sys.argv[1]));eta=Fr(d['eta']);V=[(Fr(a),Fr(b)) for a,b in d['vertices']];K=len(V)-1
assert 0<eta<1 and V[0][1]==0 and V[K][1]==0 and all(y>0 for _,y in V[1:K])
assert 1<V[0][0]<2 and 0<V[K][0]<1
assert all((V[k][0]-1)*V[k+1][1]-V[k][1]*(V[k+1][0]-1)>0 for k in range(K))
PI=arb.pi();E=arb(eta.numerator)/eta.denominator
A=lambda q:arb(q.numerator)/q.denominator
def psi_hi(t):
    if t<=0:return E
    t=A(t);return (1-(-PI*t*E).exp())/(PI*t)
st={'boxes':0,'minmargin':9.0,'maxdepth':0}
def chk(p,q,dep):
    st['boxes']+=1;st['maxdepth']=max(st['maxdepth'],dep)
    mx=(p[0]+q[0])/2;my=(p[1]+q[1])/2;r=max(abs(q[0]-p[0]),abs(q[1]-p[1]))/2
    R=arb(0,float(r)*(1+1e-12)+1e-300);Z=acb(A(mx)+R,A(my)+R)
    m=abs((PI*Z).sin())/abs(PI*Z)*(-PI*abs(Z.imag)).exp()-psi_hi(max(my-r,Fr(0)))
    if m>0:
        st['minmargin']=min(st['minmargin'],float(m.mid())-float(m.rad()));return True
    if dep>=50:return False
    c=(mx,my);return chk(p,c,dep+1) and chk(c,q,dep+1)
t0=time.time();ok=all(chk(V[k],V[k+1],0) for k in range(K));b=1/(1-eta);q=b.numerator*10**9//b.denominator
print('arb verifier:',('PASS' if ok else 'FAIL'),'segments',K,st,'time %.1fs'%(time.time()-t0))
if ok:print('Certified: C_1d >= 1/(1-eta) = %s >= %d.%09d'%(b,q//10**9,q%10**9))
