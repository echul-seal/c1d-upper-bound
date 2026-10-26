import json,sys,time
from fractions import Fraction as Fr
from mpmath import iv
iv.prec=100
d=json.load(open(sys.argv[1]));eta=Fr(d['eta']);P=[(Fr(a),Fr(b)) for a,b in d['vertices']];n=len(P)-1
if not(0<eta<1 and P[0][1]==0==P[n][1] and min(y for _,y in P[1:n])>0 and 1<P[0][0]<2 and 0<P[n][0]<1):sys.exit('geometry FAIL')
for k in range(n):
    (a,b),(c,e)=P[k],P[k+1]
    if not (a-1)*e-b*(c-1)>0:sys.exit('not star-shaped at %d'%k)
I=lambda q:iv.mpf(q.numerator)/q.denominator
E=I(eta);pi=iv.pi
def good(p,q):
    u0,u1=sorted((p[0],q[0]));v0,v1=sorted((p[1],q[1]))
    U=iv.mpf([I(u0).a,I(u1).b]);W=iv.mpf([I(v0).a,I(v1).b])
    sh=(iv.exp(pi*W)-iv.exp(-pi*W))/2
    phi2=(iv.sin(pi*U)**2+sh**2)*iv.exp(-2*pi*W)/(pi**2*(U**2+W**2))
    if v0==0:ps=E
    else:
        t=I(v0);ps=(1-iv.exp(-pi*t*E))/(pi*t)
    return phi2.a>(ps**2).b
cnt=[0]
def rec(p,q,d):
    cnt[0]+=1
    if good(p,q):return True
    if d>=50:return False
    m=((p[0]+q[0])/2,(p[1]+q[1])/2);return rec(p,m,d+1) and rec(m,q,d+1)
t0=time.time();ok=all(rec(P[k],P[k+1],0) for k in range(n))
print('mpmath-iv verifier:',('PASS' if ok else 'FAIL'),'intervals',cnt[0],'time %.1fs'%(time.time()-t0))
if ok:
    b=1/(1-eta);print('=> C_1d >= %s ~ %.9f'%(b,float(b)))
