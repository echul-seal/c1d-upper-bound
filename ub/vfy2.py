import json,sys
from fractions import Fraction as Fr
V=[int(x) for x in json.load(open(sys.argv[1]))['v_int']];n=len(V);m=2*n-1
A=max(map(abs,V));b=(n*A*A).bit_length()+2;mask=(1<<b)-1
def pk(a):
    N=0
    for x in reversed(a):N=(N<<b)|x
    return N
def up(N):
    o=[]
    for _ in range(m):o.append(N&mask);N>>=b
    return o
P=[x if x>0 else 0 for x in V];Q=[-x if x<0 else 0 for x in V];p,q=pk(P),pk(Q)
c=[x+y-2*z for x,y,z in zip(up(p*p),up(q*q),up(p*q))]
S=sum(V);M=max(map(abs,c));r=Fr(2*n*M,S*S)
print('kronecker n',n,'exact',r.numerator,'/',r.denominator)
print('floor25 0.'+str(r.numerator*10**25//r.denominator),'nactive',sum(1 for x in c if abs(x)*10**6>=M*(10**6-1)))
