import json,sys
from fractions import Fraction
V=json.load(open(sys.argv[1]))['v_int'];n=len(V);S=sum(V)
c=[0]*(2*n-1)
for i,a in enumerate(V):
    if a:
        for j,b in enumerate(V):c[i+j]+=a*b
M=max(abs(x) for x in c);r=Fraction(2*n*M,S*S)
print('n',n,'exact',r.numerator,'/',r.denominator)
print('floor 30 digits 0.'+str(r.numerator*10**30//r.denominator),'ceil15',-(-r.numerator*10**15//r.denominator))
