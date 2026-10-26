import json,numpy as np
from fractions import Fraction as Fr
d=json.load(open('lbs.json'));th=d['th'];fs=[]
def q(x):return Fr(round(x*10**6),10**6)
for r in sorted(d['res'],key=lambda r:-r['eta']):
    if r['minmarg']<=0.001:continue
    V=[[q(1+R*np.cos(a)),q(R*np.sin(a))] for a,R in zip(th,r['path'])];V[0][1]=Fr(0);V[-1][1]=Fr(0)
    fn='cP_%.3f.json'%r['eta'];json.dump({'eta':str(q(r['eta'])),'vertices':[[str(a),str(b)] for a,b in V]},open(fn,'w'));fs.append(fn)
print(' '.join(fs))
