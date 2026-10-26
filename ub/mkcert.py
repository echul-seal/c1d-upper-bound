import json,sys
d=json.load(open(sys.argv[1]));v=d['v'];s=max(abs(x) for x in v)
V=[int(round(x/s*10**12)) for x in v]
json.dump({'problem':'C_1d upper bound witness: signed step function, n equal cells on [-1/4,1/4], cell values proportional to v_int','n':len(V),'v_int':V},open(sys.argv[2],'w'))
