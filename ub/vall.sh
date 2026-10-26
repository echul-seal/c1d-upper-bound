#!/bin/sh
t=$1
python3 mkcert.py $t.json c_$t.json
python3 vfy2.py c_$t.json > vf_$t.txt
python3 -c "import json;d=json.load(open('c_'+'$t'+'.json'));json.dump({'values':d['v_int']},open('cons_'+'$t'+'.json','w'))"
python3 -I v224.py cons_$t.json > v224_$t.txt
python3 verify_int.py c_$t.json > vi_$t.txt
sha256sum c_$t.json cons_$t.json >> vf_$t.txt
