"""Analyze a finite-support JSON input. No inference from raw observations."""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from deadline_model import deadline_mask,mean_bounds,frontier
p=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('example.json')
d=json.loads(p.read_text());H=deadline_mask(d['division_times'],d['latencies'],d['deadline'])
a,b=d['division_mass'],d['latency_mass']
r={'data_status':d.get('data_status','unspecified; not validated'),'mean_bounds':mean_bounds(a,b,H)}
if 'mean' in d:
 f=frontier(a,b,H,d['mean']);r['frontier']=None if f is None else {k:v for k,v in f.items() if k not in ['z','e','row_zero']}
print(json.dumps(r,indent=2))
