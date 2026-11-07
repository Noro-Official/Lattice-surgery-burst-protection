from pathlib import Path
import csv,gzip,json,math
import numpy as np
R=Path(__file__).resolve().parent
A=json.loads((R/'RESULTS.json').read_text());V=A['decoders']['D0C'];S=json.loads((R/'study.json').read_text());E=dict(np.load(R/S['events_file']));Q=json.loads((R/'quiet_d9.json').read_text());H={}
with gzip.open(R/'per_event_counts.csv.gz','rt') as f:
 for r in csv.DictReader(f):H.setdefault((r['condition'],r['input_mode']),np.empty(600))[int(r['event_index'])]=int(r['any_output_failures'])/int(r['shots'])
B=4000;seed=20261106;rng=np.random.default_rng(seed);counts=E['pool_counts'].astype(float);pi=counts/counts.sum();idx=[np.flatnonzero(E['stratum']==h) for h in range(3)];picks=[ii[rng.integers(0,len(ii),size=(B,len(ii)))] for ii in idx];weights=rng.dirichlet(counts+1,size=B)
rng=np.random.default_rng(seed+1);quiet={}
for base in ['D0','D0C']:
 for m in ['Z','X']:
  q=Q[f'9|{m}|{base}'];quiet[base,m]=rng.beta(q['fail']+.5,q['shots']-q['fail']+.5,size=B)
def risk(c,boot=False):
 v=[]
 for m in ['Z','X']:
  x=H[c,m];v.append(sum(weights[:,h]*x[picks[h]].mean(1) for h in range(3))-quiet['D0C',m] if boot else sum(pi[h]*x[idx[h]].mean() for h in range(3))-Q[f'9|{m}|D0C']['p'])
 return np.max(v,axis=0) if boot else max(v)
n=np.array(A['protected_counts']);y=np.array([risk(f'n{x}|typed') for x in n]);yb=np.array([risk(f'n{x}|typed',True) for x in n]);target=risk('n168|shallow');tb=risk('n168|shallow',True)
def cross(r,t):
 lr=np.log(np.maximum(r,1e-7));lt=np.log(max(t,1e-7))
 if lr[0]<=lt:return float(n[0])
 for k in range(1,len(n)):
  if lr[k]<=lt:return float(n[k-1]+(n[k]-n[k-1])*(lr[k-1]-lt)/(lr[k-1]-lr[k]))
 return float(n[-1])
stars=np.array([cross(yb[:,j],tb[j]) for j in range(B)]);star=cross(y,target);ci=np.percentile(stars,[2.5,97.5]);saving_ci=np.percentile(168-stars,[2.5,97.5])
assert np.allclose(y,V['typed_excess'],rtol=0,atol=1e-15) and np.allclose(ci,V['n_star_ci'],rtol=0,atol=1e-12) and np.allclose(saving_ci,V['saving_ci'],rtol=0,atol=1e-12)
assert abs(star-V['n_star'])<1e-12
print(json.dumps({'required_count':star,'count_interval':ci.tolist(),'saving':168-star,'saving_interval':saving_ci.tolist(),'lower_boundary_fraction':float(np.mean(stars==n[0])),'verified':True},indent=2))
