"""Frozen follow-up analysis: paired-event count crossing with shared quiet draws.
No extrapolation beyond the tested domain and no forced monotonic fit.
"""
from pathlib import Path
import json,hashlib,sys
import numpy as np
import campaign as cp
import p1core as p
NAME='count_extension_20261106';B=4000;SEED=20261106

def ci(x):return np.percentile(x,[2.5,97.5]).tolist()
def crossing(n,r,target):
 lr=np.log(np.maximum(r,1e-7));lt=np.log(max(target,1e-7))
 if lr[0]<=lt:return float(n[0]),'lower_edge',[int(n[0]),int(n[0])]
 for k in range(1,len(n)):
  if lr[k]<=lt:
   t=(lr[k-1]-lt)/(lr[k-1]-lr[k]);return float(n[k-1]+t*(n[k]-n[k-1])),'bracketed',[int(n[k-1]),int(n[k])]
 return float(n[-1]),'upper_edge',[int(n[-1]),int(n[-1])]

def plan():
 folder=cp.RES/NAME;folder.mkdir(exist_ok=True)
 qf=cp.RES/'quiet_d9.json';q=json.loads(qf.read_text())
 obj={'bootstrap_replicates':B,'seed':SEED,'quiet_file':'quiet_d9.json','quiet_sha256':hashlib.sha256(qf.read_bytes()).hexdigest(),'quiet_uncertainty':'Jeffreys beta draws shared across layouts and shared D0/F or D0C/FC pairs; independent of event and source-weight resampling','crossing':'First piecewise log-risk crossing, floor1e-7, no forced monotonicity, no extrapolation. All lower/upper-edge frequencies reported. Minimum feasible declared-rule count117. Infeasible96/114 are not simulated.','posthoc_scope':'Original600-impact bank; regenerated all typed points and shallow reference under current runtime; does not add independent impact confirmation'}
 f=folder/'ANALYSIS_PLAN.json'
 if f.exists():assert json.loads(f.read_text())==obj
 else:f.write_text(json.dumps(obj,indent=2)+'\n')
 return obj

def main():
 plan();D=cp.load(NAME);D['H']={k:v for k,v in D['H'].items() if k[1]=='D0C'};D['complete']={k:v for k,v in D['complete'].items() if k[1]=='D0C'};assert all(D['complete'].values());assert len(D['complete'])==7*2
 E=cp.Estimator(D,B=B,seed=SEED);qraw=json.loads((cp.RES/'quiet_d9.json').read_text());rng=np.random.default_rng(SEED+1)
 qb={}
 for base in ['D0','D0C']:
  for m in ['Z','X']:
   r=qraw[f'9|{m}|{base}'];qb[base,m]=rng.beta(r['fail']+.5,r['shots']-r['fail']+.5,size=B)
 def risk(cond,dec,boot=False):
  base='D0C' if dec in ['D0C','FC'] else 'D0';v=[]
  for m in ['Z','X']:
   q=qb[base,m] if boot else qraw[f'9|{m}|{base}']['p'];v.append(E.mean(cond,dec,m,boot)-q)
  return np.max(v,axis=0) if boot else float(max(v))
 counts=sorted(int(c['name'].split('|')[0][1:]) for c in D['S']['conditions'] if c['name'].endswith('|typed'))
 out={'n_impacts':D['n'],'shots':D['S']['shots'],'eligible':240,'protected_counts':counts,'minimum_feasible_typed_count':117,'reference_count':168,'decoders':{},'protocol_sha256':hashlib.sha256((cp.RES/NAME/'FROZEN_PROTOCOL.json').read_bytes()).hexdigest(),'scope':'New trajectories on the archived600-impact bank, nominal paired stratified bootstrap with quiet uncertainty; no new independent impact cohort'}
 for dec in ['D0C']:
  tgt=risk('n168|shallow',dec);tb=risk('n168|shallow',dec,True);vals=np.array([risk(f'n{n}|typed',dec) for n in counts]);boots=np.array([risk(f'n{n}|typed',dec,True) for n in counts]);star,kind,bracket=crossing(counts,vals,tgt)
  crosses=[crossing(counts,boots[:,i],tb[i]) for i in range(B)];stars=np.array([v[0] for v in crosses]);kinds=[v[1] for v in crosses]
  out['decoders'][dec]={'target':tgt,'target_ci':ci(tb),'typed_excess':vals.tolist(),'typed_excess_ci':[ci(v) for v in boots],'n_star':star,'point_status':kind,'crossing_bracket':bracket,'n_star_ci':ci(stars),'saving':168-star,'saving_ci':ci(168-stars),'bootstrap_lower_edge_fraction':kinds.count('lower_edge')/B,'bootstrap_upper_edge_fraction':kinds.count('upper_edge')/B,'point_nonmonotone_segments':int(np.sum(np.diff(vals)>0))}
  print(dec,json.dumps(out['decoders'][dec]),flush=True)
 (cp.RES/NAME/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='plan':plan()
 else:main()
