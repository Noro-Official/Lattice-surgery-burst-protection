"""Second-referee count follow-up; protocol frozen before new trajectory generation.
Uses the recovered P1 law, explicit schedule, saved baseline event bank, all four decoders.
Runs valid low-grid points and all original typed/reference controls afresh; never selects source by outcome.
"""
from pathlib import Path
import sys, json, hashlib, platform, time, shutil, os
import numpy as np
import p1core as p
import campaign as cp
ORIGINAL=Path('/Users/erickin/Desktop/untitled folder/publication_package_20261031_revised/experiments/pacer_p1')
R=p.HERE;NAME='count_extension_20261106';STUDY=R/'results'/NAME

def prepare():
 STUDY.mkdir(parents=True,exist_ok=True)
 bank='events_d9_seed78001_n20000_300-180-120_fp1.npz'
 (R/'results/events').mkdir(exist_ok=True)
 if not (R/'results/events'/bank).exists():shutil.copy2(ORIGINAL/'results/events'/bank,R/'results/events'/bank)
 eligible=int(p.eligible(9).sum());G=p.cw.geometry(9);feasible=[];fails={};masks={}
 for n in range(3,241,3):
  try:
   w=p.mask_typed(9,(eligible-n)/eligible,5);actual=eligible-int(w.sum())
   if actual==n:feasible.append(n);masks[n]=w
  except ValueError as e:
   if n in [96,114]:fails[str(n)]=str(e)
 floor=min(feasible);newcounts=sorted(set([floor,floor+6,150,132,168,204]))
 conditions=[]
 for n in newcounts:
  w=masks[n];conditions.append({'name':f'n{n}|typed','mask_name':f'n{n}|typed','mask':np.nonzero(w)[0].tolist(),'tm':9})
 shallow=p.mask_blind(9,.3,int(round(6*9/13)));conditions.append({'name':'n168|shallow','mask_name':'n168|shallow','mask':np.nonzero(shallow)[0].tolist(),'tm':9})
 spec={'name':NAME,'d':9,'pool_seed':78001,'n_pool':20000,'alloc':[300,180,120],'conditions':conditions,'decoders':['D0C'],'shots':1024,'seed':78002,'block_events':4}
 contract={'date':'2026-11-06','status':'frozen before extension trajectories','source_selection':'Recovered archive with byte-identical manuscript summaries; no selection by effect size','requested_counts':[96,114,150],'infeasible_requested_counts':fails,'minimum_feasible_typed_count':floor,'tested_typed_counts':newcounts,'reference':'Shallow at 168 protected; regenerated on same bank and current runtime','schedule':p.schedule(9),'weak_XZ_counts':{str(n):{t:int((masks[n]&(G['role']==t+'-anc')).sum()) for t in ['X','Z']} for n in newcounts},'mask_hashes':{c['name']:p.mask_hash(cp._mk(9,c['mask'])) for c in conditions},'event_bank_sha256':hashlib.sha256((R/'results/events'/bank).read_bytes()).hexdigest(),'statistical_scope':'Post hoc grid refinement using the original 600-impact bank; all curves and reference regenerated, no pooled old/new trajectory estimates; nominal 4000-replicate event bootstrap with quiet uncertainty; no new independent confirmation','primary':'Count crossing for typed relative to shallow at168 under correlated background matching; correlated oracle and plain decoders secondary','analysis':'Piecewise log-risk crossing on tested grid; finite feasible-grid boundary explicitly reported; no automatic shrinking of width or change in type split','runtime':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__},'study':spec}
 f=STUDY/'FROZEN_PROTOCOL.json'
 if f.exists():
  old=json.loads(f.read_text());assert old.get('amendment',{}).get('decoders_after_amendment')==['D0C'];assert old['schedule']==contract['schedule'] and old['tested_typed_counts']==contract['tested_typed_counts'],'Protocol model/grid changed'
 else:f.write_text(json.dumps(contract,indent=2)+'\n')
 print('Feasibility',floor,'requested failures',fails,'grid',newcounts,flush=True)
 return spec

def bench(S):
 # One saved original condition, mode and event; a reproducibility diagnostic, not a selection criterion.
 cond=next(c for c in S['conditions'] if c['name']=='n168|typed').copy();cond['name']=cond['mask_name']='f0.3|typed_paper'
 src={**S,'name':'baselines','conditions':[cond]};ef=cp.events_for(src);jn,args=cp.make_jobs(src,ef,4,4)[0]
 before=time.perf_counter();out=cp.block(args);seconds=time.perf_counter()-before
 old=ORIGINAL/'results/baselines/blocks'/f'{jn}.json';stored=json.loads(old.read_text())
 same=[a['hist']==b['hist'] for a,b in zip(stored['rows'],out['rows'])]
 report={'job':jn,'seconds':seconds,'rows':len(out['rows']),'histograms_identical':all(same),'per_row_identical':same,'old_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),'new':out,'stored':stored,'scope':'One four-event original typed control in Z mode; no full-campaign reproduction claim'}
 (STUDY/'CONTROL_REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n');print('Control',seconds,'s','identical',same,flush=True)
 return seconds
if __name__=='__main__':
 S=prepare()
 if len(sys.argv)>1 and sys.argv[1]=='prepare':raise SystemExit()
 if len(sys.argv)>1 and sys.argv[1]=='bench':bench(S);raise SystemExit()
 from functools import partial
 import concurrent.futures, multiprocessing
 concurrent.futures.ProcessPoolExecutor=partial(concurrent.futures.ProcessPoolExecutor,mp_context=multiprocessing.get_context('fork'))
 # Inherit the warmed, immutable circuit/cache state; per-task seeds remain unchanged.
 raise SystemExit(cp.run(S,int(sys.argv[1]) if len(sys.argv)>1 else 6))
