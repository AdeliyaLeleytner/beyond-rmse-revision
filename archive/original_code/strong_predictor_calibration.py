"""Outcome-triggered diagnostic: apply unchanged localizer to refitted global Plain/ADME predictors."""
import os,json
from pathlib import Path
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent
os.environ['REANALYSIS_OUT']=str(HERE/'repaired');os.environ['REANALYSIS_DATA']=str(HERE/'data_repair/audited_parent_matched.csv')
from analysis import OUT,load,PROTEINS,local_intervals,global_quantile,cluster_bootstrap,dump
df=load().set_index('ligand_id');out=OUT/'strong_predictor_calibration';out.mkdir(exist_ok=True);rows=[];comparisons=[]
for seed in range(42,47):
    rd=OUT/f'runs/butina04_seed{seed}';cfg=json.loads((rd/'complete.json').read_text());raw=np.column_stack([df[PROTEINS].mean(axis=1).to_numpy(),df[cfg['selected_markers']].to_numpy(float).clip(0,1)]);z=(raw-np.array(cfg['coordinate_mu']))/np.array(cfg['coordinate_sd']);loc=pd.DataFrame(z,index=df.index)
    for model in ['Plain','ADME']:
        d=OUT/f'global_models/seed{seed}'/model;p=pd.read_csv(d/'test_predictions.csv');cal=pd.read_csv(d/'calibration_predictions.csv');q=global_quantile(cal.abs_residual.to_numpy());width,ne=local_intervals(loc.loc[p.ligand_id].to_numpy(),loc.loc[cal.ligand_id].to_numpy(),cal.abs_residual.to_numpy(),q);p['width_joint']=width;p['covered_joint']=(p.abs_residual<=width/2).astype(int);p['score_joint']=width+20*np.maximum(p.abs_residual-width/2,0);p['neff_joint']=ne;p['score_global']=2*q+20*np.maximum(p.abs_residual-q,0)
        assert np.allclose(p['width_0.9'],2*q)
        for method in ['global','joint']:
            w=p['width_0.9'] if method=='global' else p.width_joint;c=p['covered_0.9'] if method=='global' else p.covered_joint;s=p['score_'+method];ci=cluster_bootstrap(p.cluster_id,np.column_stack([c,w,s]));rows.append({'seed':seed,'model':model,'method':method,'coverage':c.mean(),'mean_width':w.mean(),'mean_interval_score':s.mean(),'coverage_ci_low':ci[0,0],'coverage_ci_high':ci[1,0]})
        diff=p.score_joint-p.score_global;ci=cluster_bootstrap(p.cluster_id,diff);comparisons.append({'seed':seed,'model':model,'score_joint_minus_global':diff.mean(),'ci_low':ci[0,0],'ci_high':ci[1,0]});p.to_csv(out/f'{model}_seed{seed}_predictions.csv',index=False)
pd.DataFrame(rows).to_csv(out/'metrics.csv',index=False);pd.DataFrame(comparisons).to_csv(out/'paired_scores.csv',index=False)
print(pd.DataFrame(rows).groupby(['model','method'])[['coverage','mean_width','mean_interval_score']].mean().to_string())
