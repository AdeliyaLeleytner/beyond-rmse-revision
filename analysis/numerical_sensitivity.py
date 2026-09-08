from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
import json
import numpy as np,pandas as pd
from analysis import OUT,load,PROTEINS,FIXED,local_intervals,global_quantile

df=load();rd=OUT/'runs/original_reinstated';man=pd.read_csv(rd/'split_manifest.csv');tr=man.loc[man.role.eq('train'),'row_index'].to_numpy();ca=man.loc[man.role.eq('calibration'),'row_index'].to_numpy();te=man.loc[man.role.eq('test'),'row_index'].to_numpy();cal=pd.read_csv(rd/'calibration_predictions.csv');pred=pd.read_csv(rd/'test_predictions.csv')
raw=np.column_stack([np.nanmean(df[PROTEINS].to_numpy(float),axis=1),df[FIXED].to_numpy(float)]);res=cal.abs_residual.to_numpy();rows=[]
for original_q,scaler,label in [(True,ca,'published_quantile_and_scaling'),(False,ca,'exact_global_order_statistic_only'),(False,tr,'exact_global_quantile_and_training_scaler')]:
    q=global_quantile(res,original_q);mu=raw[scaler].mean(axis=0);sd=raw[scaler].std(axis=0);sd[sd==0]=1;z=(raw-mu)/sd
    w,ne=local_intervals(z[te],z[ca],res,q)
    rows.append({'variant':label,'global_width':2*q,'global_coverage':float((pred.abs_residual<=q).mean()),'local_width':float(w.mean()),'local_coverage':float((pred.abs_residual<=w/2).mean())})
pd.DataFrame(rows).to_csv(OUT/'numerical_sensitivity.csv',index=False)
# Populate actual sensitivity fallback rates; this does not change predictions/widths.
for rd in sorted((OUT/'runs').iterdir()):
    if not (rd/'complete.json').exists():continue
    cfg=json.loads((rd/'complete.json').read_text());p=pd.read_csv(rd/'test_predictions.csv');m=pd.read_csv(rd/'interval_metrics.csv')
    names=[s for s in m.method if s.startswith('joint_h')]
    if not names:continue
    man=pd.read_csv(rd/'split_manifest.csv');ca=man.loc[man.role.eq('calibration'),'row_index'].to_numpy();te=man.loc[man.role.eq('test'),'row_index'].to_numpy();cal=pd.read_csv(rd/'calibration_predictions.csv')
    raw=np.column_stack([np.nanmean(df[PROTEINS].to_numpy(float),axis=1),df[cfg['selected_markers']].to_numpy(float)]);z=(raw-np.array(cfg['coordinate_mu']))/np.array(cfg['coordinate_sd'])
    for name in names:
        w,ne=local_intervals(z[te],z[ca],cal.abs_residual.to_numpy(),cfg['global_quantile'],h=float(name.replace('joint_h','')))
        assert np.allclose(w,p['width_'+name])
        p['neff_'+name]=ne;m.loc[m.method.eq(name),'fallback_fraction']=float((ne<80).mean())
    p.to_csv(rd/'test_predictions.csv',index=False);m.to_csv(rd/'interval_metrics.csv',index=False)
print(pd.DataFrame(rows).to_string(index=False))
