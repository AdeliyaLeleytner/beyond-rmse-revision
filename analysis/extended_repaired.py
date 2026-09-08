from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
"""Identity sensitivity and global feature-set comparisons on repaired data."""
import os,json,ast,pickle,re
from pathlib import Path
import numpy as np,pandas as pd
from catboost import CatBoostRegressor,Pool
from sklearn.decomposition import PCA
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import mean_squared_error,mean_absolute_error,r2_score
from scipy.stats import spearmanr
HERE=Path(__file__).resolve().parent
os.environ.setdefault('REANALYSIS_DATA',str(PROCESSED/'audited_parent_matched.csv.gz'))
os.environ.setdefault('REANALYSIS_OUT',str(RESULTS))
import analysis as a

def strict_sensitivity():
    df=pd.read_csv(PROCESSED/'strict_exact_identity.csv.gz');a.OUT=RESULTS/'strict_identity';a.OUT.mkdir(exist_ok=True)
    ids=set(df.ligand_id)
    for seed in range(42,47):
        m=pd.read_csv(RESULTS/f'runs/butina04_seed{seed}/split_manifest.csv');m=m[m.ligand_id.isin(ids)]
        groups=m.set_index('ligand_id').loc[df.ligand_id,'cluster_id'].to_numpy()
        split={role:a.get_rows(df,g.ligand_id.to_numpy()) for role,g in m.groupby('role')}
        a.run(f'butina04_seed{seed}',df,split,groups)
    a.OUT=RESULTS

def global_models(df_override=None,model_names=None):
    df=a.load() if df_override is None else df_override;out=a.OUT/'global_models';out.mkdir(exist_ok=True)
    params=json.loads((a.REPO/'results/models/experiment_config.json').read_text())['fixed_params']
    if model_names is not None:params={k:v for k,v in params.items() if k in model_names}
    adme=json.loads((SOURCE/'admet_features.json').read_text())
    fp=sorted([c for c in df if c.startswith('FP_')],key=lambda c:int(c.split('_')[1]));core=['MW, g/mol','logP']+fp
    y=df[a.TARGET].to_numpy(float);results=[]
    for seed in range(42,47):
        manifest=pd.read_csv(a.OUT/f'runs/butina04_seed{seed}/split_manifest.csv');groups=manifest.set_index('ligand_id').loc[df.ligand_id,'cluster_id'].to_numpy()
        split={role:a.get_rows(df,g.ligand_id.to_numpy()) for role,g in manifest.groupby('role')}
        for name in params:
            rd=out/f'seed{seed}'/name.replace(' ','_');rd.mkdir(parents=True,exist_ok=True)
            if (rd/'complete.json').exists():results.append(json.loads((rd/'complete.json').read_text()));continue
            ix={k:v[df.iloc[v].bbb_rule_pass.eq(1).to_numpy()] if name=='BBB pass' else v for k,v in split.items()}
            tr,ca,te=(ix[k] for k in ['train','calibration','test'])
            fit0,val0=next(GroupShuffleSplit(n_splits=1,test_size=.1,random_state=42).split(tr,groups=groups[tr]));fit,val=tr[fit0],tr[val0]
            cols=core+([] if name=='Plain' else a.PROTEINS)+(adme if name=='ADME' else [])
            raw=df[cols].apply(pd.to_numeric,errors='coerce')
            # Transform parameters fitted solely on the inner fitting part; no ES/cal/test fitting.
            med=raw.iloc[fit].median().fillna(0);x=raw.fillna(med);med.to_csv(rd/'imputation_medians.csv')
            if name=='PCA':
                pca=PCA(n_components=3,random_state=42).fit(x.iloc[fit][a.PROTEINS]);pc=pca.transform(x[a.PROTEINS]);x=x[core].copy()
                for j in range(3):x[f'PC{j+1}']=pc[:,j]
                with (rd/'pca.pkl').open('wb') as f:pickle.dump(pca,f)
            kw=dict(params[name],loss_function='RMSE',random_seed=42,verbose=False,allow_writing_files=False,thread_count=6)
            preliminary=CatBoostRegressor(**kw);preliminary.fit(x.iloc[fit],y[fit],eval_set=(x.iloc[val],y[val]),early_stopping_rounds=200,use_best_model=True)
            kw['iterations']=preliminary.tree_count_;model=CatBoostRegressor(**kw);model.fit(x.iloc[tr],y[tr]);model.save_model(str(rd/'model.cbm'))
            pd.DataFrame({'ligand_id':df.iloc[tr].ligand_id.to_numpy(),'inner_role':np.where(np.isin(tr,val),'early_stopping','fit')}).to_csv(rd/'inner_manifest.csv',index=False)
            pd.DataFrame({'feature':x.columns}).to_csv(rd/'features.csv',index=False)
            pc=model.predict(x.iloc[ca]);pt=model.predict(x.iloc[te]);res=np.abs(y[ca]-pc);err=np.abs(y[te]-pt)
            p=pd.DataFrame({'ligand_id':df.iloc[te].ligand_id.to_numpy(),'cluster_id':groups[te],'y':y[te],'prediction':pt,'abs_residual':err})
            pd.DataFrame({'ligand_id':df.iloc[ca].ligand_id.to_numpy(),'y':y[ca],'prediction':pc,'abs_residual':res}).to_csv(rd/'calibration_predictions.csv',index=False)
            rec={'seed':seed,'model':name,'n_train':len(tr),'n_calibration':len(ca),'n_test':len(te),'trees':model.tree_count_,'rmse':float(np.sqrt(mean_squared_error(y[te],pt))),'mae':float(mean_absolute_error(y[te],pt)),'r2':float(r2_score(y[te],pt)),'spearman':float(spearmanr(y[te],pt).statistic)}
            ci=a.cluster_bootstrap(groups[te],err**2);rec.update(rmse_ci_low=float(np.sqrt(ci[0,0])),rmse_ci_high=float(np.sqrt(ci[1,0])))
            for level in [.8,.9,.95]:
                k=int(np.ceil((len(res)+1)*level));q=float(np.partition(res,k-1)[k-1]);cov=err<=q
                rec[f'coverage_{level}']=float(cov.mean());rec[f'width_{level}']=2*q
                p[f'width_{level}']=2*q;p[f'covered_{level}']=cov.astype(int)
            rec['interval_score_90']=float((p['width_0.9']+20*np.maximum(err-p['width_0.9']/2,0)).mean());p.to_csv(rd/'test_predictions.csv',index=False)
            if seed==42:
                rng=np.random.default_rng(20260907);sample=np.sort(rng.choice(len(te),min(512,len(te)),replace=False));shap=model.get_feature_importance(Pool(x.iloc[te[sample]],y[te[sample]]),type='ShapValues',thread_count=6)[:,:-1]
                f=pd.DataFrame({'feature':x.columns,'mean_absolute_shap':np.abs(shap).mean(axis=0)})
                f['group']=['fingerprints' if c.startswith('FP_') else 'docking' if c in a.PROTEINS or c.startswith('PC') else 'ADME' if c in adme else 'physicochemical' for c in f.feature]
                f.sort_values('mean_absolute_shap',ascending=False).to_csv(rd/'shap_features.csv',index=False);f.groupby('group').mean_absolute_shap.sum().to_csv(rd/'shap_groups.csv');p.iloc[sample][['ligand_id']].to_csv(rd/'shap_sample.csv',index=False)
            a.dump(rd/'complete.json',rec);results.append(rec);print('GLOBAL',seed,name,'RMSE',rec['rmse'],flush=True)
    pd.DataFrame(results).to_csv(out/'metrics.csv',index=False)
    pd.DataFrame(results).groupby('model')[['rmse','mae','r2','spearman','coverage_0.9','width_0.9','interval_score_90']].agg(['mean','std','min','max']).to_csv(out/'five_split_summary.csv')
    if 'Baseline' not in params:return
    # Paired errors on shared test IDs; BBB is compared with baseline on BBB subset only.
    comparisons=[]
    for seed in range(42,47):
        base=pd.read_csv(out/f'seed{seed}/Baseline/test_predictions.csv').set_index('ligand_id')
        for name in ['PCA','ADME','Plain','BBB_pass']:
            p=pd.read_csv(out/f'seed{seed}'/name/'test_predictions.csv').set_index('ligand_id');b=base.loc[p.index]
            value=p.abs_residual.to_numpy()**2-b.abs_residual.to_numpy()**2;ci=a.cluster_bootstrap(p.cluster_id,value)
            comparisons.append({'seed':seed,'model':name,'n':len(p),'mse_minus_baseline':float(value.mean()),'ci_low':float(ci[0,0]),'ci_high':float(ci[1,0])})
    pd.DataFrame(comparisons).to_csv(out/'paired_comparisons.csv',index=False)

if __name__=='__main__':
    strict_sensitivity();global_models()
