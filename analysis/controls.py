from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
"""Controls specified without choosing parameters from their outcomes."""
from pathlib import Path
import json,itertools
import numpy as np,pandas as pd
from catboost import CatBoostRegressor
from sklearn.metrics import mean_squared_error
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from analysis import OUT,load,original_split,local_intervals,global_quantile,cluster_bootstrap,dump,PROTEINS,FIXED

def ablation(df):
    out=OUT/'analogue_ablation';out.mkdir(exist_ok=True)
    if (out/'summary.csv').exists():return
    split=original_split(df,True);tr=split['train'];te=split['test'];groups=df.Butina_clusters.to_numpy()
    retained=tr[~np.isin(groups[tr],np.unique(groups[te]))]
    assert not set(groups[retained])&set(groups[te])
    arms={'all_training':tr,'remove_test_cluster_members':retained}
    for seed in range(42,47):arms[f'random_size_matched_seed{seed}']=np.sort(np.random.default_rng(seed).choice(tr,len(retained),replace=False))
    cols=[f'FP_{i}' for i in range(2048)]+['MW','logP'];raw=df[cols].apply(pd.to_numeric,errors='coerce');y=df['-lgLD50, mol/kg'].to_numpy(float)
    preds=pd.DataFrame({'ligand_id':df.iloc[te].ligand_id.to_numpy(),'cluster_id':groups[te],'y':y[te],'original_train_cluster_present':np.isin(groups[te],groups[tr])})
    rows=[]
    for name,ix in arms.items():
        print('ABLATION',name,'train_n',len(ix),flush=True)
        med=raw.iloc[ix].median().fillna(0);x=raw.fillna(med)
        model=CatBoostRegressor(iterations=800,depth=6,learning_rate=.05,loss_function='RMSE',random_seed=42,verbose=False,allow_writing_files=False,thread_count=6)
        if name=='all_training':model.load_model(str(OUT/'runs'/('repaired_reference' if (OUT/'runs/repaired_reference/model.cbm').exists() else 'original_reproduction')/'model.cbm'))
        else:model.fit(x.iloc[ix],y[ix])
        model.save_model(str(out/(name+'.cbm')))
        pd.DataFrame({'ligand_id':df.iloc[ix].ligand_id.to_numpy(),'cluster_id':groups[ix]}).to_csv(out/(name+'_training.csv'),index=False)
        p=model.predict(x.iloc[te]);preds[name]=p
        rows.append({'arm':name,'n_train':len(ix),'rmse':float(np.sqrt(np.mean((p-y[te])**2))),'mae':float(np.mean(np.abs(p-y[te])))})
    preds.to_csv(out/'test_predictions.csv',index=False);pd.DataFrame(rows).to_csv(out/'summary.csv',index=False)
    # Bootstrap paired RMSE differences, conditioning on the fitted arms.
    e_remove=(preds.remove_test_cluster_members-preds.y).to_numpy()**2
    e_full=(preds.all_training-preds.y).to_numpy()**2
    e_random=np.column_stack([(preds[f'random_size_matched_seed{s}']-preds.y).to_numpy()**2 for s in range(42,47)])
    stats=[]
    for subset,mask in [('all_test',np.ones(len(preds),dtype=bool)),('test_clusters_represented_in_original_training',preds.original_train_cluster_present.to_numpy())]:
        g=groups[te][mask];u,inv=np.unique(g,return_inverse=True);count=np.bincount(inv)
        vals=np.column_stack([e_remove[mask],e_full[mask],e_random[mask]])
        sums=np.column_stack([np.bincount(inv,weights=vals[:,j]) for j in range(vals.shape[1])])
        rng=np.random.default_rng(20260907);boot=[]
        for start in range(0,2000,100):
            z=rng.integers(0,len(u),(100,len(u)));m=np.sqrt(sums[z].sum(axis=1)/count[z].sum(axis=1)[:,None])
            boot.append(np.column_stack([m[:,0]-m[:,1],m[:,0]-m[:,2:].mean(axis=1)]))
        ci=np.quantile(np.concatenate(boot),[.025,.975],axis=0);m=np.sqrt(vals.mean(axis=0))
        stats.append({'subset':subset,'n':int(mask.sum()),'rmse_remove':float(m[0]),'rmse_full':float(m[1]),'mean_rmse_random_removal':float(m[2:].mean()),'rmse_remove_minus_full':float(m[0]-m[1]),'full_difference_ci_low':float(ci[0,0]),'full_difference_ci_high':float(ci[1,0]),'rmse_remove_minus_random':float(m[0]-m[2:].mean()),'random_difference_ci_low':float(ci[0,1]),'random_difference_ci_high':float(ci[1,1])})
    pd.DataFrame(stats).to_csv(out/'paired_rmse_differences.csv',index=False)

def null_controls(df):
    out=OUT/'localization_controls';out.mkdir(exist_ok=True)
    rows=[];paired=[]
    for seed in range(42,47):
        name=f'butina04_seed{seed}';rd=OUT/'runs'/name
        manifest=pd.read_csv(rd/'split_manifest.csv');split={r:g.row_index.to_numpy() for r,g in manifest.groupby('role')}
        tr,ca,te=(split[k] for k in ['train','calibration','test'])
        pred=pd.read_csv(rd/'test_predictions.csv');cal=pd.read_csv(rd/'calibration_predictions.csv');cfg=json.loads((rd/'complete.json').read_text())
        model=CatBoostRegressor();model.load_model(str(rd/'model.cbm'))
        cols=[f'FP_{i}' for i in range(2048)]+['MW','logP'];x=df[cols];x=x.fillna(x.iloc[tr].median().fillna(0))
        trainpred=model.predict(x.iloc[tr]);q=cfg['global_quantile'];res=cal.abs_residual.to_numpy();err=pred.abs_residual.to_numpy()
        methods={}
        rawtr=trainpred[:,None];rawcal=cal.prediction.to_numpy()[:,None];rawtest=pred.prediction.to_numpy()[:,None]
        mu=rawtr.mean(axis=0);sd=rawtr.std(axis=0);sd[sd==0]=1
        methods['prediction_only']=local_intervals((rawtest-mu)/sd,(rawcal-mu)/sd,res,q)[0]
        phys=x[['MW','logP']].to_numpy(float);mu=phys[tr].mean(axis=0);sd=phys[tr].std(axis=0);sd[sd==0]=1
        methods['MW_logP']=local_intervals((phys[te]-mu)/sd,(phys[ca]-mu)/sd,res,q)[0]
        bio=np.column_stack([np.nanmean(df[PROTEINS].to_numpy(float),axis=1),df[cfg['selected_markers']].to_numpy(float)])
        z=(bio-np.array(cfg['coordinate_mu']))/np.array(cfg['coordinate_sd'])
        for ps in range(5):
            rng=np.random.default_rng(1000+ps)
            methods[f'permuted_biology_{ps}']=local_intervals(z[te][rng.permutation(len(te))],z[ca][rng.permutation(len(ca))],res,q)[0]
        for method,w in methods.items():
            cov=err<=w/2;score=w+20*np.maximum(err-w/2,0)
            ci=cluster_bootstrap(pred.cluster_id,np.column_stack([cov,w,score]))
            rows.append({'run':name,'method':method,'coverage':float(cov.mean()),'mean_width':float(w.mean()),'mean_interval_score':float(score.mean()),'score_ci_low':float(ci[0,2]),'score_ci_high':float(ci[1,2])})
            diff=pred.interval_score_joint.to_numpy()-score;dc=cluster_bootstrap(pred.cluster_id,diff)
            paired.append({'run':name,'comparator':method,'joint_minus_comparator_score':float(diff.mean()),'ci_low':float(dc[0,0]),'ci_high':float(dc[1,0])})
            pred['width_'+method]=w
        pred.to_csv(out/(name+'_predictions.csv'),index=False)
        print('NULL CONTROLS',name,flush=True)
    pd.DataFrame(rows).to_csv(out/'metrics.csv',index=False);pd.DataFrame(paired).to_csv(out/'paired_scores.csv',index=False)

def similarity(df):
    out=OUT/'similarity';out.mkdir(exist_ok=True)
    if (out/'summary.csv').exists():return
    d=df.sort_values('ligand_id');ids=d.ligand_id.to_numpy();bits=(OUT/'clusters/fingerprints.bin').read_bytes()
    assert len(bits)==256*len(ids)
    lookup={int(i):DataStructs.CreateFromBinaryText(bits[256*j:256*(j+1)]) for j,i in enumerate(ids)}
    rows=[]
    for rd in sorted((OUT/'runs').iterdir()):
        if not (rd/'complete.json').exists():continue
        m=pd.read_csv(rd/'split_manifest.csv');trainids=m.loc[m.role.eq('train'),'ligand_id'].to_numpy();fps=[lookup[int(i)] for i in trainids]
        values=[]
        for role in ['calibration','test']:
            query=m.loc[m.role.eq(role),'ligand_id'].to_numpy()
            vals=[]
            for i in query:
                a=np.asarray(DataStructs.BulkTanimotoSimilarity(lookup[int(i)],fps));best=int(a.argmax());vals.append(float(a[best]));values.append({'ligand_id':int(i),'role':role,'max_tanimoto_to_train':float(a[best]),'nearest_training_id':int(trainids[best])})
            v=np.array(vals);rows.append({'run':rd.name,'role':role,'n':len(v),'median':float(np.median(v)),'p90':float(np.quantile(v,.9)),'maximum':float(v.max()),'fraction_ge_0_6':float(np.mean(v>=.6)),'fraction_ge_0_8':float(np.mean(v>=.8)),'fraction_equal_1':float(np.mean(v==1))})
        pd.DataFrame(values).to_csv(out/(rd.name+'.csv'),index=False)
        print('SIMILARITY',rd.name,flush=True)
    pd.DataFrame(rows).to_csv(out/'summary.csv',index=False)

if __name__=='__main__':
    df=load();ablation(df);null_controls(df);similarity(df)
