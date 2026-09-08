"""Reproducible retrospective Reviewer 1 reanalysis. Original sources are read-only."""
from pathlib import Path
import os, argparse, hashlib, itertools, json, math, platform, time, warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from catboost import CatBoostRegressor

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('REANALYSIS_OUT',Path(__file__).resolve().parent))
REPO=ROOT/'tmp/reviewer1_audit/Beyond-RMSE'
DATA=Path(os.environ.get('REANALYSIS_DATA',REPO/'data/df_final.csv'))
TARGET='-lgLD50, mol/kg'
FIXED=['hERG','Respiratory','Neurotoxicity-DI']
CANDIDATES=['hERG','Respiratory','Neurotoxicity-DI','hERG-10um','EI','DILI','Hematotoxicity','EC','Carcinogenicity','Ototoxicity','Ames','NonGenotoxic_Carcinogenicity','Genotoxicity','H-HT','Nephrotoxicity-DI','SkinSen','Acute_Aquatic_Toxicity','Genotoxic_Carcinogenicity_Mutagenicity','Skin_Sensitization']
PROTEINS=json.loads((REPO/'results/tables/conformal_run_config.json').read_text())['proteins']

def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def load():
    d=pd.read_csv(DATA)
    assert d.ligand_id.is_unique and d[TARGET].notna().all()
    return d

def get_rows(df,ids):
    pos=pd.Series(np.arange(len(df)),index=df.ligand_id)
    assert set(ids)<=set(pos.index)
    return pos.loc[ids].to_numpy()

def original_split(df,include_excluded=False):
    r=pd.read_csv(REPO/'results/splits/split_indices.csv')
    r=r[r.split_index.eq(0)]
    ids=r.loc[r['set'].eq('train'),'ligand_id'].to_numpy().copy()
    np.random.default_rng(42).shuffle(ids)
    n=int(round(len(ids)*.2))
    te=r.loc[r['set'].eq('test'),'ligand_id'].to_numpy()
    if not include_excluded:te=te[te!=1298]
    available=set(df.ligand_id)
    return {k:get_rows(df,np.array([i for i in v if i in available])) for k,v in {'train':np.sort(ids[n:]),'calibration':np.sort(ids[:n]),'test':te}.items()}

def registry_split(df):
    r=pd.read_csv(REPO/'results/splits/split_indices.csv');r=r[r.split_index.eq(1)]
    return {('calibration' if k=='val' else k):get_rows(df,g.ligand_id.to_numpy()) for k,g in r.groupby('set')}

def group_split(df,groups,seed):
    development,test=next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=seed).split(df,groups=groups))
    tr,cal=next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=seed).split(development,groups=groups[development]))
    return {'train':development[tr],'calibration':development[cal],'test':test}

def audit(df,split,groups):
    ans={'counts':{k:len(v) for k,v in split.items()},'clusters':{k:len(np.unique(groups[v])) for k,v in split.items()},'pairs':{}}
    for a,b in itertools.combinations(split,2):
        assert not set(split[a])&set(split[b])
        common=np.intersect1d(groups[split[a]],groups[split[b]])
        ans['pairs'][a+'/'+b]={'shared_clusters':len(common),'a_molecules':int(np.isin(groups[split[a]],common).sum()),'b_molecules':int(np.isin(groups[split[b]],common).sum())}
    assigned=np.concatenate(list(split.values()))
    ans['unassigned_ids']=df.loc[~df.index.isin(assigned),'ligand_id'].astype(int).tolist()
    return ans

def screen(df,train,groups,out):
    # All selection labels are confined to proper training, never calibration/test.
    y=df.iloc[train][TARGET].to_numpy(float)
    vals=df.iloc[train][CANDIDATES].apply(pd.to_numeric,errors='coerce').to_numpy(float)
    rows=[]
    inner=list(GroupShuffleSplit(n_splits=30,test_size=.2,random_state=20260907).split(vals,groups=groups[train]))
    for j,name in enumerate(CANDIDATES):
        a=np.array([spearmanr(vals[t,j],y[t]).statistic for t,v in inner])
        b=np.array([spearmanr(vals[v,j],y[v]).statistic for t,v in inner])
        rows.append({'marker':name,'rho_training_mean':float(np.mean(a)),'rho_training_sd':float(np.std(a,ddof=1)),'rho_internal_validation_mean':float(np.mean(b)),'rho_internal_validation_sd':float(np.std(b,ddof=1)),'positive_training_fraction':float(np.mean(a>0))})
    table=pd.DataFrame(rows).sort_values('rho_training_mean',ascending=False)
    chosen=[];axes=set()
    for row in table.itertuples():
        axis='hERG' if row.marker in ['hERG','hERG-10um'] else row.marker
        if row.rho_training_mean>0 and axis not in axes:
            chosen.append(row.marker);axes.add(axis)
        if len(chosen)==3:break
    assert len(chosen)==3
    table['selected']=table.marker.isin(chosen)
    table.to_csv(out/'marker_screening.csv',index=False)
    return chosen

def global_quantile(residual,original=False):
    n=len(residual);k=int(math.ceil((n+1)*.9))
    if k>n:return float('inf')
    return float(np.quantile(residual,k/n,method='higher') if original else np.partition(residual,k-1)[k-1])

def local_intervals(query,cal,resid,global_q,h=.9):
    order=np.argsort(resid);rv=resid[order];cal=cal[order]
    widths=[];neffs=[]
    for start in range(0,len(query),256):
        dist=((query[start:start+256,None,:]-cal[None,:,:])**2).sum(axis=2)
        w=np.exp(-.5*dist/h**2);sw=w.sum(axis=1)
        ne=sw**2/np.maximum((w*w).sum(axis=1),1e-300)
        ranks=(np.cumsum(w,axis=1)<.9*sw[:,None]).sum(axis=1)
        q=rv[np.minimum(ranks,len(rv)-1)]
        q=np.where((ne<80)|(sw<=0),global_q,q)
        widths.extend(2*q);neffs.extend(ne)
    return np.array(widths),np.array(neffs)

def cluster_bootstrap(groups,values,n=2000):
    # Resample test clusters, preserving all rows and paired methods in each draw.
    unique,inv=np.unique(groups,return_inverse=True);counts=np.bincount(inv)
    vals=np.asarray(values);vals=vals[:,None] if vals.ndim==1 else vals
    sums=np.column_stack([np.bincount(inv,weights=vals[:,j]) for j in range(vals.shape[1])])
    rng=np.random.default_rng(20260907);stats=[]
    for start in range(0,n,100):
        ix=rng.integers(0,len(unique),size=(min(100,n-start),len(unique)))
        stats.append(sums[ix].sum(axis=1)/counts[ix].sum(axis=1)[:,None])
    return np.quantile(np.concatenate(stats),[.025,.975],axis=0)

def run(name,df,split,groups,original=False,thread_count=6,strict=True,reuse=None,screen_markers=True):
    out=OUT/'runs'/name;out.mkdir(parents=True,exist_ok=True)
    if (out/'complete.json').exists():
        print('CACHED',name,flush=True);return json.loads((out/'complete.json').read_text())
    start=time.time();checks=audit(df,split,groups)
    if strict:assert all(v['shared_clusters']==0 for v in checks['pairs'].values())
    dump(out/'split_audit.json',checks)
    manifest=pd.concat([pd.DataFrame({'ligand_id':df.iloc[ix].ligand_id.to_numpy(),'role':role,'cluster_id':groups[ix],'row_index':ix}) for role,ix in split.items()])
    manifest.to_csv(out/'split_manifest.csv',index=False)
    tr,ca,te=(split[x] for x in ['train','calibration','test'])
    cols=sorted([c for c in df if c.startswith('FP_')],key=lambda c:int(c.split('_')[1]))+['MW','logP']
    x=df[cols].apply(pd.to_numeric,errors='coerce');med=x.iloc[tr].median().fillna(0);x=x.fillna(med)
    y=df[TARGET].to_numpy(float)
    model=CatBoostRegressor(iterations=800,depth=6,learning_rate=.05,loss_function='RMSE',random_seed=42,verbose=False,allow_writing_files=False,thread_count=thread_count)
    if reuse:model.load_model(str(reuse))
    else:model.fit(x.iloc[tr],y[tr])
    model.save_model(str(out/'model.cbm'));med.to_csv(out/'training_medians.csv')
    pc=model.predict(x.iloc[ca]);pt=model.predict(x.iloc[te]);res=np.abs(y[ca]-pc)
    pd.DataFrame({'ligand_id':df.iloc[ca].ligand_id.to_numpy(),'y':y[ca],'prediction':pc,'abs_residual':res}).to_csv(out/'calibration_predictions.csv',index=False)
    q=global_quantile(res,original);selected=screen(df,tr,groups,out) if screen_markers else FIXED
    mean_e=np.nanmean(df[PROTEINS].to_numpy(float),axis=1)
    assert np.isfinite(mean_e).all()
    bio=df[selected].apply(pd.to_numeric,errors='coerce').fillna(0).clip(0,1).to_numpy(float)
    raw=np.column_stack([mean_e,bio]);scaler=ca if original else tr
    mu=raw[scaler].mean(axis=0);sd=raw[scaler].std(axis=0);sd[sd==0]=1
    z=(raw-mu)/sd
    pred=pd.DataFrame({'ligand_id':df.iloc[te].ligand_id.to_numpy(),'cluster_id':groups[te],'y':y[te],'prediction':pt,'abs_residual':np.abs(y[te]-pt),'mean_E_44':mean_e[te],'phenotype_burden':bio[te].sum(axis=1)})
    for j,m in enumerate(selected):pred[m]=bio[te,j]
    methods={'global':np.full(len(te),2*q)}
    for method,features in [('joint',[0,1,2,3]),('pheno',[1,2,3]),('mie',[0])]:
        w,ne=local_intervals(z[te][:,features],z[ca][:,features],res,q)
        methods[method]=w;pred['neff_'+method]=ne
    # Bandwidth sensitivity is specified before outcomes; never choose by test coverage.
    if name.endswith('seed42') or name=='legacy_groups_registry1':
        for h in [.6,1.2]:
            w,ne=local_intervals(z[te],z[ca],res,q,h=h);methods[f'joint_h{h}']=w;pred[f'neff_joint_h{h}']=ne
    rows=[]
    err=np.abs(y[te]-pt)
    for method,w in methods.items():
        cov=err<=w/2;score=w+20*np.maximum(err-w/2,0)
        pred['width_'+method]=w;pred['covered_'+method]=cov.astype(int);pred['interval_score_'+method]=score
        ci=cluster_bootstrap(groups[te],np.column_stack([cov,w,score]))
        rows.append({'method':method,'coverage':float(cov.mean()),'coverage_ci_low':float(ci[0,0]),'coverage_ci_high':float(ci[1,0]),'mean_width':float(w.mean()),'width_ci_low':float(ci[0,1]),'width_ci_high':float(ci[1,1]),'median_width':float(np.median(w)),'width_p10':float(np.quantile(w,.1)),'width_p90':float(np.quantile(w,.9)),'mean_interval_score':float(score.mean()),'interval_score_ci_low':float(ci[0,2]),'interval_score_ci_high':float(ci[1,2]),'fallback_fraction':float((pred['neff_'+method]<80).mean()) if method!='global' else 0.0})
    pd.DataFrame(rows).to_csv(out/'interval_metrics.csv',index=False)
    threshold_rows=np.concatenate(list(split.values())) if original else tr
    bth=float(np.median(bio[threshold_rows].sum(axis=1)));eth=float(np.median(mean_e[threshold_rows]))
    pred['quadrant']=np.char.add(np.char.add(np.where(pred.phenotype_burden>=bth,'Pheno+','Pheno-'),' '),np.where(pred.mean_E_44<=eth,'MIE+','MIE-'))
    quad=[]
    for quadname,g in pred.groupby('quadrant'):
        ci=cluster_bootstrap(g.cluster_id,np.column_stack([g.covered_joint,g.width_joint]))
        quad.append({'quadrant':quadname,'n':len(g),'clusters':int(g.cluster_id.nunique()),'joint_width':float(g.width_joint.mean()),'pheno_width':float(g.width_pheno.mean()),'mie_width':float(g.width_mie.mean()),'coverage':float(g.covered_joint.mean()),'coverage_ci_low':float(ci[0,0]),'coverage_ci_high':float(ci[1,0]),'interval_lt_global_n':int((g.width_joint<g.width_global).sum())})
    pd.DataFrame(quad).to_csv(out/'quadrants.csv',index=False);pred.to_csv(out/'test_predictions.csv',index=False)
    delta=np.column_stack([pred.covered_joint-pred.covered_global,pred.width_joint-pred.width_global,pred.interval_score_joint-pred.interval_score_global])
    ci=cluster_bootstrap(groups[te],delta)
    paired={key:{'mean':float(delta[:,j].mean()),'ci_low':float(ci[0,j]),'ci_high':float(ci[1,j])} for j,key in enumerate(['coverage_joint_minus_global','width_joint_minus_global','interval_score_joint_minus_global'])}
    dump(out/'paired_comparison.json',paired)
    rmse=float(np.sqrt(mean_squared_error(y[te],pt)))
    base={'run':name,'counts':checks['counts'],'selected_markers':selected,'rmse':rmse,'mae':float(mean_absolute_error(y[te],pt)),'r2':float(r2_score(y[te],pt)),'spearman':float(spearmanr(y[te],pt).statistic),'global_quantile':q,'original_quantile_rule':original,'coordinate_scaler_fit':'calibration' if original else 'proper_train','coordinate_mu':mu.tolist(),'coordinate_sd':sd.tolist(),'quadrant_thresholds':{'phenotype':bth,'mean_E':eth},'model_seed':42,'elapsed_seconds':time.time()-start,'train_fingerprint_count':len(cols)-2}
    if original and name=='original_reproduction':
        old=pd.read_csv(REPO/'results/tables/conformal_test_predictions.csv').set_index('ligand_id')
        a=pred.set_index('ligand_id').loc[old.index]
        base['reproduction']={'max_prediction_difference':float(np.max(np.abs(a.prediction-old.pred_lgld50_plain))),'max_joint_width_difference':float(np.max(np.abs(a.width_joint-old.local_interval_width_90))),'covered_count_difference':int(a.covered_joint.sum()-old.covered_90.sum())}
    dump(out/'complete.json',base)
    print('COMPLETE',name,'RMSE',round(rmse,6),'global/joint coverage',pred.covered_global.mean(),pred.covered_joint.mean(),'width',pred.width_global.mean(),pred.width_joint.mean(),'seconds',round(time.time()-start),flush=True)
    return base

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['original']);args=ap.parse_args()
    df=load();groups=df.Butina_clusters.to_numpy()
    saved=OUT/'runs/original_reproduction/model.cbm'
    run('original_reproduction',df,original_split(df),groups,original=True,strict=False,screen_markers=False,reuse=saved if saved.exists() else None)
    run('original_reinstated',df,original_split(df,True),groups,original=True,strict=False,reuse=OUT/'runs/original_reproduction/model.cbm',screen_markers=False)
