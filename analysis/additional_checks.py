from paths import ROOT, PROCESSED, RESULTS, ADDITIONAL
"""Additional checks specified in ADDITIONAL_CHECKS_PROTOCOL.md; no model fitting."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import r2_score
from catboost import CatBoostRegressor
HERE=Path(__file__).resolve().parent
PREV=ROOT/'analysis'
sys.path.insert(0,str(PREV))
from analysis import PROTEINS,cluster_bootstrap,local_intervals
OUT=ADDITIONAL;OUT.mkdir(exist_ok=True)
def write(name,rows):pd.DataFrame(rows).to_csv(OUT/(name+'.csv'),index=False)
def ci(g,v):
    a=cluster_bootstrap(g,v);return {'ci_low':float(a[0,0]),'ci_high':float(a[1,0])}
def score(err,w):return w+20*np.maximum(err-w/2,0)
def kernel_unrestricted(zte,zca,res,h=.9):
    order=np.argsort(res);res=res[order];zca=zca[order];ws=[];ns=[]
    for start in range(0,len(zte),256):
        dist=((zte[start:start+256,None,:]-zca[None,:,:])**2).sum(2)
        w=np.exp(-.5*dist/h**2);sw=w.sum(1);ne=sw**2/np.maximum((w*w).sum(1),1e-300)
        ranks=(np.cumsum(w,axis=1)<.9*sw[:,None]).sum(1)
        ws.extend(2*res[np.minimum(ranks,len(res)-1)]);ns.extend(ne)
    return np.array(ws),np.array(ns)
def main():
    df=pd.read_csv(PROCESSED/'audited_parent_matched.csv.gz').set_index('ligand_id')
    docking=[];fallback=[];correlations=[];vifs=[];deciles=[];medians=[];counts=[];bbb=[];bbbpaired=[];permutation=[]
    proteins=df[PROTEINS].to_numpy(float)
    assert np.isfinite(proteins).all()
    lookup=pd.Series(np.arange(len(df)),index=df.index)
    for seed in range(42,47):
        rd=RESULTS/f'runs/butina04_seed{seed}';config=json.loads((rd/'complete.json').read_text())
        manifest=pd.read_csv(rd/'split_manifest.csv');p=pd.read_csv(rd/'test_predictions.csv');cal=pd.read_csv(rd/'calibration_predictions.csv')
        trids=manifest.loc[manifest.role.eq('train'),'ligand_id'];tr=lookup.loc[trids].to_numpy();ca=lookup.loc[cal.ligand_id].to_numpy();te=lookup.loc[p.ligand_id].to_numpy()
        bio=df[config['selected_markers']].clip(0,1).to_numpy(float)
        means=proteins[tr].mean(0);sds=proteins[tr].std(0);assert (sds>0).all()
        norm=((proteins-means)/sds).mean(1)
        pd.DataFrame({'target':PROTEINS,'training_mean':means,'training_sd':sds}).to_csv(OUT/f'docking_normalization_seed{seed}.csv',index=False)
        raw=np.column_stack([proteins.mean(1),bio]);z=(raw-raw[tr].mean(0))/raw[tr].std(0)
        raww,ne=kernel_unrestricted(z[te],z[ca],cal.abs_residual.to_numpy())
        q=config['global_quantile'];err=p.abs_residual.to_numpy();groups=p.cluster_id.to_numpy()
        for cutoff in [40,80,120]:
            w=np.where(ne<cutoff,2*q,raww)
            if cutoff==80:assert np.allclose(w,p.width_joint,atol=1e-12)
            fallback.append({'seed':seed,'minimum_neff':cutoff,'fallback_n':int((ne<cutoff).sum()),'n_test':len(p),'fallback_fraction':float((ne<cutoff).mean()),'coverage':float((err<=w/2).mean()),'mean_width':float(w.mean()),'interval_score':float(score(err,w).mean())})
        nr=np.column_stack([norm,bio]);nz=(nr-nr[tr].mean(0))/nr[tr].std(0)
        w,nne=local_intervals(nz[te],nz[ca],cal.abs_residual.to_numpy(),q)
        ns=score(err,w)
        for name,base in [('raw_mean44',p.interval_score_joint.to_numpy()),('phenotype_only',p.interval_score_pheno.to_numpy())]:
            docking.append({'seed':seed,'comparison':name,'normalized_score':float(ns.mean()),'reference_score':float(base.mean()),'delta':float((ns-base).mean()),**ci(groups,ns-base),'normalized_coverage':float((err<=w/2).mean()),'normalized_width':float(w.mean()),'normalized_fallback_fraction':float((nne<80).mean())})
        pd.DataFrame({'ligand_id':p.ligand_id,'normalized_mean44':norm[te],'width_normalized':w,'interval_score_normalized':ns}).to_csv(OUT/f'normalized_docking_predictions_seed{seed}.csv',index=False)
        for extra in [False,True]:
            names=config['selected_markers']+(['hERG-10um'] if extra else [])
            x=df.loc[trids,names].to_numpy(float);cx=pd.DataFrame(x,columns=names).corr(method='spearman')
            for i,a in enumerate(names):
                for j,b in enumerate(names):
                    if j>i:correlations.append({'seed':seed,'diagnostic_includes_unused_hERG10':extra,'a':a,'b':b,'spearman':cx.iloc[i,j]})
                other=np.delete(x,i,axis=1);design=np.column_stack([np.ones(len(x)),other]);pred=design@np.linalg.lstsq(design,x[:,i],rcond=None)[0];rr=r2_score(x[:,i],pred)
                vifs.append({'seed':seed,'diagnostic_includes_unused_hERG10':extra,'feature':a,'vif':1/(1-rr)})
        bt=bio.sum(1);edges=np.unique(np.quantile(bt[tr],np.linspace(0,1,11)));bins=np.searchsorted(edges[1:-1],bt[te],side='right')
        for b in sorted(set(bins)):
            ix=bins==b;res=(p.y-p.prediction).to_numpy()[ix]
            deciles.append({'seed':seed,'bin':int(b+1),'n':int(ix.sum()),'btox_mean':float(bt[te][ix].mean()),'residual_variance':float(res.var()),'mse':float((res*res).mean()),'mean_residual':float(res.mean())})
        # Whole-cluster bootstrap of proper-training median boundaries.
        trg=manifest.set_index('ligand_id').loc[trids,'cluster_id'].to_numpy();uniq=np.unique(trg);members=[np.flatnonzero(trg==g) for g in uniq];rng=np.random.default_rng(20260907)
        bs=[]
        for _ in range(1000):
            ix=np.concatenate([members[j] for j in rng.integers(0,len(uniq),len(uniq))]);bs.append([np.median(bt[tr][ix]),np.median(raw[tr,0][ix])])
        bs=np.array(bs);baseb=config['quadrant_thresholds']['phenotype'];basee=config['quadrant_thresholds']['mean_E']
        orig=(bt[te]>=baseb).astype(int)*2+(raw[te,0]<=basee)
        changed=np.zeros(len(te))
        for b,e in bs:changed+=((bt[te]>=b).astype(int)*2+(raw[te,0]<=e)!=orig)
        for j,label in enumerate(['BTox','mean44']):medians.append({'seed':seed,'coordinate':label,'median':[baseb,basee][j],'ci_low':float(np.quantile(bs[:,j],.025)),'ci_high':float(np.quantile(bs[:,j],.975)),'mean_quadrant_change_probability':float((changed/len(bs)).mean()),'fraction_change_probability_above_0.1':float((changed/len(bs)>.1).mean())})
        for quad,g in p.groupby('quadrant'):
            indicator=(g.width_joint<g.width_global).to_numpy(float);n=len(g);k=int(indicator.sum());ph=k/n;zz=1.95996398454;center=(ph+zz*zz/(2*n))/(1+zz*zz/n);half=zz*np.sqrt(ph*(1-ph)/n+zz*zz/(4*n*n))/(1+zz*zz/n)
            counts.append({'seed':seed,'quadrant':quad,'n':n,'n_narrower':k,'fraction':ph,'wilson_low':center-half,'wilson_high':center+half,**ci(g.cluster_id,indicator)})
        groot=RESULTS/f'global_models/seed{seed}'
        bp=pd.read_csv(groot/'BBB_pass/test_predictions.csv').set_index('ligand_id');yy=bp.y.to_numpy();groupsb=bp.cluster_id.to_numpy();errs={}
        for model in ['BBB_pass','Plain','Baseline']:
            pp=pd.read_csv(groot/model/'test_predictions.csv').set_index('ligand_id').loc[bp.index];assert np.allclose(pp.y,yy);yp=pp.prediction.to_numpy();res=yy-yp;errs[model]=res**2
            bbb.append({'seed':seed,'model':model,'n_test':len(bp),'target_variance':float(yy.var()),'rmse':float(np.sqrt(np.mean(res**2))),'rmse_over_target_sd':float(np.sqrt(np.mean(res**2)/yy.var())),'mae':float(np.abs(res).mean()),'r2':float(r2_score(yy,yp)),'spearman':float(spearmanr(yy,yp).statistic),'ccc':float(2*np.mean((yy-yy.mean())*(yp-yp.mean()))/(yy.var()+yp.var()+(yy.mean()-yp.mean())**2))})
        for ref in ['Plain','Baseline']:bbbpaired.append({'seed':seed,'comparison':'BBB_minus_'+ref,'mse_delta':float((errs['BBB_pass']-errs[ref]).mean()),**ci(groupsb,errs['BBB_pass']-errs[ref])})
        # Feature blocks are permuted jointly. No single-feature permutation sum.
        for modelname in ['Baseline','ADME']:
            md=groot/modelname;features=pd.read_csv(md/'features.csv').feature.tolist();med=pd.read_csv(md/'imputation_medians.csv',index_col=0).iloc[:,0];pp=pd.read_csv(md/'test_predictions.csv');x=df.loc[pp.ligand_id,features].apply(pd.to_numeric,errors='coerce').fillna(med);arr=x.to_numpy();model=CatBoostRegressor();model.load_model(str(md/'model.cbm'));yp=model.predict(x,thread_count=2);assert np.allclose(yp,pp.prediction,atol=1e-10)
            y=pp.y.to_numpy();base=(y-yp)**2;gr=np.array(['Fingerprints' if c.startswith('FP_') else 'PhysChem' if c in ['MW, g/mol','MW','logP'] else 'Proteins' if c in PROTEINS else 'ADME' for c in features]);rng=np.random.default_rng(20260908)
            for group in sorted(set(gr)):
                cols=np.flatnonzero(gr==group);deltas=[]
                for rep in range(5):
                    perm=rng.permutation(len(arr));a=arr.copy();a[:,cols]=arr[perm][:,cols];out=model.predict(a,thread_count=2);deltas.append((y-out)**2-base)
                delta=np.mean(deltas,axis=0);permutation.append({'seed':seed,'model':modelname,'group':group,'n_features':len(cols),'n_permutations':5,'mean_mse_increase':float(delta.mean()),'permutation_sd':float(np.std(np.mean(deltas,axis=1),ddof=1)),**ci(pp.cluster_id,delta)})
        print('ADDITIONAL CHECKS COMPLETE',seed,flush=True)
    for name,rows in [('docking_normalization',docking),('fallback_sensitivity',fallback),('phenotype_correlations',correlations),('phenotype_vif',vifs),('btox_residual_deciles',deciles),('quadrant_threshold_bootstrap',medians),('narrower_interval_confidence',counts),('bbb_same_subset',bbb),('bbb_paired_comparison',bbbpaired),('group_permutation_importance',permutation)]:write(name,rows)
    grid=pd.read_csv(RESULTS/'tables/figure5_grid.csv');maps={'grid_cells':len(grid),'mean_E_min':float(grid.mean_E.min()),'mean_E_max':float(grid.mean_E.max()),'BTox_min':float(grid.phenotype_burden.min()),'BTox_max':float(grid.phenotype_burden.max()),'mask_neff_below':120,'panels':{c:{'masked_n':int(grid[c].isna().sum()),'masked_fraction':float(grid[c].isna().mean()),'h':.75 if c=='coverage' else .45} for c in ['mean','variance','width','coverage']}}
    (OUT/'figure5_display_metadata.json').write_text(json.dumps(maps,indent=2)+'\n')
    (OUT/'verification.json').write_text(json.dumps({'five_seeds_complete':True,'raw_kernel_matches_saved_intervals':True,'saved_model_predictions_match':True,'BBB_test_ids_and_targets_identical':True,'all_target_docking_scores_finite':True,'normalization_fit':'proper training only','no_model_refitting':True},indent=2)+'\n')
if __name__=='__main__':main()
