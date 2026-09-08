import os,json,itertools,hashlib
from pathlib import Path
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent
os.environ['REANALYSIS_DATA']=str(HERE/'data_repair/audited_parent_matched.csv');os.environ['REANALYSIS_OUT']=str(HERE/'repaired')
from analysis import OUT,DATA,load,screen,global_quantile,local_intervals,dump,TARGET

def summarize(folder):
    rows=[]
    for f in sorted((folder/'runs').glob('*/complete.json')):
        cfg=json.loads(f.read_text());m=pd.read_csv(f.parent/'interval_metrics.csv').set_index('method')
        row={k:cfg[k] for k in ['run','rmse','mae','r2','spearman']};row.update(cfg['counts'])
        for method in ['global','joint']:
            for short,col in [('coverage','coverage'),('width','mean_width'),('score','mean_interval_score')]:row[method+'_'+short]=m.loc[method,col]
        rows.append(row)
    pd.DataFrame(rows).to_csv(folder/'overview.csv',index=False)

def verify():
    df=load();checks=[]
    assert global_quantile(np.arange(100))==90
    w,ne=local_intervals(np.zeros((1,1)),np.zeros((100,1)),np.arange(100),100);assert w[0]==178 and ne[0]==100
    identity=pd.read_csv(HERE/'data_repair/identity_manifest.csv');assert len(df)==12584
    allowed=identity[identity.identity_status.isin(['exact_isomeric_identity','unique_charge_parent'])];assert set(allowed.ligand_id)==set(df.ligand_id)
    assert 1298 in set(df.ligand_id) and 9116 in set(df.ligand_id) and 11434 not in set(df.ligand_id)
    assert allowed.exact_identity.sum()==12193
    assert allowed[~allowed.exact_identity].same_uncharged_isomeric_parent.all() and allowed[~allowed.exact_identity].unique_parent.all()
    for folder,n in [(OUT,12584),(OUT/'strict_identity',12193)]:
        found=0
        for rd in sorted((folder/'runs').iterdir()):
            if not (rd/'complete.json').exists():continue
            found+=1;cfg=json.loads((rd/'complete.json').read_text());m=pd.read_csv(rd/'split_manifest.csv');p=pd.read_csv(rd/'test_predictions.csv');ca=pd.read_csv(rd/'calibration_predictions.csv')
            assert m.ligand_id.is_unique and len(m)==(n-1 if rd.name=='repaired_reference' else n)
            roles={r:set(g.ligand_id) for r,g in m.groupby('role')}
            for x,y in itertools.combinations(roles,2):assert not roles[x]&roles[y]
            assert set(p.ligand_id)==roles['test'] and set(ca.ligand_id)==roles['calibration']
            if rd.name.startswith('butina'):
                groups={r:set(g.cluster_id) for r,g in m.groupby('role')}
                for x,y in itertools.combinations(groups,2):assert not groups[x]&groups[y]
            assert np.isclose(cfg['rmse'],np.sqrt(((p.y-p.prediction)**2).mean()))
            for row in pd.read_csv(rd/'interval_metrics.csv').itertuples():
                err=np.abs(p.y-p.prediction);width=p['width_'+row.method]
                assert np.array_equal(p['covered_'+row.method].to_numpy(),(err<=width/2).astype(int).to_numpy())
                assert np.isclose(row.mean_width,width.mean()) and np.isclose(row.coverage,(err<=width/2).mean())
                assert np.isclose(row.mean_interval_score,(width+20*np.maximum(err-width/2,0)).mean())
            checks.append({'cohort':folder.name,'run':rd.name,'n':len(m),'verified':True})
        assert found==(13 if folder==OUT else 5),found
    # Holdout-label perturbation checks actual marker-screen implementation.
    rd=OUT/'runs/butina04_seed42';m=pd.read_csv(rd/'split_manifest.csv');train=m.loc[m.role.eq('train'),'row_index'].to_numpy();groups=m.set_index('ligand_id').loc[df.ligand_id,'cluster_id'].to_numpy()
    changed=df.copy();ix=np.setdiff1d(np.arange(len(df)),train);changed.loc[ix,TARGET]=np.random.default_rng(99).normal(100,20,len(ix))
    temp=OUT/'verification';temp.mkdir(exist_ok=True);screen(changed,train,groups,temp)
    old=pd.read_csv(rd/'marker_screening.csv').sort_values('marker').reset_index(drop=True);new=pd.read_csv(temp/'marker_screening.csv').sort_values('marker').reset_index(drop=True);pd.testing.assert_frame_equal(old,new);(temp/'marker_screening.csv').unlink()
    global_count=0
    for f in (OUT/'global_models').glob('seed*/*/complete.json'):
        cfg=json.loads(f.read_text());p=pd.read_csv(f.parent/'test_predictions.csv');ca=pd.read_csv(f.parent/'calibration_predictions.csv');inner=pd.read_csv(f.parent/'inner_manifest.csv')
        assert len(inner)==cfg['n_train'] and len(ca)==cfg['n_calibration'] and len(p)==cfg['n_test']
        assert not set(inner.ligand_id)&set(p.ligand_id) and not set(inner.ligand_id)&set(ca.ligand_id) and not set(ca.ligand_id)&set(p.ligand_id)
        assert np.isclose(cfg['rmse'],np.sqrt(((p.y-p.prediction)**2).mean()))
        for level in [.8,.9,.95]:
            q=np.partition(ca.abs_residual.to_numpy(),int(np.ceil((len(ca)+1)*level))-1)[int(np.ceil((len(ca)+1)*level))-1]
            assert np.allclose(p[f'width_{level}'],2*q)
            assert np.array_equal(p[f'covered_{level}'],(p.abs_residual<=q).astype(int))
        manifest=pd.read_csv(OUT/f"runs/butina04_seed{cfg['seed']}/split_manifest.csv").set_index('ligand_id');innergroups={r:set(manifest.loc[g.ligand_id,'cluster_id']) for r,g in inner.groupby('inner_role')};assert not innergroups['fit']&innergroups['early_stopping']
        global_count+=1
    assert global_count==25,global_count
    dump(OUT/'verification.json',{'confidence_runs':checks,'global_model_runs_verified':global_count,'n_primary':len(df),'identity_manifest_checked':True,'holdout_label_perturbation_does_not_change_marker_selection':True,'quantile_oracles':True,'repaired_dataset_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'note':'These checks verify recorded computation and cohort accounting, not biological correctness or prospective validity.'})
    print('VERIFICATION PASS: 18 confidence runs + 25 global feature-set fits; identity accounting; cluster partitions; interval calculations; training-only marker selection.',flush=True)

if __name__=='__main__':
    summarize(OUT);summarize(OUT/'strict_identity');verify()
    import report;report.main()
