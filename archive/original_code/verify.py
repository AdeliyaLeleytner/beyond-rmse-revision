import hashlib,itertools,json
import numpy as np,pandas as pd
from analysis import OUT,REPO,load,screen,global_quantile,local_intervals

def main():
    assert global_quantile(np.arange(100))==90
    w,ne=local_intervals(np.zeros((1,1)),np.zeros((100,1)),np.arange(100),100)
    assert w[0]==178 and ne[0]==100
    w,ne=local_intervals(np.zeros((1,1)),np.zeros((2,1)),np.array([.2,.8]),1)
    assert w[0]==2 and ne[0]==2
    df=load();records=[]
    for rd in sorted((OUT/'runs').iterdir()):
        if not (rd/'complete.json').exists():continue
        cfg=json.loads((rd/'complete.json').read_text());m=pd.read_csv(rd/'split_manifest.csv');p=pd.read_csv(rd/'test_predictions.csv');c=pd.read_csv(rd/'calibration_predictions.csv')
        assert m.ligand_id.is_unique
        roles={r:set(g.ligand_id) for r,g in m.groupby('role')}
        for a,b in itertools.combinations(roles,2):assert not roles[a]&roles[b]
        assert set(p.ligand_id)==roles['test'] and set(c.ligand_id)==roles['calibration']
        assert len(m)==(12650 if rd.name=='original_reproduction' else 12651)
        if rd.name.startswith('butina') or rd.name=='legacy_groups_registry1':
            clusters={r:set(g.cluster_id) for r,g in m.groupby('role')}
            for a,b in itertools.combinations(clusters,2):assert not clusters[a]&clusters[b]
        for row in pd.read_csv(rd/'interval_metrics.csv').itertuples():
            assert np.isclose(row.coverage,p['covered_'+row.method].mean())
            assert np.isclose(row.mean_width,p['width_'+row.method].mean())
            assert np.allclose(p['covered_'+row.method],p.abs_residual<=p['width_'+row.method]/2)
        if rd.name=='original_reproduction':
            assert cfg['reproduction']['max_prediction_difference']<1e-10
            assert cfg['reproduction']['max_joint_width_difference']<1e-10
            old=pd.read_csv(REPO/'results/tables/conformal_quadrants.csv').sort_values('Quadrant')
            new=pd.read_csv(rd/'quadrants.csv').sort_values('quadrant')
            assert old.n.tolist()==new.n.tolist()
            assert np.allclose(old['Joint width'],new.joint_width)
        records.append({'run':rd.name,'verified':True,'n_assigned':len(m)})
    assert len(records)==16,len(records)
    # Perturb all nontraining target labels; training-only marker selection must not change.
    rd=OUT/'runs/butina04_seed42';m=pd.read_csv(rd/'split_manifest.csv');train=m.loc[m.role.eq('train'),'row_index'].to_numpy()
    groups=np.empty(len(df),dtype=int)
    for r in m.itertuples():groups[r.row_index]=r.cluster_id
    perturbed=df.copy();ix=np.setdiff1d(np.arange(len(df)),train)
    perturbed.loc[ix,'-lgLD50, mol/kg']=np.random.default_rng(99).normal(100,20,len(ix))
    temp=OUT/'verification';temp.mkdir(exist_ok=True)
    selected=screen(perturbed,train,groups,temp)
    original=pd.read_csv(rd/'marker_screening.csv').sort_values('marker').reset_index(drop=True)
    altered=pd.read_csv(temp/'marker_screening.csv').sort_values('marker').reset_index(drop=True)
    pd.testing.assert_frame_equal(original,altered)
    (temp/'marker_screening.csv').unlink()
    payload={'runs':records,'exact_baseline_reproduced':True,'local_quantile_oracle_tests':True,'selection_invariant_to_calibration_and_test_targets':True,'public_dataset_sha256':hashlib.sha256((REPO/'data/df_final.csv').read_bytes()).hexdigest()}
    (OUT/'verification.json').write_text(json.dumps(payload,indent=2)+'\n')
    print('PASS:16runs, disjointness, accounting, numerical outputs, baseline replication, quantile oracles, training-only selection.')
if __name__=='__main__':main()
