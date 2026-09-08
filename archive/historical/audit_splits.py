"""Read-only audit of the public Beyond-RMSE snapshot; writes audit artifacts only."""
from pathlib import Path
import itertools, json, hashlib
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]
REPO = ROOT / 'tmp/reviewer1_audit/Beyond-RMSE'
OUT = Path(__file__).resolve().parent
cols = ['ligand_id','Pubchem CID','Butina_clusters','Cleaned SMILES','-lgLD50, mol/kg']
df = pd.read_csv(REPO/'data/df_final.csv', usecols=cols)
registry = pd.read_csv(REPO/'results/splits/split_indices.csv')
ex = pd.read_csv(REPO/'data/excluded_molecules.csv')
excluded = set(ex.loc[ex.scope.eq('all_model_splits'),'ligand_id'])
report = {'public_commit':'23b50718f8eedbfb2eb8dc33067dfb0887b159b2','dataset_rows':len(df),'unique_ligand_ids':int(df.ligand_id.nunique()),'missing_cluster_labels':int(df.Butina_clusters.isna().sum()),'exclusions':ex.to_dict('records'),'experiments':{}}
for si in [0,1]:
    sets = {k:g.ligand_id.to_numpy() for k,g in registry[registry.split_index.eq(si)].groupby('set')}
    raw_counts = {k:len(v) for k,v in sets.items()}
    # Notebook 05 applies exclusions; notebook 02 loads the unfiltered registry.
    if si == 0:
        ids = np.array([i for i in sets['train'] if i not in excluded])
        np.random.default_rng(42).shuffle(ids)
        n_cal = int(round(len(ids)*.2))
        sets = {'train':np.sort(ids[n_cal:]),'calibration':np.sort(ids[:n_cal]),'test':np.array([i for i in sets['test'] if i not in excluded])}
    joined = pd.concat([pd.DataFrame({'set':k,'ligand_id':ids}) for k,ids in sets.items()],ignore_index=True).merge(df,on='ligand_id',how='left',validate='one_to_one')
    assert joined.Butina_clusters.notna().all()
    pairs=[]
    for a,b in itertools.combinations(sets,2):
        da=joined[joined['set'].eq(a)]; db=joined[joined['set'].eq(b)]
        overlap=set(da.Butina_clusters)&set(db.Butina_clusters)
        pairs.append({'pair':f'{a}/{b}','shared_ligand_ids':len(set(sets[a])&set(sets[b])),'shared_clusters':len(overlap),'a_molecules_in_shared_clusters':int(da.Butina_clusters.isin(overlap).sum()),'b_molecules_in_shared_clusters':int(db.Butina_clusters.isin(overlap).sum())})
    report['experiments'][str(si)]={'raw_counts':raw_counts,'effective_counts':{k:len(v) for k,v in sets.items()},'unassigned_ids':sorted(set(df.ligand_id)-set(joined.ligand_id)),'pairs':pairs}
    joined[['ligand_id','set','Butina_clusters']].to_csv(OUT/f'audited_split_{si}.csv',index=False)
preds = pd.read_csv(REPO/'results/tables/conformal_test_predictions.csv')
report['prediction_table']={'rows':len(preds),'columns':list(preds.columns),'covered':int(preds.covered_90.sum()),'coverage':float(preds.covered_90.mean()),'mean_width':float(preds.local_interval_width_90.mean())}
local = pd.read_csv(ROOT/'Toxicodynamics-aware-data/processed/df_final.csv',usecols=['ligand_id','Butina_clusters'])
x = local.merge(df[['ligand_id','Butina_clusters']],on='ligand_id',suffixes=('_local','_public'))
report['local_public_comparison']={'local_rows':len(local),'same_id_set':set(local.ligand_id)==set(df.ligand_id),'same_cluster_labels':bool((x.Butina_clusters_local==x.Butina_clusters_public).all())}
for role in ['train','val','test']:
    old=pd.read_csv(ROOT/'Toxicodynamics-aware-data/splits'/f'{role}.csv')
    ids=old.ligand_id if 'ligand_id' in old else old.iloc[:,0]
    public_ids=registry.loc[registry.split_index.eq(1)&registry['set'].eq(role),'ligand_id']
    report['local_public_comparison'][role+'_matches_public_split1']=set(ids)==set(public_ids)
report['file_sha256']={str(f.relative_to(REPO)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [REPO/'data/df_final.csv',REPO/'data/excluded_molecules.csv',REPO/'results/splits/split_indices.csv',REPO/'notebooks/05_fig5_conformal_landscape.ipynb']}
(OUT/'audit_results.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
print(json.dumps(report,indent=2,ensure_ascii=False))
