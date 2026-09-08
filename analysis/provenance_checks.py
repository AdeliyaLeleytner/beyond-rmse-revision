from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
import sys,json,hashlib
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.model_selection import train_test_split,GroupShuffleSplit
from analysis import OUT,REPO,load,dump

df=load();fps=df[[f'FP_{j}' for j in range(2048)]].to_numpy(np.int8)
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
mismatches=[]
for i,smiles in enumerate(df['Cleaned SMILES']):
    fp=gen.GetFingerprint(Chem.MolFromSmiles(smiles));a=np.zeros(2048,dtype=np.int8);DataStructs.ConvertToNumpyArray(fp,a)
    if not np.array_equal(a,fps[i]):mismatches.append(int(df.iloc[i].ligand_id))
r=pd.read_csv(REPO/'results/splits/split_indices.csv');target=set(r.loc[r.split_index.eq(0)&r['set'].eq('test'),'ligand_id'])
matches=[]
for order,ids in [('csv_order',df.ligand_id.to_numpy()),('ascending_ligand_id',np.sort(df.ligand_id.to_numpy()))]:
    for seed in range(100):
        tr,te=train_test_split(ids,test_size=.2,random_state=seed)
        if set(te)==target:matches.append({'ordering':order,'generator':'sklearn.train_test_split(test_size=0.2)','seed':seed})
        perm=np.random.default_rng(seed).permutation(ids)
        if set(perm[:len(target)])==target:matches.append({'ordering':order,'generator':'numpy.default_rng(seed).permutation; first n_test','seed':seed})
groups=df.Butina_clusters.to_numpy();development,test=next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=42).split(df,groups=groups))
train,cal=next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=42).split(development,groups=groups[development]))
checks={}
for role,ix in [('train',development[train]),('val',development[cal]),('test',test)]:
    saved=set(r.loc[r.split_index.eq(1)&r['set'].eq(role),'ligand_id']);reconstructed=set(df.iloc[ix].ligand_id)
    checks[role]={'n_reconstructed':len(ix),'exact_id_set_match':saved==reconstructed,'only_saved':sorted(saved-reconstructed),'only_reconstructed':sorted(reconstructed-saved)}
info={'stored_Morgan_includeChirality_true_mismatch_count':len(mismatches),'mismatch_ids':mismatches,'random_split_reconstruction_matches':matches,'legacy_group_split_seed42_reconstruction':checks,'historical_provenance_caveat':'Matching assignments is a computational reconstruction, not proof of the original execution history.'}
if 'MW, g/mol' in df:
    diff=df['MW']-df['MW, g/mol'];info['MW_vs_MW_g_per_mol']={'max_abs_difference':float(np.abs(diff).max()),'mean_abs_difference':float(np.abs(diff).mean()),'exact_equal_fraction':float((diff==0).mean())}
dump(OUT/'provenance_reconstruction.json',info);print(json.dumps(info,indent=2))
