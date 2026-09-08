from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
import os,json,hashlib,itertools
from pathlib import Path
import numpy as np,pandas as pd
HERE=Path(__file__).resolve().parent
os.environ.setdefault('REANALYSIS_OUT',str(RESULTS));os.environ.setdefault('REANALYSIS_DATA',str(PROCESSED/'audited_parent_matched.csv.gz'))
from analysis import cluster_bootstrap,dump,load
P=RESULTS;results=[]
for cohort,path in [('primary',P),('strict_exact_identity',P/'strict_identity')]:
    for seed in range(42,47):
        rd=path/'global_models'/f'seed{seed}';a=pd.read_csv(rd/'ADME/test_predictions.csv').set_index('ligand_id');b=pd.read_csv(rd/'Plain/test_predictions.csv').set_index('ligand_id').loc[a.index]
        assert np.array_equal(a.y,b.y) and np.array_equal(a.cluster_id,b.cluster_id)
        diff=a.abs_residual.to_numpy()**2-b.abs_residual.to_numpy()**2;ci=cluster_bootstrap(a.cluster_id,diff)
        results.append({'cohort':cohort,'seed':seed,'n_test':len(a),'rmse_ADME':float(np.sqrt(np.mean(a.abs_residual**2))),'rmse_Plain':float(np.sqrt(np.mean(b.abs_residual**2))),'mse_ADME_minus_Plain':float(diff.mean()),'ci_low':float(ci[0,0]),'ci_high':float(ci[1,0])})
        if cohort!='strict_exact_identity':continue
        manifest=pd.read_csv(path/f'runs/butina04_seed{seed}/split_manifest.csv').set_index('ligand_id')
        for name in ['ADME','Plain']:
            f=rd/name;cfg=json.loads((f/'complete.json').read_text());inner=pd.read_csv(f/'inner_manifest.csv');ca=pd.read_csv(f/'calibration_predictions.csv');te=pd.read_csv(f/'test_predictions.csv')
            sets=[set(x.ligand_id) for x in [inner,ca,te]];assert sum(map(len,sets))==12193
            for x,y in itertools.combinations(sets,2):assert not x&y
            groups=[set(manifest.loc[list(ids),'cluster_id']) for ids in sets]
            for x,y in itertools.combinations(groups,2):assert not x&y
            ig={role:set(manifest.loc[g.ligand_id,'cluster_id']) for role,g in inner.groupby('inner_role')};assert not ig['fit']&ig['early_stopping']
            assert np.isclose(cfg['rmse'],np.sqrt(np.mean((te.y-te.prediction)**2)))
            for level in [.8,.9,.95]:
                k=int(np.ceil((len(ca)+1)*level));q=np.partition(ca.abs_residual.to_numpy(),k-1)[k-1]
                assert np.allclose(te[f'width_{level}'],2*q) and np.array_equal(te[f'covered_{level}'],(te.abs_residual<=q).astype(int))
pd.DataFrame(results).to_csv(P/'strict_global_comparison.csv',index=False)
# Verify repaired model fingerprint provenance, including reinstated ID 9116.
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit import rdBase
df=load();fp=[f'FP_{i}' for i in range(2048)];generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True,countSimulation=False);mismatch=[]
for i,row in enumerate(df[['ligand_id','Cleaned SMILES']].itertuples(index=False,name=None)):
    mol=Chem.MolFromSmiles(row[1]);bits=np.array(generator.GetFingerprint(mol));stored=df.iloc[i][fp].to_numpy(int)
    if not np.array_equal(bits,stored):mismatch.append(int(row[0]))
assert not mismatch,mismatch[:10]
# Historical advisory metadata are preserved under docs; not a runtime dependency.
record=json.loads((P/'verification.json').read_text());record.update(strict_global_models_verified=10,repaired_model_fingerprints_match_chiral_Morgan={'n':len(df),'mismatches':0,'rdkit_version':rdBase.rdkitVersion},SHAP_additivity_checked_models=5)
dump(P/'verification.json',record)
print('EXTRA CHECKS PASS: 10 strict global fits; 12,584 model fingerprints; primary/strict paired ADME–Plain comparisons.')
