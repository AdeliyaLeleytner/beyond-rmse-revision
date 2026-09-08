from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
"""Recover ADMET identity before row filtering; retain audited chemical identities only."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from analysis import ROOT,OUT,REPO,dump

out=REPAIR_OUTPUT;out.mkdir(exist_ok=True)
basepath=SOURCE/'master_before_admet.csv.gz';srcpath=SOURCE/'prepared_admet.csv.gz';rawpath=SOURCE/'admet_raw.csv.gz'
base=pd.read_csv(basepath);src=pd.read_csv(srcpath);pub=pd.read_csv(REPO/'data/df_final.csv.gz')
raw=pd.read_csv(rawpath,usecols=['Cleaned SMILES','MW','hERG'],low_memory=False)
assert len(base)==len(raw)==12652
invalid=raw[['MW','hERG']].eq('Invalid Molecule').any(axis=1)
valid_idx=np.flatnonzero(~invalid.to_numpy());assert len(src)==len(valid_idx)
assert np.array_equal(raw.loc[~invalid,'Cleaned SMILES'].to_numpy(),src['Cleaned SMILES'].to_numpy())
assert set(pub.ligand_id)<=set(base.ligand_id)
basekeys=base.set_index('ligand_id')['Cleaned SMILES'];assert (pub['Cleaned SMILES'].to_numpy()==basekeys.loc[pub.ligand_id].to_numpy()).all()
uncharger=rdMolStandardize.Uncharger()
def key(smiles,neutral=False):
    m=Chem.MolFromSmiles(smiles)
    if m is None:return None
    if neutral:m=uncharger.uncharge(m)
    return Chem.MolToSmiles(m,isomericSmiles=True)
base_exact=base['Cleaned SMILES'].map(key);base_parent=base['Cleaned SMILES'].map(lambda s:key(s,True))
src_exact=src['Cleaned SMILES'].map(key);src_parent=src['Cleaned SMILES'].map(lambda s:key(s,True))
base_parent_counts=base_parent.value_counts();src_parent_counts=src_parent.value_counts()
records=[]
for j,i in enumerate(valid_idx):
    exact=base_exact.iloc[i]==src_exact.iloc[j]
    neutral=base_parent.iloc[i]==src_parent.iloc[j]
    unique_parent=(base_parent_counts[base_parent.iloc[i]]==1 and src_parent_counts[src_parent.iloc[j]]==1)
    status='exact_isomeric_identity' if exact else ('unique_charge_parent' if neutral and unique_parent else 'unresolved_structure_transform')
    records.append({'ligand_id':int(base.iloc[i].ligand_id),'master_row':int(i),'admet_valid_row':int(j),'raw_admet_row':int(i),'identity_status':status,'master_smiles':base.iloc[i]['Cleaned SMILES'],'admet_smiles':src.iloc[j]['Cleaned SMILES'],'exact_identity':exact,'same_uncharged_isomeric_parent':neutral,'unique_parent':unique_parent})
audit=pd.DataFrame(records);audit.to_csv(out/'identity_manifest.csv.gz',index=False)
# Verify that the published block is precisely the unkeyed sequence of source outputs.
shared=[c for c in src if c in pub and c!='Cleaned SMILES']
assert len(shared)==120
assert np.allclose(pub[shared].to_numpy(float),src[shared].to_numpy(float),equal_nan=True)
assert np.array_equal(pub['Cleaned SMILES.1'],src['Cleaned SMILES'])
attached_ids=base.iloc[valid_idx].ligand_id.to_numpy()
pub_audit=pd.DataFrame({'published_ligand_id':pub.ligand_id,'actual_source_ligand_id':attached_ids,'wrong_association':pub.ligand_id.to_numpy()!=attached_ids})
pub_audit.to_csv(out/'published_association_audit.csv.gz',index=False)
# Reconstruct by stable ID. No numerical imputation for missing identities.
source=src.copy();source.insert(0,'ligand_id',attached_ids);source=source.rename(columns={'Cleaned SMILES':'ADMET_source_SMILES'})
core=base.iloc[valid_idx].copy()
# Retain extra non-ADMET descriptors present only in the public matrix by ligand_id.
extra=[c for c in pub if c not in core and c not in shared and c!='Cleaned SMILES.1']
if extra:core=core.merge(pub[['ligand_id']+extra],on='ligand_id',how='left',validate='one_to_one')
core=core.drop(columns=[c for c in shared if c in core])
recovered=core.merge(source,on='ligand_id',validate='one_to_one').merge(audit[['ligand_id','identity_status']],on='ligand_id',validate='one_to_one')
for name,allowed in [('audited_parent_matched',['exact_isomeric_identity','unique_charge_parent']),('strict_exact_identity',['exact_isomeric_identity'])]:
    d=recovered[recovered.identity_status.isin(allowed)].sort_values('ligand_id').reset_index(drop=True)
    d.to_csv(out/(name+'.csv.gz'),index=False)
recovered.to_csv(out/'candidate_full_recovery_not_primary.csv.gz',index=False)
recovered[recovered.identity_status.eq('unresolved_structure_transform')][['ligand_id','Cleaned SMILES','ADMET_source_SMILES']].to_csv(out/'unresolved_identity.csv.gz',index=False)
meta={'premerge_master_n':len(base),'raw_admet_n':len(raw),'invalid_admet_ligand_ids':base.loc[invalid.to_numpy(),'ligand_id'].astype(int).tolist(),'published_master_missing_ids':sorted(set(base.ligand_id)-set(pub.ligand_id)),'published_wrong_ADMET_association_n':int(pub_audit.wrong_association.sum()),'published_wrong_association_first_last_id':[int(pub_audit.loc[pub_audit.wrong_association,'published_ligand_id'].min()),int(pub_audit.loc[pub_audit.wrong_association,'published_ligand_id'].max())],'identity_status_counts':audit.identity_status.value_counts().to_dict(),'primary_n':int(audit.identity_status.isin(['exact_isomeric_identity','unique_charge_parent']).sum()),'strict_n':int(audit.exact_identity.sum()),'extra_core_columns':extra,'raw_and_master_alignment_caveat':'Original API request manifest absent. Raw/master row pairing is supported by original cleaning code, matching lengths, exact structural anchors, and independently verified exact or unique uncharged isomeric identity for every retained row. Unresolved transformations are excluded, not presumed correct.','source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [basepath,srcpath,rawpath]}}
dump(out/'repair_summary.json',meta);print(json.dumps(meta,indent=2))
