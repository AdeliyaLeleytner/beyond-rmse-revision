from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
from pathlib import Path
import hashlib,json,platform,time
import numpy as np,pandas as pd
from rdkit import Chem,DataStructs,rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from sklearn.metrics import adjusted_rand_score
from analysis import OUT,REPO,DATA,load,dump

def main():
    out=OUT/'clusters';out.mkdir(exist_ok=True)
    if (out/'provenance.json').exists():print('CACHED clustering');return
    df=load().sort_values('ligand_id').reset_index(drop=True)
    cols=sorted([c for c in df if c.startswith('FP_')],key=lambda c:int(c.split('_')[1]))
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False,countSimulation=False)
    fps=[];mismatch=[]
    for i,row in enumerate(df[['ligand_id','Cleaned SMILES']].itertuples(index=False,name=None)):
        mol=Chem.MolFromSmiles(row[1])
        if mol is None:raise ValueError(f'Invalid cleaned SMILES ligand_id={row[0]}')
        fp=gen.GetFingerprint(mol);fps.append(fp)
        bits=np.zeros(2048,dtype=np.int8);DataStructs.ConvertToNumpyArray(fp,bits)
        if not np.array_equal(bits,df.iloc[i][cols].to_numpy(np.int8)):mismatch.append(int(row[0]))
    print('Parsed',len(fps),'mismatched stored fingerprints',len(mismatch),flush=True)
    n=len(fps);dist=np.empty(n*(n-1)//2,dtype=np.float64)
    start=0
    for i in range(1,n):
        dist[start:start+i]=1-np.asarray(DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i]));start+=i
    with (out/'fingerprints.bin').open('wb') as f:
        for fp in fps:f.write(DataStructs.BitVectToBinaryText(fp))
    meta={'rdkit_version':rdBase.rdkitVersion,'input_order':'ascending ligand_id','input_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'n':n,'fingerprint':{'type':'Morgan bit vector','radius':2,'fpSize':2048,'includeChirality':False,'countSimulation':False},'distance':'1 - Tanimoto similarity','reordering':False,'stored_fingerprint_mismatch_count':len(mismatch),'stored_fingerprint_mismatch_ids':mismatch,'cutoffs':{},'historical_threshold_status':'unknown; these are newly documented retrospective reanalysis partitions'}
    for cutoff in [.3,.4,.5]:
        raw=Butina.ClusterData(dist,n,cutoff,isDistData=True,reordering=False)
        labels=np.empty(n,dtype=int)
        for members in raw:
            # Stable cluster IDs: minimum ligand_id among members, no dependence on enumeration.
            labels[list(members)]=int(df.iloc[list(members)].ligand_id.min())
        counts=pd.Series(labels).value_counts()
        pd.DataFrame({'ligand_id':df.ligand_id,'cluster_id':labels}).to_csv(out/f'butina_distance_{cutoff:.1f}.csv',index=False)
        meta['cutoffs'][str(cutoff)]={'n_clusters':len(raw),'singleton_clusters':int((counts==1).sum()),'max_cluster_size':int(counts.max()),'ari_with_legacy':float(adjusted_rand_score(df.Butina_clusters,labels))}
        print('CLUSTERED',cutoff,meta['cutoffs'][str(cutoff)],flush=True)
    dump(out/'provenance.json',meta)
if __name__=='__main__':main()
