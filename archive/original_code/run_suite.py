import argparse,json,time
import numpy as np,pandas as pd
from analysis import OUT,load,run,group_split,registry_split,dump

def cluster_labels(df,cutoff):
    r=pd.read_csv(OUT/'clusters'/f'butina_distance_{cutoff:.1f}.csv').set_index('ligand_id')
    return r.loc[df.ligand_id,'cluster_id'].to_numpy()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--family',choices=['legacy','new','random','sensitivity','all'],default='all');args=ap.parse_args()
    df=load();jobs=[]
    legacy=df.Butina_clusters.to_numpy()
    if args.family in ['legacy','all']:
        sp=registry_split(df)
        jobs.append(('legacy_groups_registry1',sp,legacy,True))
        pool=np.concatenate([sp['train'],sp['calibration']]);pool=np.sort(pool);np.random.default_rng(42).shuffle(pool)
        n=len(sp['calibration']);control={'train':np.sort(pool[n:]),'calibration':np.sort(pool[:n]),'test':sp['test']}
        jobs.append(('legacy_same_test_random_calibration',control,legacy,False))
    if args.family in ['new','random','all']:
        groups=cluster_labels(df,.4)
        for seed in range(42,47):
            split=group_split(df,groups,seed)
            if args.family in ['new','all']:jobs.append((f'butina04_seed{seed}',split,groups,True))
            if args.family in ['random','all']:
                ids=np.random.default_rng(seed).permutation(len(df));nt=len(split['test']);nc=len(split['calibration'])
                random={'train':np.sort(ids[nt+nc:]),'calibration':np.sort(ids[nt:nt+nc]),'test':np.sort(ids[:nt])}
                jobs.append((f'random_size_matched_seed{seed}',random,groups,False))
    if args.family in ['sensitivity','all']:
        for cutoff in [.3,.5]:
            groups=cluster_labels(df,cutoff)
            jobs.append((f'butina0{int(cutoff*10)}_seed42',group_split(df,groups,42),groups,True))
    for name,split,groups,strict in jobs:
        print('START',name,{k:len(v) for k,v in split.items()},flush=True)
        run(name,df,split,groups,strict=strict)
if __name__=='__main__':main()
