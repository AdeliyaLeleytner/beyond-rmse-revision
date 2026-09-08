"""Post-reanalysis claim check: joint vs phenotype-only localization, no model refitting."""
from pathlib import Path
import pandas as pd
from analysis import cluster_bootstrap
HERE=Path(__file__).resolve().parent
rows=[]
for seed in range(42,47):
    pred=pd.read_csv(HERE/f'repaired/runs/butina04_seed{seed}/test_predictions.csv')
    delta=(pred.interval_score_joint-pred.interval_score_pheno).to_numpy()
    ci=cluster_bootstrap(pred.cluster_id,delta)
    rows.append({'seed':seed,'joint_minus_pheno_score':float(delta.mean()),'ci_low':float(ci[0,0]),'ci_high':float(ci[1,0])})
pd.DataFrame(rows).to_csv(HERE/'repaired/docking_incremental_localization.csv',index=False)
