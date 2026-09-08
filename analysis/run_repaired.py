from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
"""End-to-end repaired-cohort analysis, independent of published-matrix diagnostics."""
import json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];HERE=Path(__file__).resolve().parent
os.environ.setdefault('REANALYSIS_DATA',str(PROCESSED/'audited_parent_matched.csv.gz'))
os.environ.setdefault('REANALYSIS_OUT',str(RESULTS))
from analysis import OUT,load,run,original_split
OUT.mkdir(exist_ok=True)
from prepare_clusters import main as prepare
prepare()
df=load()
# Same historical assignments, restricted to auditable identities, before re-splitting.
run('repaired_reference',df,original_split(df,True),df.Butina_clusters.to_numpy(),original=True,strict=False,screen_markers=False)
for family in ['new','random','sensitivity']:
 subprocess.run([sys.executable,str(HERE/'run_suite.py'),'--family',family],check=True,env=os.environ.copy())
from controls import ablation,null_controls,similarity
ablation(df);null_controls(df);similarity(df)
print('REPAIRED SUITE COMPLETE',flush=True)
