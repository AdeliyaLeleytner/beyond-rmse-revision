"""End-to-end repaired-cohort analysis, independent of published-matrix diagnostics."""
import json,os,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
os.environ['REANALYSIS_DATA']=str(HERE/'data_repair/audited_parent_matched.csv')
os.environ['REANALYSIS_OUT']=str(HERE/'repaired')
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
