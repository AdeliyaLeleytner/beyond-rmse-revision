"""Compare recomputed numeric artifacts with the committed corrected evidence."""
from pathlib import Path
import argparse,json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('reproduction',type=Path);args=ap.parse_args()
    fresh=args.reproduction.resolve();checks=[];missing=[]
    for family in ['repaired','additional']:
        old=ROOT/'results'/family;new=fresh/family
        for source in sorted(old.rglob('*.csv')):
            target=new/source.relative_to(old)
            if not target.exists():missing.append(str(source.relative_to(ROOT)));continue
            a=pd.read_csv(source);b=pd.read_csv(target)
            pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=False,rtol=1e-8,atol=1e-10)
            checks.append({'path':str(source.relative_to(ROOT)),'rows':len(a),'numeric_and_ID_values_match':True})
    for name in ['audited_parent_matched','strict_exact_identity','identity_manifest','published_association_audit']:
        a=pd.read_csv(ROOT/f'data/processed/{name}.csv.gz');b=pd.read_csv(fresh/f'data_repair/{name}.csv.gz')
        pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=False,rtol=1e-9,atol=1e-11)
        checks.append({'path':f'data/processed/{name}.csv.gz','rows':len(a),'numeric_and_ID_values_match':True})
    report={'source':'Fresh independent output directory, complete pipeline execution','status':json.loads((fresh/'run_state.json').read_text())['status'],'compared_artifacts':len(checks),'checks':checks,'not_regenerated_by_pipeline':missing,'comparison_tolerance':{'rtol':1e-8,'atol':1e-10}}
    (ROOT/'provenance/portability_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Matched',len(checks),'artifacts; not regenerated:',missing)
if __name__=='__main__':main()
