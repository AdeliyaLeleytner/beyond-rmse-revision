"""Recompute the corrected study in a new directory; never overwrite committed results."""
from pathlib import Path
import argparse,datetime,hashlib,json,os,subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
def signature():
    paths=sorted([p for p in (ROOT/'analysis').glob('*.py') if p.name != 'notebook_support.py']+list((ROOT/'data').rglob('*')))
    return hashlib.sha256(''.join(str(p.relative_to(ROOT))+hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()).encode()).hexdigest()
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=ROOT/'reproduction'/datetime.datetime.now().strftime('%Y%m%d-%H%M%S'))
    ap.add_argument('--resume',action='store_true',help='Resume only if source and data hashes match.')
    args=ap.parse_args();out=args.output.resolve()
    if out==ROOT or out.is_relative_to(ROOT/'results') or out.is_relative_to(ROOT/'data'):ap.error('Use a separate reproduction/ directory.')
    sig=signature();state=out/'run_state.json'
    if out.exists():
        if not args.resume or not state.exists():ap.error('Output exists; select a new directory or explicitly --resume.')
        info=json.loads(state.read_text())
        if info['source_data_signature']!=sig:ap.error('Code/data changed: cached outputs cannot be reused.')
    else:
        out.mkdir(parents=True);info={'source_data_signature':sig,'stages_completed':[],'status':'running'}
    env=os.environ.copy();env.update(REANALYSIS_OUT=str(out/'repaired'),REVISION_REPAIR_OUT=str(out/'data_repair'),REVISION_ADDITIONAL_OUT=str(out/'additional'),REVISION_FIGURES_OUT=str(out/'figures'))
    def stage(name,script):
        if name in info['stages_completed']:return
        state.write_text(json.dumps(info,indent=2)+'\n')
        print('\nSTAGE:',name,flush=True)
        subprocess.run([sys.executable,str(ROOT/'analysis'/script)],cwd=ROOT,env=env,check=True)
        info['stages_completed'].append(name);state.write_text(json.dumps(info,indent=2)+'\n')
    stage('identity_repair','repair_data.py')
    env['REVISION_PROCESSED_DIR']=str(out/'data_repair')
    env['REANALYSIS_DATA']=str(out/'data_repair/audited_parent_matched.csv.gz')
    for name,script in [('structural_splits_and_confidence','run_repaired.py'),('primary_global_and_strict_confidence','extended_repaired.py'),('strict_global','strict_global.py'),('verification_and_maps','finalize.py'),('global_figures_and_SHAP','global_figures.py'),('strong_predictor_calibration','strong_predictor_calibration.py'),('strict_and_fingerprint_checks','extra_checks.py'),('incremental_docking','incremental_docking_check.py'),('reviewer_additional_checks','additional_checks.py'),('readable_figure4','figure4_readable.py')]:stage(name,script)
    info['status']='complete';state.write_text(json.dumps(info,indent=2)+'\n')
    print('Completed:',out)
if __name__=='__main__':main()
