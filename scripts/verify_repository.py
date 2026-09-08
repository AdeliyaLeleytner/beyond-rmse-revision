"""Verify distributed checksums, notebook execution and key scientific accounting."""
from pathlib import Path
import hashlib,json
import nbformat
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
def main():
    hashes=ROOT/'SHA256SUMS.txt';checked=0
    if not hashes.exists():raise FileNotFoundError('Missing SHA256SUMS.txt')
    for line in hashes.read_text().splitlines():
        expected,relative=line.split('  ',1);path=ROOT/relative
        assert path.is_file(),relative
        assert hashlib.sha256(path.read_bytes()).hexdigest()==expected,relative
        checked+=1
    nbs=sorted((ROOT/'notebooks').glob('*.ipynb'));assert len(nbs)==13
    code_cells=figures=0
    for p in nbs:
        nb=nbformat.read(p,as_version=4);nbformat.validate(nb)
        for c in nb.cells:
            if c.cell_type=='code':
                assert c.execution_count is not None,p.name
                code_cells+=1
                for o in c.outputs:
                    assert o.output_type!='error',(p.name,o)
                    figures+=int('image/png' in o.get('data',{}))
    repaired=ROOT/'results/repaired';pair=pd.read_csv(repaired/'strict_global_comparison.csv')
    assert len(pair)==10 and (pair.rmse_ADME<pair.rmse_Plain).all()
    for cohort,n in [('',12584),('strict_identity',12193)]:
        for seed in range(42,47):
            m=pd.read_csv(repaired/cohort/f'runs/butina04_seed{seed}/split_manifest.csv')
            assert len(m)==n and m.ligand_id.is_unique
            groups={r:set(g.cluster_id) for r,g in m.groupby('role')}
            assert not groups['train']&groups['calibration']
            assert not groups['train']&groups['test']
            assert not groups['calibration']&groups['test']
    report={'files_verified':checked,'executed_notebooks':len(nbs),'executed_code_cells':code_cells,'embedded_figures':figures,'primary_and_strict_partitions_verified':10,'corrected_ADME_ranking_verified':True}
    (ROOT/'provenance/repository_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
