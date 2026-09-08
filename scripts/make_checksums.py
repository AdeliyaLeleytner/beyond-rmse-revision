from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[1]
IGNORED={'.git','.venv','.runtime','reproduction','__pycache__','.ipynb_checkpoints'}
EXCLUDED={'SHA256SUMS.txt','provenance/repository_validation.json'}
paths=sorted(p for p in ROOT.rglob('*') if p.is_file() and not IGNORED.intersection(p.relative_to(ROOT).parts) and str(p.relative_to(ROOT)) not in EXCLUDED and p.name!='.DS_Store')
(ROOT/'SHA256SUMS.txt').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(ROOT))+'\n' for p in paths))
print('Checksummed',len(paths),'files.')
