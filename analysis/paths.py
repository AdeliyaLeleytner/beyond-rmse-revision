"""Portable paths. Recomputed outputs must be separate from the committed evidence."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/source'
PROCESSED=Path(os.environ.get('REVISION_PROCESSED_DIR',ROOT/'data/processed')).resolve()
RESULTS=Path(os.environ.get('REANALYSIS_OUT',ROOT/'results/repaired')).resolve()
ADDITIONAL=Path(os.environ.get('REVISION_ADDITIONAL_OUT',ROOT/'results/additional')).resolve()
REPAIR_OUTPUT=Path(os.environ.get('REVISION_REPAIR_OUT',ROOT/'.runtime/data_repair')).resolve()
for directory in (RESULTS,ADDITIONAL,REPAIR_OUTPUT): directory.mkdir(parents=True,exist_ok=True)
