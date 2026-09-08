from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
from paths import RESULTS, PROCESSED
import extended_repaired as e
import pandas as pd
e.a.OUT=RESULTS/'strict_identity'
e.global_models(pd.read_csv(PROCESSED/'strict_exact_identity.csv.gz'),['Plain','ADME'])
