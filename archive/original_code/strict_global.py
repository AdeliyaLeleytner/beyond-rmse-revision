"""Check the new ADME/Plain result without any charge-parent matching."""
import pandas as pd
import extended_repaired as e
e.a.OUT=e.HERE/'repaired/strict_identity'
e.global_models(pd.read_csv(e.HERE/'data_repair/strict_exact_identity.csv'),['Plain','ADME'])
