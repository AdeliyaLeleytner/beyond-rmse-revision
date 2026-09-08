# Changes for the reviewer revision

## Data integrity

- Replace positional ADMET association with recovered IDs and independent structural identity checks.
- Restore valid ligand9116; exclude invalid ADMET ligand11434 and67 unresolved transformations.
- Retain ligand1298 in corrected experiments; preserve the original exclusion only in historical reproduction.
- Separate source master MW from ADMET MW, preserve their actual model definitions and record source hashes.

## Evaluation and preprocessing

- Document newly generated Butina clusters and whole-cluster proper-training/calibration/test assignment across five seeds.
- Fit corrected coordinate scaling, marker selection and quadrant boundaries on proper training.
- Use the corrected global conformal order statistic; retain the historical quantile convention only in original reproduction.
- Fit imputation/PCA in the inner fitting partition; select tree count with inner grouped early stopping and refit on proper training.
- Report individual predictions, observed coverage and interval score alongside widths; cluster-bootstrap uncertainty is conditional on fitted artifacts.

## Additional experiments

- Exact-identity sensitivity, cutoff/bandwidth/fallback sensitivity, size-matched molecular-random partitions.
- Fixed-test training-analogue removal with five equally sized random reductions.
- Prediction-only, MW/logP and shuffled biological-coordinate localizers.
- ADME versus Plain on corrected primary and strict datasets; unchanged localizer applied to both stronger global predictors.
- Matched BBB subgroup comparison, grouped feature permutations and SHAP additivity.
- Target-wise docking normalization, incremental mean44 beyond phenotype-only localization, VIF/correlation/residual-decile checks.
- Whole-cluster bootstrap of quadrant thresholds and narrower-than-global interval proportions; display support/masking inventory.

## Delivery

- Portable repository-relative code and compressed source/processed tables.
- Thirteen executed notebooks with English narrative comments, figures and explicit limits.
- Full fresh-run reproduction command, source/data cache signature, notebook executor, checksum verification and independent artifact comparison.
- Historical public-repository disagreements documented without treating them as corrected results.

No new docking, external ADMET API calls, Optuna search, weighted biological target model or independent external validation was performed.
