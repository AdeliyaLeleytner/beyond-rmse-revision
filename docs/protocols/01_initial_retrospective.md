# Retrospective reanalysis protocol — fixed before corrected-model results

Date: 2026-09-07. Scope: Reviewer 1, received part of comment 1. Preserve public snapshot and original outputs. No external submission. Full dataset: 12,651 molecules, including ligand_id 1298. Original replication alone reproduces its exclusion; a second evaluation restores it without retraining because it belongs only to the original test set.

Point model remains the published Plain CatBoost: Morgan bits + MW + logP; 800 iterations, depth 6, learning rate 0.05, seed 42, CPU. Median imputation is fit on proper training. No test-based hyperparameter optimization. Package versions locked separately.

## Evaluation sequence

1. Reproduce original molecule assignments, predictor and interval calculations. Compare individual predictions and widths to public CSVs. Evaluate the excluded molecule with that same predictor.
2. Refit confidence model on public split_index=1, with val used as calibration. All clusters are disjoint; its historical cluster-generation cutoff remains unknown.
3. Primary reproducible structural robustness experiment: regenerate Morgan fingerprints from Cleaned SMILES, radius 2, 2,048 bits, chirality off and count simulation off, ascending ligand_id. RDKit Butina, distance = 1 - Tanimoto, distance cutoff 0.4, reordering=False. The value is an explicitly NEW analysis parameter, not a claim about historical code. GroupShuffleSplit: 20% groups test, then 20% remaining groups calibration. Split seeds 42,43,44,45,46; fixed model seed42. No selection among seeds by outcomes.
4. Matched-size molecule-random controls for those five seeds, using exactly the structural split's train/calibration/test counts. These characterize split difficulty; different test membership prevents interpreting their metric difference as a pure causal leakage estimate.
5. Cutoff sensitivity: distance0.3 and0.5, seed42, unchanged model and interval settings. These outcomes will not select the primary cutoff.
6. Same-test calibration-allocation control: retain public split1 test and development pool; reassign the same number of development molecules to calibration individually with default_rng42. Compare with public split1 clusterwise allocation, acknowledging that proper-training membership changes too.

## Biological coordinates and uncertainty

For corrected runs, phenotype screening uses only proper-training labels. Candidate list is exactly the19 supplementaryS3 markers. Thirty internal group-holdout80/20 partitions (seed20260907) estimate train and internal-validation Spearman correlations. Select the three highest mean training correlations greater than0, retaining at most one of hERG/hERG-10um. Internal validation is descriptive and cannot select features. This is a documented training-only rule, not claimed to recover undocumented historical selection. If selected markers differ, record them; do not force the published trio post hoc. This remains retrospective, since the candidate list/hypothesis was informed by the previous full-data study.

Fit biological-coordinate means/SDs and quadrant thresholds on proper training for corrected runs. Calibration labels are used only for residual quantiles. Original replication retains original calibration-based coordinate scaling and original whole-table quadrant thresholds.

Global interval: corrected implementation uses exact order statistic rank ceil((n_cal+1)*0.9). Original replication retains NumPy's original quantile(method='higher') convention; report any difference. Local interval: unchanged Gaussian weighting, h0.9, weighted90th residual quantile, fallback to global if neff<80. Report it as empirical local residual calibration, not a proved local conditional-coverage guarantee. Seed42 runs additionally report h0.6/h1.2 sensitivity without selecting h from test outcomes.

Compare global/joint/pheno-only/MIE-only intervals on identical test molecules within every run: coverage, full width, interval score at alpha0.1, and paired differences. Report test-cluster bootstrap95% intervals (2,000 resamples), conditional on fitted model and calibration sample. Across five partition seeds report all outcomes and descriptive mean/SD; do not treat overlapping test sets as independent experiments for a significance test. These resampling summaries do not establish conformal exchangeability or universal coverage guarantees.

Save all assigned IDs and cluster IDs, model files, calibration/test predictions, selected features, screening tables, scaling parameters, metrics, checksums and exact package versions. Verify disjoint IDs and clusters and no unaccounted compound in each corrected cluster run. Assess cross-partition nearest-neighbour Tanimoto directly; different Butina labels do not bound every pair's similarity.

## Additional controls, specified before their outcomes

The independent consultation recommended two controls directly relevant to the user's training-effect question and the biological interpretation. On original split0's complete test set (2,531 molecules), compare: original train; original train with all test-cluster members removed; five random training subsamples of that same reduced size (seeds42–46). Keep model settings fixed, evaluate paired RMSE differences, and bootstrap whole test clusters. This targets the effect of training analogues beyond reducing training-set size; it does not turn a retrospective analysis into prospective validation.

For each of the five primary cluster splits, compare the joint biological localizer with prediction-only and MW/logP localizers, plus five permutations of biological coordinates independently within calibration and test (seeds1000–1004). Preserve point predictions and global quantile. Compare interval scores, coverage and width; do not select a model using these test results.

Screening implementation correction before interpretation: the four supplementary candidates encoded as ordinal values above1 must be screened in their original scale, not clipped to probabilities. All screenings are recomputed using the raw candidate values; the three selected phenotype probabilities remain clipped only when constructing biological coordinates. If selection changes, dependent results must be recomputed.
