# Beyond RMSE — corrected analysis and reviewer experiments

Executed Jupyter notebooks for the Scientific Reports revision of **Beyond RMSE: biology-informed uncertainty estimation for toxicological QSAR**.

This repository contains the corrected compound–ADMET associations, fitted models, train/calibration/test manifests, individual predictions, all additional reviewer experiments, and figures with editable SVG versions. Notebook comments and figure interpretations are in English; Russian author notes are in `docs/`.

## Start here

Open [00_start_here.ipynb](notebooks/00_start_here.ipynb). All notebooks include saved outputs; GitHub can display their tables and figures without running Python.

| Notebook | What it demonstrates |
|---|---|
| [00 Overview](notebooks/00_start_here.ipynb) | Updated results, scope and reading order |
| [01 Identity repair](notebooks/01_compound_identity_and_data_repair.ipynb) | The positional-join error, keyed reconstruction and chemical checks |
| [02 Structural partitions](notebooks/02_structural_splits_and_analogue_controls.ipynb) | No shared clusters; fixed-test analogue-removal control; nearest-neighbour similarity |
| [03 Global models](notebooks/03_global_models_and_revised_ranking.ipynb) | Direct prediction checks for 35 fitted models and revised ADME–Plain ranking |
| [04 Prediction intervals](notebooks/04_conformal_prediction_and_stronger_predictor.ipynb) | Global/local calibration, coverage, width and interval score, including ADME |
| [05 SHAP and permutations](notebooks/05_SHAP_and_group_permutations.ipynb) | Model reliance, grouped feature permutation and readable SHAP plots |
| [06 Matched BBB comparison](notebooks/06_BBB_matched_comparison.ipynb) | Specialist versus Plain/Baseline on identical test compounds |
| [07 Biological coordinates](notebooks/07_biological_coordinates_and_docking_normalization.ipynb) | BTox, VIF, target-wise normalization and incremental mean44 contribution |
| [08 Robustness](notebooks/08_permuted_coordinates_and_sensitivity.ipynb) | Shuffled coordinates, alternative localizers, cutoffs, bandwidth and fallback |
| [09 Applicability maps](notebooks/09_applicability_maps_and_threshold_stability.ipynb) | Display masking, median bootstrap, quadrant counts and uneven coverage |
| [10 Historical audit](notebooks/10_historical_and_related_repository_audits.ipynb) | Original reproduction and differences from the related public repository |
| [11 Reproduction](notebooks/11_reproduce_and_reviewer_map.ipynb) | Full rerun, verification and the 25-point reviewer map |
| [12 Figure gallery](notebooks/12_complete_figure_gallery.ipynb) | All final manuscript figures |

## Revised results

- The primary cohort contains **12,584** structurally checked associations; the exact-identity sensitivity contains **12,193**. The original public matrix had **1,678** shifted ADMET associations. No new ADMET API predictions were requested.
- Across five whole-cluster partitions, mean RMSE is **0.4930 for ADME** and **0.5122 for Plain**. ADME wins in all five primary and five strict partitions. These configurations retain different published hyperparameters; this is not an isolated causal feature-addition estimate.
- Biological localization also helps the stronger ADME predictor: mean interval score **2.3345 → 2.2571**, observed coverage **90.63% → 89.62%**. The stronger-predictor follow-up is explicitly retrospective/outcome-triggered.
- On the **same BBB-test compounds**, Plain has mean RMSE **0.4628**, versus **0.4759** for the BBB specialist. Specialist training is not supported as an improvement.
- Local intervals do **not** have a demonstrated uniform conditional-coverage guarantee. Docking and ADMET coordinates are computational proxies, not measured mechanisms.

The global-comparison Plain model uses master-table MW. The original fixed confidence predictor uses ADMET MW and different fixed hyperparameters; their results are labeled separately. This revision retains the article reanalysis docking values. The separate Toxicodynamics-aware repository clips positive docking scores to zero and is not numerically interchangeable with this pipeline.

## Install and replay

Python **3.12** was used. CPU execution is sufficient; no API credentials or GPU service is needed. Install the recorded environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/verify_repository.py
python scripts/execute_notebooks.py
```

Alternatively, open the notebooks in Jupyter, VS Code, or another notebook editor and choose the installed environment. Each notebook finds the repository root automatically. No file paths outside the repository are required.

`execute_notebooks.py` starts a fresh kernel for each notebook and writes standalone HTML previews to `.runtime/html/`. It verifies saved predictions, recomputes selected interval/permutation calculations and draws figures; it does **not** claim to retrain every model. Notebook-generated PNG/SVG graphs are also distributed in `results/notebook_figures/`.

## Full recomputation

```bash
python scripts/reproduce.py --output reproduction/my-run
python scripts/compare_reproduction.py reproduction/my-run
```

This reconstructs the audited data from the included source tables; computes Butina clusters; fits confidence and global models; runs exact-identity, random-size-matched and analogue controls; recalculates intervals, SHAP and grouped permutations; and renders figures. Pairwise clustering materializes approximately 633 MB of distances for the primary cohort; several GB of available RAM and additional disk space are useful. Runtime depends on the machine.

A new output directory is required. `--resume` is allowed only when code and source-data signatures match. Never reuse cached runs after altering the protocol, code or data. Committed `results/` are the evidence snapshot and are not overwritten by this runner. Historical flawed-matrix experiments are preserved in `archive/historical/`, clearly separated from revised results; Notebook 10 inspects their saved reproduction checks.

## Repository layout

- `notebooks/`: executed narrative analyses with tables and embedded figures.
- `analysis/`: portable executable research code used by the notebooks and full runner.
- `data/source/`: losslessly compressed source snapshots and original model settings.
- `data/processed/`: corrected/strict datasets, identity manifest, exclusion records and clearly named diagnostic-only candidate recovery.
- `results/repaired/`: primary, strict and sensitivity models, transformations, individual predictions, split manifests, tables and figures.
- `results/additional/`: normalization, VIF, fallback, threshold, BBB and permutation experiments.
- `results/figures/`: final readable Figure 4 panels; `results/notebook_figures/`: additional diagnostic figures.
- `docs/`: reviewer-to-experiment map, protocols, draft responses and manuscript replacements.
- `archive/`: historical results and untouched original research scripts for provenance.
- `provenance/`: source import hashes, execution records and fresh-run comparison.

## Scope and status

The study concerns the supplied experimental mouse intravenous LD50 collection. No new docking, external ADMET prediction, Optuna optimization or independent external validation is claimed. Historical clustering cutoff, complete study-source metadata and the exact upstream ADMET model/training provenance remain unresolved. Sixty-seven ambiguous structural transformations remain excluded.

The primary and additional protocols distinguish decisions made before their corrected results from subsequent diagnostics. Five overlapping partitions are summarized descriptively, not treated as independent studies. Cluster bootstrap intervals condition on the fitted model and calibration set.

The English/Russian author responses and manuscript replacement files are **drafts for integration**, not evidence of journal submission. The historical reports/protocols retain their original dates and chronology. Use this README and the executed notebooks for the current repository layout and status. Full third-party reviewer reports and Google Doc account metadata are not included.

## Provenance and integrity

Source snapshots: [Beyond-RMSE, 23b5071](https://github.com/chemagents/Beyond-RMSE/tree/23b50718f8eedbfb2eb8dc33067dfb0887b159b2) and [Toxicodynamics-aware-computational-toxicology, 5a3b2a2](https://github.com/chemagents/Toxicodynamics-aware-computational-toxicology/tree/5a3b2a257782e40538eea0d2d6ea940fa87ae822).

`provenance/import_manifest.json` traces copied artifacts to their source paths and original hashes. `SHA256SUMS.txt` covers the distributed repository files. `provenance/portability_validation.json` records comparisons against a fresh full execution. Preserve the Git commit when citing a version. No Zenodo DOI is claimed.

See [LICENSE](LICENSE) for code and [DATA_LICENSE.md](DATA_LICENSE.md) for source-data terms.
