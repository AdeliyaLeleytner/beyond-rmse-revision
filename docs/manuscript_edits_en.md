# Text for the revised manuscript

These replacements are prepared for author integration. They have not been inserted into the manuscript source. Use them together with `reanalysis_methods_results_en.md`, which contains the complete data-accounting and reanalysis Methods. The new data and conclusions must replace the old results consistently throughout the Abstract, main text, supplement and captions.

## Proposed title

Beyond RMSE: biology-informed uncertainty estimation for toxicological QSAR

## Introduction opening and motivation

Acute toxicity is a whole-organism outcome that can arise through several molecular and physiological processes. The median lethal dose, LD₅₀, summarizes the dose associated with mortality in half of an experimental population under specified conditions. Acute-toxicity measurements also inform hazard classification; for example, OECD Test Guideline 423 describes classification using acute oral toxicity outcomes. That context motivates the study of toxicity prediction, but does not make different species or exposure routes interchangeable. Here we study experimental mouse intravenous LD₅₀. Our results are not a validation for human toxicity, chronic exposure or oral hazard classification. [OECD Test Guideline 423](https://www.oecd.org/en/publications/test-no-423-acute-oral-toxicity-acute-toxic-class-method_9789264071001-en.html).

Quantitative structure–activity relationship (QSAR) models estimate toxicological outcomes from molecular descriptors. Their average error is useful for comparing predictors, but can conceal substantial differences between compounds. A model may be relatively accurate for one part of a collection and much less reliable elsewhere, even when a single overall metric appears satisfactory. We use local heterogeneity to describe these differences in the distribution of prediction errors across descriptor space.

Biological context may help with both prediction and uncertainty estimation. In this study, that context is computational: predicted absorption, distribution, metabolism, excretion and toxicity (ADMET) descriptors and docking scores on a defined protein panel. These variables are not measurements of the corresponding phenotypes or proof of a molecular initiating event (MIE). We also examine a subgroup defined by physicochemical rules associated with blood–brain barrier (BBB) passage, motivated by possible differences in central nervous system (CNS) exposure. Passing those rules is not an experimental permeability measurement.

## Introduction final paragraph

We ask two related questions. First, do computational biological descriptors improve LD₅₀ point prediction relative to structural descriptors in the evaluated model configurations? Second, do they help characterize heterogeneous prediction errors and construct more useful prediction intervals? We compare structural, docking and ADME feature configurations, examine a BBB-rule-pass subgroup, and assess global and local residual calibration across repeated whole-cluster partitions. The Background introduces the modelling and uncertainty concepts, the Methods describe data assembly and evaluation, and the Results examine point prediction, feature reliance and interval behaviour before discussing their toxicological scope.

## Background transition and uncertainty

QSAR provides the link between molecular representation and the toxicological endpoint. Machine-learning methods can model complex relationships within that representation, but their predictive accuracy does not by itself describe the uncertainty of a new estimate. For a compound being considered for further testing, both the estimated toxicity and the plausible magnitude of error matter.

Split conformal prediction uses errors on a separate calibration sample to construct prediction intervals. Under the relevant exchangeability conditions, its standard coverage statement is marginal over calibration and new observations; it does not state that every individual compound or subgroup has that coverage. Heterogeneous sources do not automatically violate exchangeability, but changes in source composition, experimental conditions or dependence can invalidate an assumed transfer of coverage. [Angelopoulos and Bates](https://arxiv.org/html/2107.07511v6).

We compare a global residual interval with a local, kernel-weighted residual procedure in computational biological coordinates. The latter is evaluated empirically; no distribution-free conditional-coverage theorem is claimed for this implementation. Local residual patterns and the availability of relevant calibration observations offer a way to examine applicability-domain behaviour. They provide evidence about where additional validation may be needed, rather than a stand-alone assurance of safety.

## Evaluation metrics

We report several metrics because they answer different questions. Mean absolute error (MAE) measures the average absolute deviation. Mean squared error (MSE) averages squared deviations, and root mean squared error (RMSE) expresses its square root in the target's units, with greater sensitivity to large errors than MAE. R² compares squared prediction error with the variability of the observed endpoint. Spearman's ρ describes rank ordering, while Lin's concordance correlation coefficient (CCC) describes numerical agreement, including differences in location and scale. RMSE and MSE contain the same ranking information; the additional metrics help distinguish error magnitude, ordering and agreement.

The endpoint is y = −log₁₀(LD₅₀ [mol/kg]). A difference of 0.5 in y corresponds to a dose ratio of approximately 3.16. This is an illustration of scale; RMSE of 0.5 does not mean that every prediction has a 3.16-fold error. Comparisons across different test populations also require attention to endpoint range and variance.

Intervals are evaluated using observed coverage, full width and the interval score at nominal 90%: full width plus twenty times the amount by which the absolute error exceeds the half-width, when positive. Narrowing an interval is useful only when the resulting misses are also considered. Reported cluster-bootstrap intervals condition on the fitted model and calibration sample. The five structural partitions overlap and are summarized descriptively, not treated as independent experiments.

## Descriptor and BBB Methods clarification

Descriptor provenance is separated into source-table values, local calculations and external ADMETlab outputs. In the BBB procedure, existing master-table molecular weight (MW, g/mol) and logP are retained where present. Missing values are filled using RDKit MolWt and Crippen MolLogP. TPSA, hydrogen-bond donor and acceptor counts, and rotatable-bond counts use CalcTPSA, CalcNumHBD, CalcNumHBA and CalcNumRotatableBonds, respectively. For these local calculations, multi-fragment SMILES are represented by the valid fragment with the largest heavy-atom count. The recorded source environment specifies RDKit 2025.3.5; the repeat-analysis environment uses 2025.03.6. Upstream source-table generation settings and the precise external ADMETlab model version could not be independently recovered.

The BBB-rule-pass subgroup satisfies all of MW≤500 g/mol, 0≤logP≤5, TPSA≤90 Å², hydrogen-bond donors≤2 and rotatable bonds≤8. Acceptor count is computed but is not a classification condition. Nonfinite MW, logP or TPSA fails the rule. These rules reproduce all 12,584 saved classifications in the audited cohort. They define a computational subgroup, not measured BBB permeability or a separately validated BBB classifier.

Global-comparison models use the saved master-table MW column; the original confidence-predictor configuration uses the saved ADMET MW column. Both definitions are retained and documented. Existing model fingerprints were verified as chirality-aware Morgan radius-2, 2,048-bit fingerprints for all audited compounds; the separately specified clustering fingerprints are achiral.

## Normalization and display-coordinate sensitivity

The primary four-dimensional localization space contains raw panel-mean docking score and three separately standardized predicted toxicity descriptors. A target-wise-normalized sensitivity variant first centers and scales each of the 44 docking-score columns using proper training, averages those standardized scores, and standardizes the resulting aggregate on proper training before localization. All fitted predictors, partitions, residuals and other settings are held fixed. Across five partitions, mean interval score is 2.397 for raw mean44 and 2.409 for the normalized variant; coverage is 89.37% and 89.36%. This comparison does not establish a universally preferable aggregation, and it does not turn either coordinate into an affinity or mechanistic scale.

BTox is the sum of hERG, Respiratory and Neurotoxicity-DI probabilities for display. Calibration uses their three separate values. Equal display weights do not imply equal biological effects; different combinations can produce the same BTox. VIFs of the actual three coordinates range from 1.15 to 1.37 across proper-training partitions. The correlated hERG-10um endpoint is excluded from this input set. Held-out residual variance is not strictly monotonic across training-defined BTox deciles, so BTox is not described as a validated ordinal toxicity scale.

## Results on matched BBB compounds

To separate subgroup composition from specialist-model performance, the full-training-set Plain and Baseline models were evaluated on exactly the same held-out compounds as the BBB-rule-pass model. Mean RMSE across partitions is 0.4628 for Plain, 0.4767 for Baseline and 0.4759 for BBB-rule-pass. Plain has lower RMSE than the specialist model in all five partitions, with paired cluster-bootstrap MSE intervals excluding zero in four. The subgroup therefore does not demonstrate a benefit from specialist training. Its smaller absolute error relative to the full test population should be interpreted alongside its different endpoint distribution.

## Results on feature reliance

The corrected ADME configuration has mean RMSE 0.4930 compared with 0.5122 for Plain and performs better in all five partitions; the direction is retained in the exact-identity sensitivity cohort. This is a comparison of the specified configurations rather than a controlled causal attribution to a feature block.

Group permutation importance complements SHAP. All columns within a block are permuted jointly, with five permutations per block and partition. Mean test MSE increases in the ADME model are 0.317 for ADME, 0.036 for fingerprints, 0.029 for physicochemical descriptors and 0.015 for docking. For Baseline, the docking-block increase is 0.061. These diagnostics describe model reliance. Cross-group relationships are disturbed by permutation, and differing model hyperparameters preclude interpreting the comparison as a matched retraining ablation. Neither permutation nor SHAP establishes a causal toxicity mechanism or suppression of biological information.

## Figure 4 replacement caption

Feature contributions in the Baseline, BBB-rule-pass and ADME configurations for split seed 42. Each panel shows the 20 features with the largest mean absolute SHAP contributions in a held-out sample of 512 compounds. Horizontal position is the contribution to predicted −log₁₀ LD₅₀; colour shows the feature value scaled within that feature between its second and 98th sample percentiles and clipped to this range. Vertical jitter separates overlapping points. Feature attribution describes the fitted model and is not a measurement of biological effect. Panels are provided separately at 180 mm width with 11-point feature labels and vector versions.

## Figure 5 caption additions

The 80×80 display grid covers proper-training first–99th percentiles of mean44 (−8.934 to −0.926) and BTox (0.081 to 2.969). Display scaling and quadrant thresholds are estimated on proper training; the maps summarize held-out test observations. The two-dimensional kernel bandwidth is 0.45 for endpoint mean, variance and interval width, and 0.75 for coverage. Cells with effective support below 120 are masked: 38.08% for the first three maps and 9.34% for coverage. This display smoothing is separate from the four-dimensional interval calibration, which uses bandwidth 0.9 and global fallback at effective support below 80. No explicit boundary-bias correction is applied; low-support masks should not be interpreted as such a correction or as biological exclusion boundaries.

## Discussion opening

Computational biological descriptors help address two practical questions: what toxicity the model predicts and how its errors vary across compounds. After correcting descriptor-to-compound associations, the ADME configuration improves point prediction across all five structural partitions. Local calibration also improves the trade-off between interval width and missed outcomes. This benefit persists when applied to the stronger ADME predictor: its mean interval score decreases from 2.3345 to 2.2571, with mean coverage changing from 90.63% to 89.62% and mean width from 1.5248 to 1.4276.

These findings support using point predictions together with empirical uncertainty and local-support information when selecting compounds for further experimental assessment. They do not establish a mechanism of toxicity. The phenotype descriptors are external model predictions, and the docking scores depend on a particular target panel and protocol. The incremental contribution of raw mean44 over phenotype-only localization is approximately 0.95% in mean interval score, with uncertain paired differences in three of five partitions. Specialist BBB training did not improve on Plain for the same subgroup. The useful result is therefore the measured behaviour of the combined prediction and uncertainty workflow, with each component judged against its appropriate comparator.

## Limitations additions

The repeat analyses are retrospective. Original hypotheses, candidate markers and published hyperparameters were developed using this collection, even though new preprocessing, screening and early stopping were restricted to proper training. External predictor training sets were unavailable for overlap checks. Study-level source and year metadata and the original ADMET request manifest were not recovered. Unresolved structural associations remain excluded. Whole-cluster partitions and cluster-bootstrap uncertainty do not establish coverage under a new laboratory, time period or exposure protocol; the local weighted procedure has no claimed distribution-free conditional-coverage guarantee here. External validation with known study provenance would be needed to assess such transfer.

## Conclusion replacement

Computational biological descriptors complement structural QSAR modelling by supporting both LD₅₀ prediction and characterization of heterogeneous errors. Across repeated structural partitions, local calibration improved the trade-off between interval width and missed outcomes, including with the more accurate ADME predictor. The additional contribution of panel-mean docking was modest, and specialist BBB training did not improve on Plain within the same subgroup. These findings support evaluating biological descriptors for both prediction and uncertainty estimation, while keeping their computational interpretation distinct from experimentally established mechanisms. The conclusions apply to the audited mouse intravenous LD₅₀ collection and the evaluated target panel, docking protocol and model configurations.
