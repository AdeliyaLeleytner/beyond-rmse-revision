# Response to the reviewers

Thank you for the careful reading of our manuscript. The comments helped us distinguish three questions more clearly: how well the models predict LD₅₀, whether computational biological descriptors help characterize prediction errors, and how much biological interpretation those descriptors support. We repeated the main analyses and carried out the additional comparisons described below.

The reproducibility audit also identified a data-assembly issue: a block of ADMET outputs was shifted relative to compound identifiers. We reconstructed the associations from the source tables and checked them against molecular structures. The new primary analyses use 12,584 verified associations, with a stricter sensitivity analysis restricted to 12,193 exact structural matches. The record accounting is described in our response to Reviewer 1, Comment 1.

This correction changes one result: the ADME configuration now improves point prediction across all five structural partitions. Local uncertainty estimation remains useful, including when applied to this more accurate predictor. The revised interpretation therefore needs to recognize both uses of these descriptors. Below, we distinguish completed analyses from the corresponding text and materials prepared for incorporation into the revision.

## Reviewer 1

### Comment 1 Sample accounting and structural partitions

The one-compound discrepancy has a specific explanation. Compound ligand_id 1298 was present in the original test partition but excluded from the analysis reported in Table 2. We could not establish an independent data-quality reason for this exclusion and reinstated it. With the original predictor otherwise unchanged, reinstatement changes RMSE from 0.5006 to 0.5016.

The audit also confirmed that the original Table 2 calibration set was sampled at the molecule level. The saved cluster assignments contain 613 shared clusters between training and calibration, 724 between training and test, and 343 between calibration and test. Compound identifiers themselves do not overlap. This concerns the uncertainty-analysis partition; the global model comparison used a different partition. The original description did not make that distinction sufficiently clear.

We repeated the analyses with whole clusters assigned to proper training, calibration and test. Clustering uses radius-2, 2,048-bit Morgan fingerprints, without chirality for clustering, distance 1−Tanimoto, Butina distance cutoff 0.4 and reordering=False. GroupShuffleSplit allocates 20% of clusters to test and then 20% of the remaining clusters to calibration. At seed 42, the partitions contain 7,910/2,142/2,532 compounds and 4,548/1,138/1,422 clusters, with no cluster intersections. Seeds 43–46 provide additional partitions; distance cutoffs 0.3 and 0.5 were examined separately. These are documented reanalysis settings, not a claim to have recovered the historical cutoff.

We would also clarify the statistical premise. Structural similarity between training and calibration molecules does not, by itself, violate exchangeability. For ordinary split conformal prediction, the relevant requirement concerns calibration and new observations conditional on the fitted predictor. Close analogues can make an assessment optimistic for transfer to new chemical series, which is a separate concern. Whole-cluster evaluation addresses structural transfer; it does not establish exchangeability or automatically prove coverage for this sampling design.

The wider audit detected 1,678 shifted compound–ADMET associations. Of the 12,652 premerge records, one had an invalid ADMET output and 67 lacked independently confirmed structural correspondence. These records were excluded without using LD₅₀ or prediction errors. The resulting 12,584-compound cohort contains 12,193 exact canonical-isomeric-SMILES matches and 391 unique matches to the same uncharged isomeric parent in both source tables. The main result directions also hold in the exact-match sensitivity cohort. Inclusion, exclusion and split manifests have been prepared for the revision.

### Comment 2 Heterogeneous sources and coverage uncertainty

We agree that coverage should not be assumed to transfer unchanged between data sources. However, pooling independent studies does not by itself make split conformal prediction statistically unsound. Observations sampled from the same heterogeneous mixture can be exchangeable; changes in that mixture or unaccounted dependence create a different problem. Our endpoint already fixes species and route—mouse, intravenous administration—which removes some sources of heterogeneity, though laboratory and protocol differences remain.

We calculated uncertainty using 2,000 bootstrap resamples of whole test clusters and evaluated five structural partitions. In the primary partition, local empirical coverage is 89.42%, with a 95% cluster-bootstrap interval of 87.86–90.91%. Across partitions, mean coverage is 89.37%, ranging from 87.47% to 91.41%; mean global coverage is 89.87%. The bootstrap conditions on the fitted model and calibration set. Inclusion of 90% in a bootstrap interval is not treated as proof of validity.

The cited paper by Oliveira, Orenstein, Ramos and Romano (2024) also studies ordinary split conformal prediction under specified non-exchangeable settings; it does not prescribe a replacement algorithm for every mixed-source dataset. We have not established its conditions for this collection. Our conclusions therefore concern empirical coverage in the tested partitions, without claiming source-shift protection or exact conditional coverage in each local region. [Oliveira et al., JMLR, 2024](https://jmlr.org/papers/v25/23-1553.html).

Reliable study, laboratory and publication-year identifiers are not available for individual measurements in the working table. We cannot reconstruct a defensible source-stratified analysis from these records. This is an explicit limitation, and structural bootstrap is not presented as a substitute. Validation on an independently sourced dataset with study metadata is needed to assess transfer between sources.

### Comment 3 Bandwidths and local support

The comment identifies a distinction that needs to be clearer. The intervals are estimated in four dimensions: mean docking score and three separate predicted toxicity descriptors. That calibration uses h=0.9. The value h=0.45 applies to the two-dimensional display of already computed quantities in Figure 5; the coverage display uses h=0.75. The display bandwidth does not set the calibration bandwidth.

We have prepared separate descriptions of these steps and removed Scott's rule as a justification for calibration. The display bandwidth is a visualization setting, not an established optimum. Calibration sensitivity was also evaluated at h=0.6 and 1.2, without selecting a setting by test performance.

The n_eff=80 fallback is a practical minimum-support rule, not a theoretical constant. It is triggered for 31 of 2,532 test compounds in the primary split (1.22%), and for 0.54–2.94% across the five splits. We additionally evaluated thresholds 40 and 120. Primary-split coverage at thresholds 40/80/120 is 89.49/89.42/89.26%, and interval score is 2.308/2.310/2.318. These results describe sensitivity within the tested range; they do not establish an optimal threshold. Falling back to a global interval need not widen an interval that would otherwise be local.

At the separate display-support threshold n_eff<120, 38.08% of grid cells are masked in the mean, variance and width panels and 9.34% in the coverage panel. The latter uses a broader kernel. Masking identifies inadequate local support but does not fully correct boundary bias.

### Comment 4 BTox and dependence among phenotype descriptors

BTox is not the single phenotype input to calibration. The algorithm uses hERG, Respiratory and Neurotoxicity-DI separately, standardized using proper training. Their sum is a compact display coordinate. The algorithm therefore does not require the sum to be a calibrated quantitative toxicity scale.

Equal weights make the display transparent and avoid fitting additional coefficients. They are not estimates of equal biological contributions. We propose replacing “ordinal-like” with “summed predicted-toxicity coordinate” and stating that different descriptor combinations can have the same sum. Replacing this display with a weighted projection or PCA would change the visualization, but would not itself change intervals already computed from the three separate coordinates.

The requested dependence analysis gives VIFs of 1.15–1.37 for the three actual inputs across the training partitions. hERG-10um is not included alongside hERG. In a separate diagnostic, their Spearman correlation is 0.80–0.81, consistent with retaining only one of these related endpoints. Standardization does not eliminate correlation, which is why these diagnostics are reported separately.

We also examined held-out residual variance in ten BTox ranges defined on proper training. Variance is generally lower at higher BTox in the primary partition, but the relationship is not strictly monotonic. This supports investigating heterogeneous errors, not validating BTox as an ordinal biological scale. Monotonicity of the display sum is not required by the local calibration procedure.

### Comment 5 Averaging docking scores across 44 targets

We agree that mean44 should not be interpreted as a summed binding affinity or a quantitative measure of toxic action. Its role here is more limited: it is one coordinate of a molecule's computational profile on a fixed target panel. Calibration uses neighbouring profiles to estimate prediction errors; it does not convert their average score into a binding probability or assume the same pharmacological action at all targets.

We tested the proposed normalization directly. Each target's mean and standard deviation were estimated on proper training, and the resulting standardized scores were averaged. The predictor, partitions and other calibration settings were held fixed. Across five structural partitions, mean interval score is 2.409 with target-wise normalization and 2.397 with the original mean44; lower is better. Coverage is 89.36% and 89.37%, respectively. Thus, normalization did not provide a consistent improvement in this comparison. We retain the original aggregation as the previously specified primary analysis and report normalization as a sensitivity analysis, rather than selecting an aggregation by test performance.

The additional contribution of raw mean44 over phenotype-only localization is modest: about a 0.95% decrease in mean interval score, with paired 95% intervals including zero in three of five partitions. We therefore do not claim that mean44 is necessary or represents a common toxicity mechanism. We propose replacing “MIE-proxy burden” with “panel-mean docking score.” Target-specific biological weights would require independent evidence that is not available from these scores alone.

### Comment 6 The BBB-pass comparison

We performed the requested comparison on exactly the same held-out BBB-pass compounds. Across five partitions, mean RMSE is 0.4628 for the full-training-set Plain model and 0.4759 for the BBB-pass model. Plain has lower RMSE in all five partitions; in four, the paired cluster-bootstrap interval for the MSE difference excludes zero. The comparison table also includes target variance, RMSE divided by target standard deviation, R², Spearman's ρ and CCC.

These results do not support a benefit from training a separate BBB-pass model. Comparing its RMSE with full-test RMSE cannot establish improved skill because test membership and endpoint range differ. Lower R², ρ or CCC alone also cannot attribute every difference exclusively to range compression; the matched comparison provides the relevant evidence.

We retain this analysis as a characterization of model behaviour in a subgroup defined by computational rules. The interpretation is that the filter identifies a subgroup with a different error distribution, while specialist training did not improve on Plain within that subgroup.

### Comment 7 SHAP and feature-group contributions

A shift in SHAP contributions describes a fitted model; it does not establish that ADME descriptors suppress real mechanisms or that modality bias is the cause. We have prepared revised wording that retains the observed attribution pattern without that causal interpretation.

We also performed grouped permutation importance for Baseline and ADME across all five partitions, with five permutations per group. Columns within each group were permuted jointly to preserve within-group relationships. In the ADME model, the mean MSE increase is 0.317 for the ADME block, 0.036 for fingerprints, 0.029 for physicochemical descriptors and 0.015 for docking descriptors. Permuting the docking block in Baseline increases MSE by 0.061. This provides evidence of model reliance that is independent of the SHAP calculation.

Permutation breaks relationships between groups, and the published model configurations have different hyperparameters. We do not describe this as a fully controlled retraining ablation or a causal biological analysis. Feature-set comparisons, permutation sensitivity and SHAP address complementary questions, with these limitations stated explicitly.

Following the data correction, the ADME configuration has mean RMSE 0.4930 versus 0.5122 for Plain and performs better in all five splits. The supported finding is therefore that adding ADME changes feature reliance and improves prediction in the evaluated configurations. It does not show that biologically meaningful information has been obscured.

### Comment 8 Implementation and reproducibility

The recovered optimization code uses 20 Optuna trials and five internal group folds. Search ranges are iterations 100–800, depth 4–8, learning rate 0.01–0.3 on a logarithmic scale, l2_leaf_reg 1–10, random_strength 0–2 and bagging_temperature 0–5. The TPE sampler and CatBoost use seed 42. For the reanalysis, published hyperparameters were frozen; tree count was selected using an internal group holdout confined to proper training. Calibration and test were not used for this selection. The evaluation remains retrospective because the original choices were developed on this collection.

Figure 5 uses an 80×80 grid spanning proper-training first–99th percentiles: −8.934 to −0.926 for mean44 and 0.081 to 2.969 for BTox. Displayed outcomes come from held-out test observations; coordinate scales and quadrant boundaries come from training. Bandwidths, support thresholds and masked fractions are given in response 3. No explicit boundary-bias correction was applied; range restriction and support masking should not be described as such a correction.

We bootstrapped whole proper-training clusters 1,000 times to assess median thresholds. In the primary split, BTox has median 1.7854 and a 95% interval of 1.7462–1.8235; mean44 has median −6.8159 and an interval of −6.8943 to −6.7455. On average, 1.86% of test compounds change quadrant under bootstrap threshold variation; the corresponding range across splits is 1.26–1.86%. These are descriptive sample boundaries, not biological cutoffs.

The “interval < global” proportions now have both Wilson and whole-test-cluster bootstrap intervals. For example, Pheno+/MIE+ contains 777/786 narrower intervals (98.85%): the Wilson interval is 97.84–99.40%, and the cluster-bootstrap interval is 97.73–99.63%. The latter better reflects the clustered evaluation design; full results cover every quadrant and partition.

The original public snapshot is pinned to commit 23b50718f8eedbfb2eb8dc33067dfb0887b159b2. Scripts, environment versions, manifests, predictions and checksums have been prepared for the new analysis. An updated public release and permanent reference still need to accompany the revised submission; the historical commit identifies the original analysis, not the new results.

## Reviewer 2

### Comment 1 Introduction and objectives

Thank you for these specific suggestions. We have prepared replacement introductory text that explains LD₅₀ as an integrated endpoint of acute toxicity under defined experimental conditions. It is useful for benchmarking QSAR because it reflects a whole-organism outcome rather than a single mechanism. Acute-toxicity measurements also inform hazard classification, as illustrated by OECD 423 for oral testing. We state the scope of this dataset—mouse intravenous toxicity—so this context is not mistaken for validation under that guideline, in humans or for other exposure routes. [OECD Test Guideline 423](https://www.oecd.org/en/publications/test-no-423-acute-oral-toxicity-acute-toxic-class-method_9789264071001-en.html).

The accompanying metric explanation distinguishes error magnitude (MAE and RMSE), squared error (MSE), performance relative to endpoint variance (R²), rank ordering (Spearman's ρ), and numerical agreement including bias and scale (CCC). RMSE is the square root of MSE and weights large errors more strongly than MAE. This prepares readers for apparently different messages from these metrics in a restricted subgroup.

Local heterogeneity is introduced through the question of whether the same predictor is equally reliable for compounds in different parts of the dataset. The proposed final paragraph asks whether computational biological descriptors improve point prediction and whether they help characterize that variation in errors, then outlines the feature comparisons, interval analysis and repeated structural evaluation. Abbreviations are defined on first use, including CNS and MIE, and unnecessary abbreviations are avoided. The Introduction sets out the questions rather than opening with their answers.

### Comment 2 Background and interdisciplinary transitions

We have prepared a clearer transition from structure–activity modelling to machine learning, followed by the reason for estimating uncertainty. A low average prediction error does not tell a reader how large the error might be for the next compound.

Conformal prediction is then introduced through calibration errors and the assumptions behind its standard marginal-coverage result. The local procedure is distinguished from that guarantee. Its connection with applicability-domain analysis is explained in terms of observed error patterns and the availability of relevant calibration examples. This makes the motivation accessible without presenting an interval as a certificate of compound safety.

### Comment 3 Methods

We appreciate this assessment. We retain the overall structure and have prepared additional implementation details on data provenance, structural partitions and calibration, together with the record-level audit and repeat-analysis protocol.

### Comment 4 Results and figures

We have prepared expanded interpretation for Table 1 and Figures 1–2. For example, a difference of 0.5 on the −log₁₀ LD₅₀ scale corresponds to an approximately 3.16-fold dose ratio; this illustrates scale, rather than implying that every prediction has that error. RMSE comparisons are considered alongside MAE, ranking, agreement and paired uncertainty estimates.

The numerical interpretation also reflects the corrected data. ADME now has the best RMSE across five partitions, while the matched BBB comparison does not support specialist-model superiority. We discuss these as empirical differences between evaluated configurations rather than attributing them to an untested biological mechanism.

For Figure 2, the prepared text explains why width must be read together with coverage. Across five partitions of the original confidence-predictor configuration, local calibration reduces mean width from 1.606 to 1.516 while coverage changes from 89.87% to 89.37%. Interval score, which penalizes both width and misses, decreases from 2.478 to 2.397. The conclusion therefore does not rely on narrower intervals alone. Regional variation is discussed in terms of observed residuals and local data support.

Figure 4 has been regenerated as separate panels with readable feature labels and vector versions. The replacement text consistently uses “ADMETlab.” We have followed the reviewer's suggestion to explain each question and observation before presenting its interpretation.

### Comment 5 Discussion

The replacement Discussion opens with the practical finding: computational biological descriptors help characterize heterogeneous prediction errors, and, in the corrected data, direct ADME inclusion also improves point prediction. Implementation details follow that finding.

The practical implication is that prediction intervals and local-support information can help prioritize compounds for additional experimental assessment. We state the dataset's endpoint and exposure limits and distinguish that potential use from a validated regulatory decision procedure. Computational descriptors are not presented as experimentally established toxicity pathways.

### Comment 6 Conclusions

We agree with starting from the positive contribution. The proposed conclusion begins: “Computational biological descriptors complement structural QSAR modelling by supporting both LD₅₀ prediction and characterization of heterogeneous errors. Across repeated structural partitions, local calibration improved the trade-off between interval width and missed outcomes; this benefit also persisted with the more accurate ADME predictor.”

The scope then follows: the additional contribution of mean44 is modest, and specialist BBB training did not outperform Plain on the same subgroup. We therefore cannot adopt the suggested claim of improved BBB-subdomain reliability without that qualification. This retains a positive conclusion while aligning it with the direct comparisons performed during revision.

## Reviewer 3

### Comment 1 Logarithm notation

We propose using log₁₀ consistently in Equation 1 and related labels: y = −log₁₀(LD₅₀ [mol/kg]). Specifying the base removes ambiguity without changing the endpoint transformation.

### Comment 2 Physicochemical descriptor generation

The prepared Methods text distinguishes source-table descriptors, locally calculated descriptors and ADMETlab outputs. TPSA, hydrogen-bond donor and acceptor counts, and rotatable-bond counts are calculated with RDKit CalcTPSA, CalcNumHBD, CalcNumHBA and CalcNumRotatableBonds. For these calculations, a multi-fragment SMILES is represented by its valid fragment with the largest heavy-atom count.

The BBB procedure uses existing master-table MW and logP where available; missing values are replaced by MolWt and Crippen MolLogP calculations. It would therefore be inaccurate to describe all input MW/logP values as newly calculated in RDKit. The global models and original confidence predictor also use different stored molecular-weight columns; that distinction is documented rather than silently changing the input.

The source environment file specifies RDKit 2025.3.5; the reanalysis environment uses 2025.03.6. Stored model fingerprints were independently verified as chirality-aware, radius-2, 2,048-bit Morgan fingerprints for all 12,584 compounds. For upstream MW/logP generation and the external ADMETlab version, we report available provenance without inventing unrecovered settings.

### Comment 3 BBB-pass definition

BBB-pass is a computational physicochemical filter, not an experimental permeability label. The implemented conditions are MW≤500 g/mol, 0≤logP≤5, TPSA≤90 Å², hydrogen-bond donors≤2 and rotatable bonds≤8, all applied together. Acceptor count is calculated but has no separate classification cutoff. Nonfinite MW, logP or TPSA fails the filter.

We checked these rules against the corrected cohort and reproduced the saved classifications for all 12,584 compounds. We propose “BBB-rule-pass subgroup” to avoid equating passage through heuristic rules with measured permeability or the output of a separately validated BBB model.

### Comment 4 Confidence terminology

“Biology-informed uncertainty estimation” is more precise for the task. The output is a prediction interval for LD₅₀, not a probability that a particular prediction is correct. We propose the title “Beyond RMSE: biology-informed uncertainty estimation for toxicological QSAR” and consistent references to prediction intervals and empirical coverage in the text.

### Comment 5 Predicted versus measured phenotypes

These coordinates are outputs of external models; their use does not imply that the endpoints were measured for the compounds in this study. We have prepared consistent wording using “ADMETlab-derived predictions” and “predicted toxicity descriptors.” We also state that external predictor training sets were unavailable for checking possible overlap with this collection.

### Comment 6 Interpretation of docking scores

The proposed wording describes docking scores as computational interaction-related proxies under the selected structures and protocol. They do not establish experimental affinity, activity direction or toxicity mechanism. Biological context comes from the target panel and predicted endpoints; its usefulness for describing prediction errors is evaluated separately from any mechanistic interpretation.

### Comment 7 Aggregation across targets

We performed the suggested target-wise normalization using proper-training means and standard deviations before aggregation. With predictor and partitions unchanged, mean interval score is 2.409 for normalized scores and 2.397 for raw mean44 across five partitions; coverage is 89.36% and 89.37%. Normalization did not provide a consistent improvement and is reported as an aggregation sensitivity analysis.

Mean44 is retained as a simple descriptor of this fixed computational panel, without interpreting it as a combined biological effect. Target-specific biological weights cannot be inferred from docking scores alone and would require independent evidence and a separate selection procedure, rather than tuning against the current test outcomes.

### Comment 8 BTox and alternative representations

The three predicted endpoints are already used separately in calibration. BTox is their display projection. Equal weights make that projection reproducible without fitted coefficients, but do not imply equal biological effects. Substituting PCA or a weighted projection would change the display, not intervals calculated in the original three-coordinate space.

VIFs of 1.15–1.37 across training partitions provide a dependence diagnostic for the actual inputs. They do not validate their sum as a toxicity scale. The proposed interpretation therefore keeps BTox as a visualization coordinate and bases conclusions about interval quality on the separate descriptors and comparisons with global calibration.

### Comment 9 Partitioning and optimization

The reanalysis documents Morgan settings, Butina distance cutoff 0.4, 7,108 clusters and two-stage allocation of whole clusters to training, calibration and test. Seed 42 gives 7,910/2,142/2,532 compounds; seeds 43–46 provide repeats. All cluster intersections were checked and are zero. Cluster separation is not described as an upper bound on every cross-cluster similarity; nearest-neighbour similarity was measured directly.

The historical optimization settings were recovered as 20 trials, five group folds, specified parameter ranges and seed 42. Published hyperparameters were frozen in the repeat analysis, with early stopping confined to an internal proper-training split. Test outcomes were not used to select a new configuration. Because the original hypotheses and configurations were developed on this collection, this is a retrospective reanalysis, not independent external validation. Full settings are provided in responses 1 and 8 to Reviewer 1.

### Comment 10 Scope of conclusions about docking

We agree that the result should be tied to the experiment. We do not draw a general conclusion that docking is unhelpful for toxicology. This study evaluates a particular endpoint, 44-target panel, docking protocol and set of modelling configurations. Within that setting, panel-mean docking makes a modest additional contribution to localization, with uneven statistical support across partitions.

The positive finding is that computational biological context helps characterize heterogeneous errors and adapt prediction intervals. That benefit persists with the more accurate ADME predictor, supporting the joint use of point prediction and uncertainty estimation.

### Comment 11 Mechanistic claims

We distinguish a useful descriptor space from evidence of a mechanism. The descriptors have specified biological referents, but their names do not demonstrate activation of a molecular initiating event or a toxicity pathway. The proposed conclusions concern observed error structure and prediction-interval behaviour. Independent experimental measurements of interactions and downstream effects would be needed for a mechanistic conclusion.
