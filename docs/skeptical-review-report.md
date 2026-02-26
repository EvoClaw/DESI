# Skeptical Review Report — Nature Genetics
## DESI / PCHMM Deep-Time Demographic Inference

**Reviewer persona:** Senior reviewer, Nature Genetics, human population genetics & demographic inference  
**Date:** 2025-02-27  
**Status:** Pre-paper review — identify weaknesses and experiment supplements

---

## Executive Summary

The core empirical finding — within-AFR mean TMRCA (1251 Ka) > within-EAS (875 Ka), difference 376 Ka — is **statistically robust** but **not novel**. Africans having deeper coalescence than non-Africans is the expected outcome of Out-of-Africa (OoA) demography and is already well established via heterozygosity, PSMC, and other methods. The paper's claim that "ancient structure is needed" is **undermined** by the fact that a calibrated OoA model (Ne_b=473) can analytically produce the same 376 Ka difference. Several planned analyses were not executed (formal Bayes factors, baselines on same samples, structure fraction curve). **Recommendation:** Major revisions required; several experiment supplements are REQUIRED before submission to Nature Genetics.

---

## PART A — Vulnerability Analysis

### 1. Alternative Explanations for Key Findings

#### 1.1 AFR > EAS TMRCA Difference (Main Finding)

| Alternative | Plausibility | Evidence For/Against |
|-------------|--------------|----------------------|
| **Out-of-Africa alone** | **HIGH** | OoA model with Ne_b=473 predicts diff≈376 Ka analytically. No structure required. |
| **Recombination rate heterogeneity** | **MEDIUM** | K>1 preferred in ALL groups (AFR, EAS, EUR, AMR, SAS). Low-r regions have longer TMRCA; high-r regions shorter. Mean differences could reflect different r landscapes across populations. |
| **Background selection (BGS)** | **MEDIUM** | BGS reduces diversity in genic regions; AFR has more ancestral diversity. If BGS affects EAS more (e.g., different gene density in regions with different ancestry), mean TMRCA could be biased. No BGS masking was performed. |
| **Phasing errors** | **LOW–MEDIUM** | 1kGP phasing quality varies; EAS may have different phasing error rates. Would need to quantify. |
| **Reference genome bias** | **LOW** | GRCh38 is not African; could bias TMRCA estimates. Unlikely to explain 376 Ka difference. |
| **Sample composition** | **LOW** | 240 samples, 15 groups. Within-AFR n=11,100 pairs vs within-EAS n=4,900 — different sampling variance. Jackknife SE=5.1 Ka suggests stability. |

**Critical point:** The OoA alternative is not merely plausible — it is **analytically demonstrated**. P(coalesce in bottleneck)=31.4% × 55 Ka + 68.6% × 1251 Ka ≈ 875 Ka for EAS. The paper cannot claim "structure is needed" without formally rejecting OoA.

#### 1.2 K-Mixture (K>1) in All Groups

| Alternative | Plausibility | Evidence |
|-------------|--------------|----------|
| **Recombination rate heterogeneity** | **HIGH** | Different genomic regions have different r; TMRCA distributions are mixtures of regions with different r. K>1 does NOT imply structure. |
| **Ancestral structure** | **MEDIUM** | Could be structure, but K>1 in ALL groups (including within-AFR, within-EAS) suggests a universal confounder. |
| **Window-size artifact** | **MEDIUM** | 100 kb windows may straddle recombination hotspots; creates bimodality. |

**Conclusion:** K>1 alone cannot be interpreted as evidence for ancestral structure. The planned "structure fraction over time" curve was meant to test this — it was not produced as designed.

---

### 2. Weakly Supported Claims

| Claim | Support Level | Gap |
|-------|---------------|-----|
| "Ancient structure is needed to explain the data" | **WEAK** | OoA with Ne_b=473 explains diff≈376 Ka. No formal model comparison. |
| "Model D (Cobraa-like) best explains observed difference" | **WEAK** | Model D predicts 370 Ka, OoA predicts 464 Ka with Ne_b=200 — but OoA with Ne_b=473 predicts 376 Ka. Parameter tuning was not systematic. |
| "K>1 indicates ancestral structure" | **WEAK** | K>1 in ALL groups; confounded by recombination. |
| "PCHMM is validated for real data" | **MODERATE** | E1 validated on 50 Mb simulations; real data has BGS, phasing, different r. |
| "AFR/EAS ratio is robust" | **STRONG** | E7 shows ratio=1.429 invariant to μ. |
| "Per-chromosome consistency" | **STRONG** | E6: all 22 chr AFR>EAS, Z=73. |

---

### 3. Gaps vs. Original Plan

| Planned | Done? | Impact |
|---------|-------|--------|
| **Formal Bayes factor comparison** (4 models) | **NO** | **FATAL** — Core contribution was model selection. Instead: mean comparison only. |
| **Structure fraction over time curve** (100 Ka windows, 100 Ka–2 Ma) | **NO** | **MAJOR** — K-mixture was done but not as planned (no time-windowed curve). |
| **SINGER ARGs + gLike likelihoods** | **NO** | **MAJOR** — Method pivot to PCHMM; no ARG-based inference. |
| **FitCoal on same samples** | **NO** | **REQUIRED** — Evaluation protocol: must_run. |
| **PSMC on same samples** | **NO** | **REQUIRED** | |
| **MSMC2 on same samples** | **NO** | **REQUIRED** | |
| **Cobraa on same samples** | **NO** | **REQUIRED** (or compare reported) | |
| **PHLASH on same samples** | **NO** | **REQUIRED** | |
| **Prior sensitivity (Bayes factor stability)** | **NO** | N/A — no Bayes factors. |
| **Simulation model recovery rate** | **PARTIAL** | E2 compares means, not model recovery. |

---

### 4. Artifact Possibilities

| Artifact | Risk | Mitigation |
|----------|------|------------|
| **PCHMM bias at >1 Ma** | **MEDIUM** | E1: within 3.5% on 50 Mb. Real data: 3 Gb, different r, BGS. No validation at 1+ Ma on real-like simulations. |
| **100 kb window size** | **MEDIUM** | May create artificial bimodality at recombination boundaries. No ablation (500 kb, 1 Mb). |
| **Poisson emission assumption** | **LOW–MEDIUM** | Real data may have overdispersion (BGS, selection). |
| **Sample size (240)** | **LOW** | 11,100 AFR pairs, 4,900 EAS pairs — adequate for mean estimates. |
| **BGS** | **MEDIUM** | No masking of genic regions. Literature: BGS strongly affects coalescent methods. |
| **Mutation rate** | **LOW** | Ratio invariant; absolute values depend on μ. |

---

### 5. If I Had to REJECT — Primary Reasons

1. **Lack of novelty:** AFR having deeper TMRCA than non-Africans is the expected outcome of OoA. Heterozygosity, PSMC, and numerous studies already show this. The paper does not establish that PCHMM adds information beyond known facts.

2. **Failure to reject OoA:** The calibrated OoA model (Ne_b=473) explains the 376 Ka difference. The paper claims structure is needed but does not formally compare OoA vs. structured models. This is a **logical gap**.

3. **No baseline comparison:** Evaluation protocol explicitly required FitCoal, PSMC, MSMC2, PHLASH, Cobraa on same samples. None were run. A new method paper without comparison to existing methods is not publishable at Nature Genetics.

4. **Method pivot without justification:** Planned: SINGER + gLike (ARG-based). Delivered: PCHMM (pairwise HMM). The pivot is not motivated; readers will ask why ARGs were abandoned.

5. **K-mixture overinterpretation:** K>1 in all groups is presented as evidence for structure, but recombination heterogeneity is an equally plausible explanation. No positive control (simulation where K=1 is ground truth) to show PCHMM/K-mixture can distinguish.

---

### 6. CRITICAL NOVELTY CHECK

| Finding | Novel? | Domain Expert Prediction |
|---------|--------|--------------------------|
| AFR mean TMRCA > EAS mean TMRCA | **NO** | **Expected.** Heterozygosity: π_AFR > π_EAS. PSMC: Ne_AFR(t) > Ne_EAS(t) at all t. TMRCA is a related statistic. |
| AFR/EAS ratio = 1.43 invariant to μ | **MINOR** | Nice robustness result; not paradigm-shifting. |
| K>1 in all groups | **AMBIGUOUS** | Could be structure OR recombination. Without disambiguation, not novel. |
| Model D predicts diff≈376 Ka | **NO** | Model D was parameterized to match Cobraa; circular. |
| Per-chromosome consistency (Z=73) | **SUPPORTING** | Confirms robustness; not a discovery. |
| ACB≈AFR, PEL≈EAS (E8) | **EXPECTED** | Admixed populations reflect ancestry proportions. |

**Verdict:** There is **no finding** that a domain expert would not have predicted. The paper describes a new method (PCHMM) and applies it to 1kGP, but the biological conclusions are either (a) already known, or (b) insufficiently supported (structure vs. OoA).

---

### 7. CRITICAL CAVEAT ON E2 — Does OoA Undermine "Structure Needed"?

**YES.**

The paper's logic: Model A (panmictic+bottleneck) predicts diff≈0 → ruled out. Model D (structured) predicts diff≈370 → matches observed. Therefore structure is needed.

**Flaw:** The OoA model is a **panmictic** model with a bottleneck. With Ne_b=200, it predicts diff≈464 Ka (too large). But with Ne_b=473, it predicts diff≈376 Ka — **identical to observed**. The OoA model has one free parameter (Ne_b). The paper did not:
- Systematically vary Ne_b to find the best-fitting value
- Compare OoA (best-fit Ne_b) vs. Model D via likelihood or Bayes factor
- Establish that Model D fits *better* than OoA

**Conclusion:** The claim "ancient structure is needed" is **not supported**. OoA with a moderate bottleneck (Ne_b≈473) can explain the data. The paper must either:
- (a) Run formal model comparison (Bayes factors) and show Model D >> OoA, or
- (b) Retract the "structure needed" claim and reframe as "PCHMM recovers expected OoA signal; structure is one possible explanation among others."

---

## PART B — Experiment Supplement Recommendations

### Summary Table

| ID | Supplement | Addresses | Criticality | Scope |
|----|-------------|-----------|-------------|-------|
| S1 | Run FitCoal on same 240 samples | Baseline gap | **REQUIRED** | Medium (days) |
| S2 | Run PSMC on same samples (10 diploids/group) | Baseline gap | **REQUIRED** | Medium (days) |
| S3 | Run MSMC2 on same samples | Baseline gap | **REQUIRED** | Medium (days) |
| S4 | Run Cobraa or compare to reported Cobraa | Baseline gap | **REQUIRED** | Medium (days) |
| S5 | OoA model with Ne_b scanned; formal BF vs Model D | E2 caveat | **REQUIRED** | Medium (days) |
| S6 | BGS masking: exclude genic + 100 kb; rerun PCHMM | BGS artifact | **STRONGLY RECOMMENDED** | Medium (days) |
| S7 | Structure fraction over time curve (as planned) | Planned analysis | **STRONGLY RECOMMENDED** | Medium (days) |
| S8 | Recombination-stratified analysis | K-mixture confound | **STRONGLY RECOMMENDED** | Medium (days) |
| S9 | Window size ablation (500 kb, 1 Mb) | Method validation | **RECOMMENDED** | Small (hours) |
| S10 | PCHMM validation at 1+ Ma on BGS-like simulations | Deep-time calibration | **RECOMMENDED** | Small (hours) |

---

### Detailed Recommendations

#### S1–S4: Baseline Comparison (REQUIRED)

**What:** Run FitCoal, PSMC, MSMC2, and Cobraa (or use reported Cobraa with compatibility note) on the **same 240 samples** (or representative subset: e.g., 10 diploids per continental group for PSMC).

**Why:** Nature Genetics reviewers will ask: "How does PCHMM compare to FitCoal/PSMC?" Without this, the paper cannot be evaluated. The evaluation protocol explicitly required these baselines.

**Scope:** Medium (days). FitCoal and PSMC have published pipelines. MSMC2 requires 4–8 haplotypes per group. Cobraa may have compute constraints.

**Deliverable:** Table comparing Ne(t) or TMRCA summaries from each method. Figure: PCHMM vs. FitCoal vs. PSMC trajectories (if comparable).

---

#### S5: OoA vs. Model D Formal Comparison (REQUIRED)

**What:** 
1. Simulate OoA with Ne_b ∈ {200, 300, 400, 473, 500, 600}.
2. For each, compute predicted AFR–EAS diff.
3. Find Ne_b* such that predicted diff = 376 Ka.
4. Run formal model comparison: marginal likelihood (or composite likelihood) for OoA(Ne_b*) vs. Model D on real data (or on simulations matching real data summary).
5. Report Bayes factor or likelihood ratio.

**Why:** The paper claims structure is needed. The only way to support this is to show that Model D fits significantly better than OoA. Currently, OoA can explain the data — the claim is unfounded.

**Scope:** Medium (days). Requires simulation infrastructure (already have msprime) and likelihood computation. If full likelihood is infeasible, use approximate (e.g., summary-statistic likelihood).

**Deliverable:** Bayes factor (or LR) OoA vs. Model D. If OoA wins or is tied, the paper must reframe.

---

#### S6: BGS Masking (STRONGLY RECOMMENDED)

**What:** Mask genic regions + 100 kb flanking (or use published BGS-masked bed). Rerun PCHMM on neutral regions only. Compare AFR/EAS means and diff to full-genome results.

**Why:** BGS strongly affects coalescent-based methods (Cousins et al. 2024). If the 376 Ka difference shrinks or reverses in neutral regions, the finding is an artifact.

**Scope:** Medium (days). Need BGS mask (available from literature). Rerun pipeline on subset of windows.

**Deliverable:** Table: full genome vs. neutral-only AFR, EAS, diff. If robust, strengthens claim. If not, must report and interpret.

---

#### S7: Structure Fraction Over Time Curve (STRONGLY RECOMMENDED)

**What:** Implement the originally planned analysis: for each 100 Ka window from 100 Ka to 2 Ma, fit 2-component Gaussian mixture to TMRCA distribution, extract w_minor(t) (structure fraction). Report with 95% bootstrap CI (200 genome block resamples).

**Why:** This was a primary metric in the evaluation protocol. It directly tests "does the 930 Ka window show structure?" K-mixture was done but not in time bins — the curve is the key deliverable.

**Scope:** Medium (days). Requires time binning logic and mixture fitting per bin. Code structure may already support this.

**Deliverable:** Figure: w_minor(t) from 100 Ka to 2 Ma with CI. Main text: interpretation at 930 Ka.

---

#### S8: Recombination-Stratified Analysis (STRONGLY RECOMMENDED)

**What:** Stratify 100 kb windows by recombination rate (low / medium / high tertiles, from genetic map). Compute mean TMRCA per stratum, per population. Test whether AFR–EAS diff is consistent across strata.

**Why:** If K>1 is driven by recombination heterogeneity, the diff should vary by stratum. If it's robust across strata, structure is more plausible.

**Scope:** Medium (days). Need recombination map (1kGP or HapMap). Stratify windows, rerun or subset existing results.

**Deliverable:** Table: AFR, EAS, diff by r stratum. If diff is constant → supports structure. If diff varies → confounds K-mixture interpretation.

---

#### S9: Window Size Ablation (RECOMMENDED)

**What:** Rerun PCHMM with 500 kb and 1 Mb windows (in addition to 100 kb). Compare AFR, EAS, diff.

**Why:** 100 kb may be too small (noisy) or create boundary artifacts. Method design mentioned 500 kb / 1 Mb / 2 Mb as ablation.

**Scope:** Small (hours). Parameter change + rerun.

**Deliverable:** Table: AFR, EAS, diff at 100 kb, 500 kb, 1 Mb. If robust, strengthens method.

---

#### S10: Deep-Time PCHMM Validation (RECOMMENDED)

**What:** Simulate data with BGS-like features (reduced diversity in some regions) and known TMRCA at 1 Ma, 1.5 Ma. Run PCHMM. Report bias and RMSE at these depths.

**Why:** E1 validated on 50 Mb neutral simulations. Real data has BGS and different scale. Validation at 1+ Ma on more realistic simulations would strengthen confidence.

**Scope:** Small (hours). Extend existing simulation code.

**Deliverable:** Calibration curve: PCHMM vs. true TMRCA at 1 Ma, 1.5 Ma.

---

## WRITING FIXES (No New Experiments)

| Issue | Writing Fix |
|-------|-------------|
| "Structure is needed" | **Remove or weaken** until S5 is done. Use "consistent with" not "demonstrates." |
| K>1 = structure | Add explicit caveat: "K>1 could reflect recombination heterogeneity; see Discussion." |
| Method pivot | Add "Methods" subsection: "We implemented PCHMM instead of SINGER+gLike due to [X]; PCHMM offers [Y]." Justify the pivot. |
| No baselines | Acknowledge in Limitations: "We did not run FitCoal/PSMC on same samples; comparison to literature values only." (But this is weak — S1–S4 are strongly preferred.) |
| Novelty framing | Reframe: "PCHMM provides a fast, window-based alternative to ARG methods; we recover the expected AFR>EAS TMRCA gradient and show it is robust to μ." Do not oversell. |

---

## Final Verdict

| Criterion | Status |
|-----------|--------|
| **Statistical robustness of main finding** | PASS (Z=73, ratio invariant) |
| **Novelty** | FAIL (AFR>EAS expected from OoA) |
| **Claim "structure needed"** | FAIL (OoA explains it) |
| **Baseline comparison** | FAIL (none run) |
| **Planned analyses completed** | FAIL (Bayes factors, structure curve, baselines) |
| **Method justification** | WEAK (pivot not explained) |

**Recommendation:** Do not submit to Nature Genetics in current form. Address S1–S5 as REQUIRED. Consider S6–S8 strongly. Reframe claims to match evidence. If S5 shows OoA fits as well as Model D, the paper becomes a method paper (PCHMM) with a robustness check on known demography — still publishable at a strong journal, but not as a "structure vs. bottleneck" resolution.
