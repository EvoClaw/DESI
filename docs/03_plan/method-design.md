# DESI Method Design
# Deep-time Evolutionary Structure Inference
# Phase 3 -- Redesigned: user selected Option C (full standalone new method)

## Core Concept

DESI is a standalone new method for deep-time structured demographic inference.
It does NOT assume panmixia (K=1) like FitCoal/PHLASH, and does NOT assume a
specific K=2 structured model like Cobraa. Instead it infers K from data.

KEY INSIGHT (coalescent theory):
Under a structured ancestry model with K ancestral populations, the genome-wide
distribution of pairwise coalescence times at each time depth t is a K-component
mixture:

    P(T | t) = sum_k  w_k(t) * Gamma(T; alpha_k(t), beta_k(t))

where:
  w_k(t) = fraction of genome coalescing in ancestral population k at time t
            (structure fraction statistic, now formally model-grounded)
  alpha_k(t), beta_k(t) = Gamma shape/rate encoding Ne_k(t):
    E[T | component k] = 2 * Ne_k(t)  [Kingman coalescent]

Under panmixia (K=1): distribution is unimodal.
Under ancestral structure (K>=2): it is bimodal / multi-modal.
DESI infers K and all parameters jointly -- no prior commitment to K.

## Innovation Points (exactly 2)

### Innovation 1: Pairwise-TMRCA Mixture Model for Structured Demographic Inference

No existing method uses the DISTRIBUTION of pairwise T across the genome as the
summary statistic for inferring ancestral population structure and Ne_k(t):
  - FitCoal, PHLASH: assume K=1, use global SFS
  - MSMC2: computes pairwise coalescence rates but averages them (no mixture)
  - Gamma-SMC (Genome Research 2023): infers pairwise T efficiently, no mixture
  - Cobraa: K=2 fixed, uses haplotype HMM not TMRCA distributions

DESI is the first to jointly infer K + Ne_k(t) + w_k(t) from TMRCA distributions.

### Innovation 2: GPU-Accelerated Variational Bayes for Coalescent Mixture Inference

Full Bayesian posterior over (K, w_k(t), Ne_k(t)) via variational inference in JAX.
Exploits 8xL20 GPU parallelism. Outputs calibrated credible intervals.
Scales to thousands of samples.

## Value Proposition vs. Existing Methods

| Method  | Input | K assumption | Output | UQ |
|---------|-------|-------------|--------|----|
| FitCoal | Global SFS | K=1 fixed | Ne(t) | None |
| PSMC    | 1 diploid | K=1 fixed | Ne(t) | None |
| MSMC2   | 4-8 haplotypes | K=1 or K=2 (cross-coal.) | Ne(t), cross-coal. | None |
| PHLASH  | WGS | K=1 fixed | Ne(t) | Posterior CI |
| Cobraa  | Haplotype HMM | K=2 fixed | split+admix params | CI |
| DESI    | Phased VCF | K inferred from data | K, Ne_k(t), w_k(t) | Full posterior |

SFS-based methods cannot distinguish panmictic+bottleneck from structured+admixture
because the global SFS conflates all genomic windows. DESI uses the DISTRIBUTION of
pairwise T values across windows -- sensitive to ancestral structure because
structured ancestry creates bimodal T distributions. Strictly more informative than
the global SFS for detecting structure.

## Full Method Pipeline

### Stage 0 -- Data Preparation

1. Select ~240 phased samples from 1kGP (60 AFR, 50 EUR, 50 EAS, 50 SAS, 30 AMR)
2. Extract autosomes (chr1-chr22) from
   ~/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/
3. Keep biallelic SNVs only; remove multiallelic and indels (bcftools)

### Stage 1 -- Pairwise TMRCA Estimation (fast preprocessing)

For each pair of haplotypes (i,j) and each 1 Mb non-overlapping window w:
  Fast approach: T_ij(w) = het_count(i,j,w) / (2 * mu * L_w)
    mu = 1.2e-8 per bp per generation, L_w = window length
  Accurate approach: run Gamma-SMC per pair for per-site T (used in final analysis)

Parameters:
  n_pairs = 240 * 239 / 2 = 28,560 pairs
  n_windows ~= 3,000 (1 Mb windows, autosomes)
  Total T estimates: ~85 million
  Runtime (het-count): hours on 128-thread CPU (parallelized by chromosome)

### Stage 2 -- Time Binning

Convert T (generations) to calendar time: T_yr = T_gen * 28 yr/gen
Log-scale time bins (9 bins):
  100Ka, 200Ka, 300Ka, 500Ka, 700Ka, 930Ka, 1.25Ma, 1.5Ma, 2.0Ma

For each time bin b:
  Collect all T_ij(w) where T_ij(w) falls within +-30% of bin center
  -> Empirical distribution D_b (hundreds of thousands to millions of values)

### Stage 3 -- DESI Core: K-Component Mixture Inference (NEW ALGORITHM)

For each time bin b, fit K-component Gamma mixture to D_b:

  Generative model:
    z_i | w(b) ~ Categorical(w_1(b), ..., w_K(b))
    T_i | z_i=k ~ Gamma(alpha_k(b), beta_k(b))
    where mean(Gamma) = alpha/beta = 2 * Ne_k(b)  [coalescent]

  Prior:
    K ~ Poisson(2), truncated at K_max = 5
    w(b) ~ Dirichlet(1, ..., 1)
    Ne_k(b) ~ LogNormal(log(10000), 1)  [1K-1M range]

  Inference: Variational Bayes (mean-field VI) via JAX
    ELBO = E_q[log P(data|z,theta)] - KL[q(z,theta) || P(z,theta)]
    Optimization: Adam, 5000 iterations, batch size 50K samples/GPU
    GPU parallelism: 8xL20 -- 8 chains / 8 K-values in parallel
    Output: posterior q*(K, {w_k(b)}, {Ne_k(b)}) for each time bin b

  K selection: run K=1,2,3,4 in parallel -> compare ELBO + penalized BIC
    K_hat = model with highest ELBO minus complexity penalty

### Stage 4 -- Posterior Summary Statistics

From posterior for each time bin b:
  1. K_hat(b) = inferred number of ancestral populations
  2. Ne_k(b) = posterior mean + 95% CI for each component k
  3. w_k(b) = structure fraction -- mixture weights with 95% CI
  4. Trajectory: Ne_major(t), Ne_minor(t), w_minor(t) over time

Key summary curves (main figures):
  - w_minor(t): structure fraction from 100Ka to 2Ma
    if w_minor(930Ka) > 0 with high posterior probability -> evidence vs panmixia
  - Ne_k(t): population-specific effective size trajectories
  - K_hat(t): number of distinct ancestral components per time window

### Stage 5 -- Biological Interpretation

DESI provides direct model comparison internally (no separate gLike step needed):
  K=1 favored at 930Ka -> panmictic bottleneck (FitCoal-consistent)
  K=2 favored at 930Ka with w_minor ~0.20 -> structured ancestry (Cobraa-consistent)
  K=2 with large Ne split -> deep divergence at ~1.5Ma (Cobraa-consistent)

### Stage 6 -- Simulation Validation

msprime simulations under 4 model classes x 50 replicates = 200 simulations:
  Model A: panmictic + bottleneck at 930Ka (FitCoal parameters)
  Model B: panmictic + no bottleneck (smooth Ne(t))
  Model C: two populations, diverge 1.5Ma, admix 300Ka, no within-lineage bottleneck
  Model D: two populations, diverge 1.5Ma, major lineage bottleneck, admix 300Ka (Cobraa)

Per simulation: apply same T estimation pipeline -> run DESI -> compare inferred K
and w_k(t) to ground truth.
Report: model recovery rate (% correct K) per model class, calibration curves.

### Stage 7 -- Robustness and Sensitivity

1. Continental stratification: AFR-only vs. non-AFR
2. Prior sensitivity: vary Ne prior (LogNormal mean +-1 SD)
3. Window size: 500kb vs. 1Mb vs. 2Mb
4. T estimation: het-count vs. Gamma-SMC (accuracy comparison)
5. Mutation rate: mu = 1.0e-8, 1.2e-8, 1.4e-8
6. BGS correction: mask genic regions + 100kb flanking; rerun on neutral regions only

## Ablation Plan

| Component | Ablation | Expected impact | Purpose |
|-----------|----------|-----------------|---------|
| VI vs. EM | EM (point estimate) vs. VI (full posterior) | Wider CI with VI | Show UQ value |
| K inferred vs. fixed K=2 | Force K=2 at every bin | Biased K at older times | Validate data-driven K |
| Window size | 500kb / 1Mb / 2Mb | Power vs. coverage | Validate 1Mb choice |
| T estimation | het-count vs. Gamma-SMC | ~10-20% accuracy gain | Motivate Gamma-SMC |
| Sample size | 50 / 120 / 240 samples | Power curve | Justify N=240 |
| Time bins | 5 / 9 / 15 bins | Resolution trade-off | Validate bin choice |
| BGS masking | with/without genic regions | Effect of selection | Quantify BGS bias |

## Theoretical Proposition (Optional)

Proposition 1: Under any panmictic Ne(t), the TMRCA distribution at any t is unimodal.
Under any K>=2 population structured model, it is K-component with distinct means.
Therefore: K_hat > 1 at any t is a necessary and sufficient condition for structured
ancestry at that time depth (given sufficient power).

Formal proof via coalescent theory: planned for paper, increases publishability.

## Software Tool: desi-pop

Language: Python 3.10+
Core dependencies: JAX (GPU), NumPy, SciPy, msprime (simulation), bcftools (VCF)
Input: phased VCF, sample list, generation time, mutation rate
Output: CSV (K_hat(t), Ne_k(t), w_k(t) with CI), PNG/SVG figures
Release: open-source (MIT), tutorial notebook, documentation

## Predicted Failure Modes and Backup Plans

| Failure Mode | Probability | Backup Plan |
|---|---|---|
| K not identifiable from TMRCA alone | Medium | Add SFS as second statistic (DESI+) |
| het-count T estimates too noisy at >1Ma | Medium | Switch to Gamma-SMC for full analysis |
| VI does not converge | Low-Medium | Use EM algorithm (point estimate) |
| K=1 always wins (no power) | Low | Validate detection threshold on simulations; report power curve |
| All methods agree (no controversy) | Low | Focus on UQ contribution -- DESI provides calibrated CI where FitCoal/Cobraa do not |
