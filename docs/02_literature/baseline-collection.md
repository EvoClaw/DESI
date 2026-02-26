# Baseline & Alternative Methods Collection
# Phase 1 — Type H Project

## Must-Include Baselines (Reviewers Will Definitely Ask)

| Method | Paper | Type | Deep-time capability | Status |
|--------|-------|------|---------------------|--------|
| FitCoal | Hu et al. 2023 (Science) | SFS composite likelihood | ~2Ma, panmictic | GitHub: YunDeng98/bottleneck_demography |
| PSMC | Li & Durbin 2011 (Nature) | Pairwise HMM, single diploid | ~3Ma | Public |
| MSMC2 | Schiffels & Wang 2020 (Nat Methods) | Multiple sequences HMM | ~1-2Ma | Public |
| PHLASH | Terhorst 2025 (Nature Genetics) | Bayesian GPU coalescent | ~2Ma | GitHub: jthlab/phlash |
| Cobraa | Cousins et al. 2025 (Nature Genetics) | Structured coalescent HMM | ~2Ma | GitHub: available |
| SMC++ | Terhorst et al. 2017 (Nat Genetics) | Many unphased genomes | ~500Ka | Public |

## Important Context Methods

| Method | Paper | Notes |
|--------|-------|-------|
| SINGER | Deng et al. 2025 (Nature Genetics) | ARG inference — input for GALAXY |
| gLike | Ephraim et al. 2025 (Nature Genetics) | Genealogical likelihood — used in GALAXY |
| momi3 | Kamm et al. 2023 | Multi-pop SFS — used in FOSSIL-5 |
| msprime | Kelleher et al. | Coalescent simulator for validation |
| tskit/tsinfer | Kelleher et al. | Tree sequence toolkit |

## Pre-Computed Results (from literature)

| Method | Dataset | Key metric | Source |
|--------|---------|-----------|--------|
| FitCoal | 3,154 individuals (including 1kGP) | Ne min = ~1,280 at 930Ka; bottleneck 813-930Ka | Hu et al. 2023 |
| Cobraa | 1kGP data | 2 ancestral populations, split ~1.5Ma, admixture ~300Ka, 80:20 ratio | Cousins et al. 2025 |
| Cousins & Durvasula | SFS-based | No-bottleneck model substantially better fit than bottleneck | Cousins & Durvasula 2025 |
| PHLASH | Simulations | Lower error than FitCoal, MSMC2, SMC++ | Terhorst 2025 |
