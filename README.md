# DESI — Deep-time Evolutionary Structure Inference

> **⚠️ Automated Research Demo**
> This project and its accompanying paper were generated entirely by **Claude Sonnet 4.6** using the **[Amplify](https://evoclaw.github.io/amplify) agentic research framework** — from literature review and method design through experiment execution and paper writing — with no human editing of scientific content. Results have not been independently validated. Please treat with appropriate caution.

---

## 📄 Paper

**Genome-wide Pairwise TMRCA Excess in African Populations Is Inconsistent with FitCoal's Severe ~930 ka Bottleneck Parameterization**

> The full manuscript is available in [`paper/main.pdf`](paper/main.pdf).

---

## What is DESI?

DESI (**Deep-time Evolutionary Structure Inference**) is a computational method that estimates within-group mean pairwise coalescence times (TMRCA) from windowed aggregate heterozygosity, and uses those estimates to directly test specific ancestral demographic parameterizations against empirical multi-population data.

### The core question

A recent high-profile study (Hu et al. 2023, *Science*) inferred a catastrophic human population bottleneck at ~930,000 years ago, reducing our ancestors to ~1,280 individuals for ~117,000 years. DESI asks: *is this specific parameterization consistent with the observed distribution of within-group genealogical depths across modern continental populations?*

### How DESI works

DESI estimates per-window mean pairwise TMRCA using a **Pairwise Windowed Likelihood (PWL)** — a composite likelihood approach (not an HMM) treating each 100 kb window independently:

1. **Per-window heterozygosity aggregation** — counts heterozygous sites across all haplotype pairs within each population group and 100 kb window
2. **Poisson MLE** — estimates mean pairwise TMRCA per window from aggregate counts under the infinite-sites coalescent
3. **Block jackknife uncertainty** — uses 1 Mb genomic blocks for robust standard errors
4. **Bottleneck parameter scan** — simulates 112 demographic scenarios (varying Ne_bottleneck and duration) using `msprime`/`tskit`, computing the same window-level statistics as the empirical data for direct comparison

### Key findings

Applied to 240 phased whole genomes from the [1000 Genomes Project](https://www.internationalgenome.org/):

| Comparison | Mean TMRCA | Jackknife SE |
|---|---|---|
| Within-AFR (pooled) | 1,229.5 ka | ± 6.4 ka |
| Within-YRI (single sub-population) | 1,228.6 ka | — |
| Within-EAS | 858.8 ka | ± 5.4 ka |
| Within-CHB (single sub-population) | 863.3 ka | — |
| **AFR − EAS** | **370.7 ka** | **± 8.4 ka** |
| **YRI − CHB** | **365.3 ka** | — |

The AFR > EAS asymmetry replicates within single sub-populations (within-YRI vs within-CHB = 365.3 ka), confirming it is not an artifact of pooling multiple African sub-populations.

#### Three-way model comparison

| Model | Predicted AFR−EAS | Direction vs. observed |
|---|---|---|
| Model OoA (no bottleneck, panmictic ancestry + OoA split) | **463.9 ka** | Over-predicts (+25%) |
| FitCoal (Ne=1,280, 117 ka bottleneck) | — | Under-predicts AFR depth by 10.2% |
| Model D (reduced ancestral Ne=12,000) | **370.1 ka** | Matches within 0.2% |
| **Observed** | **370.7 ka** | — |

The observed value sits between the no-bottleneck OoA prediction (too high) and FitCoal's severe-bottleneck prediction (too low), consistent with a genuine but moderate ancestral event.

#### Direct test of FitCoal's parameterization

We simulated FitCoal's exact parameters (Ne=1,280, 117 ka, ancestral Ne=50,000) and measured the same window-level statistic as the data:

| Statistic | FitCoal simulation | Observed | Z | p |
|---|---|---|---|---|
| P(window-mean AFR TMRCA > 930 ka) | **0.615** | **0.707** | 4.6 | <0.0001 |
| Mean within-AFR TMRCA | 1,105 ka | 1,229 ka | — | — |

FitCoal's parameterization systematically under-predicts the depth of African genealogies. A parameter scan of 112 bottleneck scenarios identifies a compatible severity of **Ne ≈ 2,000** — substantially weaker than FitCoal's Ne = 1,280.

#### Robustness
- Signal unchanged in putatively neutral windows >50 kb from genes: AFR−EAS = 361.1 ka (2.6% reduction)
- AFR/EAS TMRCA ratio (1.432) is invariant to mutation rate by construction
- AMR admixed sub-populations replicate the pattern according to known ancestry proportions (within-ACB ≈ within-AFR; within-PEL ≈ within-EAS)

**Conclusion:** FitCoal's specific 930 ka bottleneck parameterization (Ne=1,280 for 117 ka) is quantitatively inconsistent with the multi-population TMRCA distribution. A weaker bottleneck (Ne ≈ 2,000) or partial ancestral population structure provides a better fit.

---

## Repository Structure

```
DESI/
├── paper/
│   ├── main.pdf              ← Full manuscript
│   ├── main.tex              ← LaTeX source
│   ├── sections/             ← Individual section files (abstract, intro, results, discussion, methods)
│   ├── figures/              ← All figures (4 main figures)
│   ├── tables/               ← LaTeX tables
│   └── references.bib        ← Bibliography
├── code/
│   ├── desi_pchmm.py         ← Core PWL implementation
│   ├── desi_run_chr.py       ← Per-chromosome TMRCA estimation
│   ├── desi_winhist_analyze.py ← Genome-wide aggregation & jackknife
│   ├── desi_sim_validate.py  ← Simulation-based calibration
│   ├── desi_subgroup_analysis.py ← Single sub-population analyses (YRI, CHB, etc.)
│   ├── desi_sensitivity.py   ← Mutation rate & AMR analysis
│   ├── desi_bgs_neutral.py   ← Background selection robustness
│   ├── bottleneck_scan.py    ← 112-scenario FitCoal parameter scan
│   ├── plot_bottleneck_scan.py ← Parameter scan heatmaps & figures
│   └── perpair_tmrca.py      ← Per-pair TMRCA extraction from VCF
├── results/
│   ├── final/                ← Full genome TMRCA results & subgroup data
│   ├── winhist/              ← Per-window histogram data & scan figures
│   ├── bottleneck_scan/      ← Parameter scan results (pkl + log)
│   └── perpair/              ← Per-pair TMRCA for chr21/22
└── data/
    └── annotations/          ← UCSC RefGene annotations for BGS masking
```

---

## Method Overview

### Pairwise Windowed Likelihood (PWL)

For each 100 kb genomic window and each population group, DESI counts total heterozygous sites across all haplotype pairs and estimates mean pairwise TMRCA via MLE:

$$\hat{T}_w = \frac{n_w \cdot G}{2\,\mu\,L_w \cdot 1000}$$

where $n_w$ is the aggregate heterozygous count, $L_w$ the callable window length, $\mu = 1.2 \times 10^{-8}$ per bp per generation, and $G = 28$ years per generation.

### Bottleneck Parameter Scan

To directly test FitCoal's parameterization, we simulated 112 combinations of:
- **Ne_bottleneck**: 300, 500, 800, 1,280, 2,000, 4,000, 8,000, 20,000
- **Duration**: 30, 50, 80, 117, 150, 200, 300 ka
- **Ancestral Ne**: 20,000 or 50,000

Each scenario uses the same OoA structure as the empirical data (AFR Ne=15,000, EAS Ne=3,000 post-OoA at 65 ka). Per-window TMRCAs are computed via `tskit` branch-mode diversity, directly comparable to the empirical estimates.

### Calibration

Validated on `msprime` coalescent simulations across three demographic models. PWL recovers exact tree-sequence TMRCA within **<3.5% error** in all cases.

---

## Dependencies

```bash
pip install numpy scipy matplotlib msprime tskit
```

---

## Data

Whole-genome phased VCF data from the [1000 Genomes Project](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) (30× GRCh38, Byrska-Bishop et al. 2022). Data files are not included in this repository due to size.

---

## Key References

- **Hu et al. 2023** — FitCoal bottleneck inference ([*Science* 381:979](https://doi.org/10.1126/science.abq7487))
- **Cousins, Scally & Durbin 2025** — Cobraa structured coalescent ([*Nature Genetics*](https://doi.org/10.1038/s41588-025-02117-1))
- **Cousins & Durvasula 2025** — Critique of FitCoal bottleneck ([*Mol Biol Evol*](https://doi.org/10.1093/molbev/msaf041))
- **Byrska-Bishop et al. 2022** — 1000 Genomes 30× data ([*Cell* 185](https://doi.org/10.1016/j.cell.2022.08.004))

---

## Generated with Amplify

This research project — including literature review, method design, all experiments, and the manuscript — was autonomously conducted by **Claude Sonnet 4.6** using **[Amplify](https://evoclaw.github.io/amplify)**, an open-source agentic framework for automated scientific research.

[![Amplify](https://img.shields.io/badge/Generated%20with-Amplify-blueviolet?style=flat-square)](https://evoclaw.github.io/amplify)

Amplify structures the full research workflow across 7 phases and 4 mandatory human-approval gates, enforcing scientific rigor through discipline skills (metric locking, anti-cherry-picking, claim-evidence alignment, reference verification) and multi-agent deliberation panels.

> **Note:** This repository demonstrates Amplify's capabilities. The scientific content has not been independently peer-reviewed or verified by domain experts. All results should be treated as preliminary and unvalidated.

---

## License

MIT
