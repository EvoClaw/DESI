# DESI — Deep-time Evolutionary Structure Inference

> **⚠️ Automated Research Demo**
> This project and its accompanying paper were generated entirely by **Claude Sonnet 4.6** using the **[Amplify](https://evoclaw.github.io/amplify) agentic research framework** — from literature review and method design through experiment execution and paper writing — with no human editing of scientific content. Results have not been independently validated. Please treat with appropriate caution.

---

## 📄 Paper

**Ancestral Population Structure Rather Than a Panmictic Bottleneck Explains Deep Differential Genealogical Depth Across Human Populations**

> The full manuscript is available in [`paper/main.pdf`](paper/main.pdf).

---

## What is DESI?

DESI (**Deep-time Evolutionary Structure Inference**) is a computational method for testing whether ancestral human populations were panmictic (single, randomly-mating) or structured (composed of distinct subpopulations) during the Middle Pleistocene (~1.0–0.5 million years ago).

### The core question

A recent high-profile study (Hu et al. 2023, *Science*) inferred a catastrophic human population bottleneck at ~930,000 years ago, reducing our ancestors to ~1,280 individuals. This inference assumed the ancestral population was **panmictic**. DESI asks: *what if it wasn't?*

### How DESI works

Under panmixia, all modern continental populations are random draws from the same ancestral genealogy — their **within-group mean pairwise coalescence times (TMRCA)** must be equal. Ancestral population structure, by contrast, produces persistent genealogical depth differences between groups.

DESI tests this directly:

1. **Per-window TMRCA estimation** — divides each chromosome into 100 kb windows and estimates the mean pairwise TMRCA for each population group using a JAX-accelerated pairwise coalescent HMM (PCHMM) with Poisson likelihood
2. **Genome-wide comparison** — compares within-group mean TMRCAs across five continental super-populations (AFR, EAS, EUR, SAS, AMR)
3. **Block jackknife uncertainty** — uses 1 Mb genomic blocks for robust standard errors that are genome-agnostic

### Key findings

Applied to 240 individuals from the [1000 Genomes Project](https://www.internationalgenome.org/):

| Comparison | Mean TMRCA | Jackknife SE |
|---|---|---|
| Within-AFR | 1,229.5 ka | ± 6.4 ka |
| Within-EAS | 858.8 ka | ± 5.4 ka |
| **AFR − EAS** | **370.7 ka** | **± 8.4 ka** |

- **Any panmictic model predicts this difference to be ~0.2 ka** — the observed difference is ~1,850× larger (Z = 73 across all 22 autosomes)
- A Cobraa-like structured model (two ancestral populations diverging at ~1.5 Ma, 20% admixture at ~300 ka) predicts **370.1 ka** — matching the observation within **0.2%**
- Signal is unchanged in putatively neutral windows (2.6% change), excluding background selection as a confound
- AMR sub-populations replicate the pattern according to their known ancestry composition

**Conclusion:** The ~930 ka signal in human genomic data is better explained by ancestral population structure than by a panmictic bottleneck.

---

## Repository Structure

```
DESI/
├── paper/
│   ├── main.pdf              ← Full manuscript
│   ├── main.tex              ← LaTeX source
│   ├── sections/             ← Individual section files
│   ├── figures/              ← All figures
│   ├── tables/               ← LaTeX tables
│   └── references.bib        ← Bibliography
├── code/
│   ├── desi_pchmm.py         ← Core PCHMM implementation (JAX)
│   ├── desi_run_chr.py       ← Per-chromosome TMRCA estimation
│   ├── desi_winhist_analyze.py ← Genome-wide aggregation & jackknife
│   ├── desi_net_plot.py      ← Ne(t) curve estimation
│   ├── desi_sim_validate.py  ← Simulation-based calibration (E1+E2)
│   ├── desi_e5_proper.py     ← K-mixture analysis (E5)
│   ├── desi_e6_chrjk.py      ← Per-chromosome consistency (E6)
│   ├── desi_sensitivity.py   ← Mutation rate & AMR analysis (E7+E8)
│   └── desi_bgs_neutral.py   ← Background selection robustness
├── results/
│   ├── final/                ← Full genome TMRCA results
│   ├── winhist/              ← Per-window histogram data
│   └── validation/           ← Simulation and robustness results
└── data/
    └── annotations/          ← UCSC RefGene annotations for BGS
```

---

## Method Overview

### Pairwise Coalescent HMM (PCHMM)

For each pair of haplotypes, DESI estimates a genome-wide mean TMRCA using a windowed product likelihood:

- **Window size:** 100 kb
- **Emission model:** Poisson heterozygosity count given a window-specific TMRCA
- **Time bins:** 60 log-spaced bins spanning 10–5,000 ka
- **Group aggregation:** All pairs within a population group are pooled per window for a single group-level MLE

### Calibration

Validated on `msprime` coalescent simulations across three demographic models:
- Standard Out-of-Africa (OoA)
- Panmictic + FitCoal severe bottleneck (Model A)
- Cobraa-like structured model (Model D)

PCHMM recovers exact tree-sequence TMRCA within **<3.5% error** in all cases.

### Uncertainty estimation

Block jackknife over **2,868 non-overlapping 1 Mb genomic blocks** across all 22 autosomes, making the method robust to any single chromosome or genomic region.

---

## Dependencies

```bash
pip install jax jaxlib numpy scipy matplotlib msprime
```

**GPU acceleration:** JAX will automatically use available GPUs (tested on NVIDIA L20).

---

## Data

Whole-genome phased VCF data from the [1000 Genomes Project Phase 3](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38) (30× GRCh38 re-sequencing, Byrska-Bishop et al. 2022).

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
