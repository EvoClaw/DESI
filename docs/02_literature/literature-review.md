# Literature Review — Deep Reading
# Phase 1, Step 4

## Summary Table

| Paper | Access | Problem | Method | Key Results | Critical Limitations |
|-------|--------|---------|--------|-------------|---------------------|
| **Hu et al. 2023** (FitCoal) | abstract_only | Infer deep-time demographic history of ancestral modern humans | Composite likelihood on SFS; FitCoal = fast infinitesimal time coalescent; models Ne(t) as piecewise constant | 930Ka bottleneck, Ne ~1,280 for ~117Ka; 65.85% genetic diversity loss; timing matches MPT, fossil gap, and speciation | Assumes PANMICTIC ancestral population; SFS confounded by BGS; no formal model comparison vs. structured ancestry; no Bayesian UQ |
| **Cousins & Durvasula 2025** | full_text (PMC) | Test whether the 930Ka bottleneck signal is real | SFS-based model comparison; test panmictic bottleneck model vs. panmictic no-bottleneck model | Panmictic model WITHOUT bottleneck fits SFS substantially BETTER than model with bottleneck; insufficient statistical evidence for the severe bottleneck | Only tests panmictic models; does not test structured ancestry (that's cobraa) |
| **Cobraa (Cousins et al. 2025)** | abstract_only | Infer deep ancestral population structure | Structured coalescent HMM; models 2 ancestral populations, their split, and admixture; identifies genomic segments from each ancestry | Two ancestral populations diverged ~1.5 Ma; admixed ~300 Ka in 80:20 ratio; bottleneck in major pop after initial divergence; majority ancestry ancestral to Neanderthals/Denisovans | Does not test panmixia as null; model complexity; requires phased whole genomes |
| **PHLASH (Terhorst 2025)** | abstract_only | Faster, better Bayesian Ne(t) inference | Random projections of coalescent intensity function; score function gradient = same cost as likelihood; GPU-accelerated Python | Lower error than MSMC2, SMC++, FitCoal on simulations; automatic UQ via Bayesian posterior | Still assumes PANMICTIC history; does not test structured ancestry; unclear how it handles BGS |
| **SINGER (Deng et al. 2025)** | abstract_only | Fast ARG inference for hundreds of genomes | MCMC ARG sampling 100x faster than ARGweaver; samples posterior over ARGs | Applied to 1kGP British + African; detects HLA ancient polymorphism, archaic introgression; handles model misspecification better | Not specifically designed for deep-time Ne(t); requires phased WGS |
| **gLike (2025)** | abstract_only | Full-likelihood demographic inference from genealogical trees | Graph structure summarizing lineage relationships; exact likelihood under parameterized models | Outperforms SFS-based methods; infers split times + admixture simultaneously; dozens of parameters | Requires pre-inferred ARG/genealogy; accuracy depends on ARG quality |
| **Deep coalescent hominins (Cousins et al. 2024)** | abstract_only | Push coalescent inference beyond 2 Ma | Extended PSMC/MSMC to deep time using human + chimpanzee data | Ne peak at 5-7 Ma in humans AND chimpanzees but NOT gorilla/orangutan; likely real ancestral signal not artifact | Interpretation uncertain at these time scales; confounding by ancestral recombination rate changes |
| **Background selection (Cousins et al. 2024)** | abstract_only | Account for BGS in demographic inference | Identifies which genomic regions are most affected by BGS; corrects coalescent-based estimates | BGS strongly affects PSMC/MSMC/FitCoal-type methods; LD-based estimates (HapNe-LD, Ne-LD) are NOT affected | Correction is approximate; region identification is uncertain |
| **HapNe (Fournier et al. 2023)** | abstract_only | Infer recent Ne (~2000 years) | IBD + LD; no phasing required; works on ancient DNA | Applied to 1kGP + UK Biobank; detects recent fluctuations | Only reaches ~2,000–10,000 years; NOT applicable to deep time |
| **SFS limits (Myers et al. 2015)** | full_text | Theoretical limits of SFS-based demographic inference | Information-theoretic analysis; minimax rates for SFS-based Ne(t) estimation | Minimax error ≥ O(1/log s) regardless of n (sample size); more individuals do NOT help; this is a fundamental barrier | Fundamental = cannot be overcome by better methods using SFS alone |

## Key Quantitative Facts

- FitCoal bottleneck: Ne ~1,280 (from ~100,000), lasting ~117,000 years (813-930 Ka)
- Cobraa: ancestral divergence ~1.5 Ma (two populations), admixture ~300 Ka
- PHLASH: outperforms FitCoal on simulated data
- Cousins & Durvasula: no-bottleneck model has substantially higher likelihood than bottleneck model
- Deep coalescent: Ne peak ~5-7 Ma, consistent human + chimp, absent in gorilla/orangutan
- LD-based Ne: NOT affected by background selection (Wiuf et al. 2021)
- SFS-based Ne: STRONGLY affected by background selection
- SFS information-theoretic limit: O(1/log s) — adding samples does not help

## Critical Contradiction (The Central Controversy)

**FitCoal story (Hu et al. 2023):** Panmictic ancestral human population → 930Ka severe bottleneck (Ne ~1,280) → population recovery → Out-of-Africa → modern populations

**Cobraa story (Cousins et al. 2025):** NOT panmictic → TWO ancestral populations diverged 1.5 Ma → major pop (80%) ancestral to Neanderthals + Denisovans → minor pop (20%) divergent, slightly deleterious → ~300 Ka admixture → Out-of-Africa → modern populations

**Cousins & Durvasula 2025:** The SFS data itself doesn't support a severe bottleneck even under the panmictic assumption.

→ These are fundamentally different pictures of the same 600K-1.5M year window. The field is at an impasse.
