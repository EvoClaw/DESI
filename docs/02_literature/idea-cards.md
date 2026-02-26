# Five Candidate Research Ideas
# Phase 1, Step 5c

═══════════════════════════════════════
Idea 1: ARCHON — ARG-based Competitive Hypothesis testing Of deep-time Nested demographic models
═══════════════════════════════════════
Core question:   Can we formally resolve the FitCoal vs. Cobraa controversy by computing
                 Bayes factors between four competing deep-time demographic models
                 (panmictic+bottleneck, panmictic+no-bottleneck, structured+no-bottleneck,
                 structured+bottleneck) using SINGER-inferred ARGs from 1kGP + gLike
                 likelihoods, and which model best explains the data?

Novelty source:  Contradiction between Hu et al. 2023 (panmictic bottleneck) and Cousins
                 et al. 2025 (structured ancestry). gLike + SINGER now make full-likelihood
                 model comparison feasible for the first time (Insight 1 + Insight 4).

Why it matters:  Directly resolves the hottest active controversy in human evolutionary
                 genomics. The answer changes our understanding of: (a) origin of
                 Neanderthals/Denisovans, (b) whether humans nearly went extinct 930Ka,
                 (c) which ancestral model should be used as a null for all future
                 population genetics analyses. A clear resolution = high-profile publication.

Feasibility:     medium-high — SINGER + gLike are published open-source tools. Running
                 SINGER on ~200-500 1kGP samples (African + Eurasian) requires ~weeks on
                 8x L20 GPUs. gLike likelihood computation is fast. Main risk: ARG
                 quality at deep timescales (SINGER was validated to ~2Ma).

Risk:            (1) SINGER ARG accuracy degrades for events >1Ma — need validation on
                 simulated data at these timescales. (2) Model comparison may not yield a
                 clear winner if models are nearly non-identifiable from present-day data.
                 (3) gLike may not handle the specific structured model of cobraa natively.

Competition:     The Durbin group (Cousins, Scally, Durbin) is clearly active in this space
                 (both cobraa and the challenge to FitCoal came from them). However they
                 have NOT done the formal likelihood-based model comparison. This is ours
                 to do.

Estimated scope: Nature Genetics / Nature — flagship, if we definitively resolve the controversy

═══════════════════════════════════════
Idea 2: PRISM — multi-PRoxy Inference for deep-time Structure and Migration
═══════════════════════════════════════
Core question:   Can jointly fitting multiple independently-confounded genomic statistics
                 (SFS, cross-population LD decay, relative cross-coalescence rates, f4
                 statistics) via simulation-based inference (SBI) provide substantially
                 better resolution of deep-time events and explicitly correct for confounders
                 that plague single-statistic methods?

Novelty source:  Cross-domain transfer from paleoclimate multi-proxy science (Insight 5);
                 BGS differential sensitivity (Insight 2); SFS theoretical limit bypass
                 (Insight 3). No existing method combines these statistics.

Why it matters:  (1) Each statistic has a different confounder profile — combining them
                 creates internal consistency checks (like paleoclimate multi-proxy). (2)
                 The discrepancy between LD-based and SFS-based Ne is a BGS diagnostic
                 that no current method exploits. (3) Can detect deep structure even without
                 assuming a specific structured model.

Feasibility:     medium — requires (a) computing 4 types of statistics from 1kGP data
                 (feasible), (b) simulating all 4 simultaneously under both panmictic and
                 structured scenarios (computationally intensive but feasible on 8x L20),
                 (c) training a neural posterior estimator. Main challenge: running enough
                 simulations to cover the structured model space.

Risk:            (1) Statistics may be too correlated to provide independent information.
                 (2) Training neural estimators for high-dimensional structured models is
                 hard — many parameters, slow convergence. (3) LD statistics from LD-pruned
                 data are unavailable — would need phased VCF or alternative LD statistics.

Competition:     SBI for demographic inference exists (Stift et al. 2025) but not for
                 deep-time structured models. Multi-proxy concept is novel.

Estimated scope: Nature Genetics / Genome Research

═══════════════════════════════════════
Idea 3: FOSSIL — Full-population cOalescence Structure Inference using aLL 26 Lineages
═══════════════════════════════════════
Core question:   Can jointly modeling all 26 1kGP populations under a hierarchical
                 coalescent (shared deep history + population-specific recent history)
                 using neural posterior estimation provide dramatically improved resolution
                 of pre-split (>100Ka) demographic events compared to any existing
                 one- or two-population method?

Novelty source:  SFS theoretical limit bypass via multi-population sharing (Insight 3 +
                 Insight 6 / seismic tomography analogy). 26 populations constrain pre-split
                 Ne(t) with N=3,202 total individuals — the shared signal is far above the
                 SFS resolution ceiling for any single population.

Why it matters:  (1) Exploits the full statistical power of 1kGP — no existing method does.
                 (2) Pre-split demographic history is constrained by ALL 26 populations,
                 dramatically improving resolution for exactly the 100Ka-2Ma window.
                 (3) Naturally integrates population-specific recent history, giving a
                 complete human population history simultaneously.

Feasibility:     medium — 26-population joint inference is computationally hard. Key: we
                 can decompose the problem: (a) shared deep history = shared parameter
                 estimated from all 26 populations jointly, (b) population-specific
                 recent history = population-specific parameters. Neural posterior
                 estimation scales to this. momi3 (published 2023) is a multi-population
                 SFS method that provides a technical foundation.

Risk:            (1) Model complexity: 26 populations × split times × admixture = huge
                 parameter space. Need strong priors. (2) Computation: simulating the
                 full 26-population tree is expensive. (3) LD-pruned data loses IBD;
                 need to use SFS + f-statistics only.

Competition:     momi3 does multi-population inference but not with this deep-time focus;
                 Relate/tsdate applied to 1kGP but not for deep-time Ne inference across all
                 26 populations jointly.

Estimated scope: Nature Genetics / PLOS Biology

═══════════════════════════════════════
Idea 4: DEIMOS — DEep-time Inference with Mixture of Structure distributions
═══════════════════════════════════════
Core question:   Can the genome-wide distribution of TMRCA (inferred from SINGER ARGs)
                 across the 100Ka-2Ma window be decomposed via mixture modeling to
                 simultaneously estimate: (a) the degree of ancestral population structure
                 over time (structure proportion as a function of time), (b) Ne within
                 each ancestral subpopulation, and (c) the timing and fraction of
                 admixture events — all WITHOUT assuming panmixia or a specific structured model?

Novelty source:  Counterfactual reasoning (Insight 8): what does the DISTRIBUTION of
                 coalescence times (not just mean) tell us? TMRCA bimodality = structure.
                 TMRCA unimodality = panmixia. This directly reads off structure from data.

Why it matters:  This is model-agnostic: it doesn't assume panmixia OR a specific structure.
                 It DISCOVERS structure from the data. If the 930Ka window shows bimodal
                 TMRCA distribution, that's evidence for structured ancestry, not a
                 bottleneck. This is the most elegant resolution to the controversy.

Feasibility:     medium — requires (a) SINGER ARG inference on diverse 1kGP samples
                 (feasible), (b) TMRCA extraction at each genomic window (straightforward
                 from ARG), (c) Gaussian mixture model fitting over time (simple statistics).
                 Main challenge: SINGER ARG uncertainty propagation into the mixture model.

Risk:            (1) TMRCA estimation at >1Ma may have large uncertainty — the ARG
                 posterior could be very diffuse. (2) Mixture modeling assumes specific
                 distributional forms. (3) Bimodality signal at deep times may be weak.

Competition:     No existing method uses TMRCA mixture modeling for panmixia vs. structure
                 testing at deep timescales. Most direct competition: gLike (parametric,
                 not mixture-based).

Estimated scope: Nature Genetics / PLoS Genetics / Genome Biology

═══════════════════════════════════════
Idea 5: COMET — Correction Of deep-time deMography using Evolution of sTatistics
═══════════════════════════════════════
Core question:   Can the observed systematic discrepancy between LD-based and SFS-based
                 Ne estimates at different time windows be used to build a quantitative
                 correction model for BGS contamination in deep-time demographic inference,
                 and does applying this correction resolve the contradiction between FitCoal
                 (bottleneck) and model-without-bottleneck?

Novelty source:  BGS differential sensitivity insight (Insight 2): LD-based Ne estimates
                 are NOT affected by BGS; SFS-based estimates ARE strongly affected.
                 The discrepancy = BGS signal. Exploit this as a diagnostic and correction.

Why it matters:  If the 930Ka bottleneck signal in FitCoal is partly or entirely BGS
                 contamination, correcting for it would change a major published Science
                 finding. This is a methodologically important contribution regardless of
                 the outcome: it either confirms or refutes the bottleneck with a principled
                 correction framework.

Feasibility:     medium — requires phased data for LD statistics (we have LD-pruned SNPs,
                 which is suboptimal; would need phased VCF). However, LD-based Ne at long
                 range (HapNe) can be computed even from LD-pruned data if we un-prune
                 the relevant LD range. Alternative: use f-statistics as a BGS-robust
                 statistic instead of LD.

Risk:            (1) LD-pruned data loses the LD information needed for this approach.
                 (2) The BGS effect at deep time scales (>500Ka) is theoretically expected
                 but not empirically validated. (3) Risk of incremental result ("confirms
                 that BGS matters") rather than a striking new finding.

Competition:     Cousins et al. 2024 already does BGS-correction for PSMC/MSMC but not
                 for deep time (>500Ka) specifically.

Estimated scope: Genome Research / PLOS Genetics (more limited; unless correction produces
                 striking new biology)
