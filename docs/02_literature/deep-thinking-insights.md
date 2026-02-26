# Deep Thinking Insights
# Phase 1, Step 5b — Six Strategy Analysis

---

## Insight 1
Strategy:    Contradiction Mining
Observation: Hu et al. 2023 (FitCoal, Science) reports a 930Ka bottleneck in a PANMICTIC ancestral population. Cousins et al. 2025 (Cobraa, Nature Genetics) reports NO panmixia — instead finds TWO ancestral populations diverged 1.5 Ma and admixed 300 Ka. Cousins & Durvasula 2025 (Mol Biol Evol) directly challenges FitCoal, showing a panmictic model WITHOUT any bottleneck fits the SFS better than FitCoal's bottleneck model.
Implication: These three papers present mutually incompatible pictures of the same 300Ka-1.5Ma window. The field has NO method that can adjudicate between them. Developing a method that formally tests (via likelihood ratio or Bayes factor) "panmictic + bottleneck" vs. "structured + no bottleneck" would resolve the most heated active controversy in human evolutionary genomics.
Evidence:    Hu et al. 2023 (Science); Cousins et al. 2025 (Nature Genetics cobraa); Cousins & Durvasula 2025 (Mol Biol Evol)
Strength:    strong — backed by three recent high-profile papers in direct disagreement
Verified:    yes

---

## Insight 2
Strategy:    Assumption Challenging
Observation: EVERY SFS-based demographic inference method (FitCoal, PHLASH, MSMC, SMC++) assumes that Ne(t) explains ALL the variation in coalescence rates across the genome. But background selection (BGS) causes local reduction in Ne due to linked purifying selection against deleterious mutations. Crucially, LD-based Ne estimates are VIRTUALLY UNAFFECTED by BGS (Wiuf et al. 2021 PLOS Genetics), while coalescent/SFS-based estimates are STRONGLY affected. Yet no method explicitly uses this differential sensitivity as an internal calibration.
Implication: If LD-based and SFS-based Ne trajectories diverge in the deep-time window (>300Ka), the discrepancy TELLS YOU how much of the "inferred demography" is actually BGS signal. A method that computes both and uses the discrepancy as a BGS correction could significantly change (and improve) deep-time Ne estimates — potentially resolving whether the 930Ka signal is demographic or selective.
Evidence:    Wiuf et al. 2021 (PLOS Genetics); Cousins et al. 2024 (background selection paper); FitCoal paper limitations
Strength:    strong — theoretically grounded and empirically supported by two independent sources
Verified:    yes

---

## Insight 3
Strategy:    Assumption Challenging
Observation: All current deep-time inference methods assume that adding more samples (increasing n) improves inference resolution. Myers et al. 2015 (PNAS) proves mathematically that SFS-based Ne(t) estimation has a minimax error of O(1/log s) independent of n — more individuals do NOT help. The 1kGP has 3,202 individuals and 26 populations, yet current methods only use 2-8 samples per run.
Implication: The path to better deep-time resolution is NOT more individuals from the SAME population — it is using INFORMATION FROM DIFFERENT POPULATIONS that all share the same deep history. Pre-split history is a shared parameter constrained by all populations simultaneously. A method that jointly leverages all 26 populations for the shared deep history could bypass the SFS resolution ceiling.
Evidence:    Myers et al. 2015 (PNAS); standard practice in FitCoal, MSMC uses 2-8 samples
Strength:    strong — mathematically rigorous theoretical grounding
Verified:    yes

---

## Insight 4
Strategy:    Limitation-to-Opportunity Conversion
Observation: FitCoal was published in 2023 and limited to SFS because ARG inference was too slow for large sample sizes. The paper explicitly notes that haplotype/genealogy-based methods would be more informative but infeasible. In September 2025, SINGER was published in Nature Genetics — ARG inference is now 100x faster, scalable to hundreds of whole-genome sequences. Also, gLike (March 2025, Nature Genetics) provides full likelihood computation under parameterized demographic models from genealogical trees.
Implication: The specific computational bottleneck that forced FitCoal to use only SFS is NOW REMOVED. It is now feasible to: (1) infer ARGs from hundreds of 1kGP samples using SINGER, (2) compute full likelihoods under competing demographic models (panmictic bottleneck vs. structured ancestry) using gLike, (3) formally compare models. FitCoal's limitation is now an opportunity.
Evidence:    FitCoal 2023 limitations section; SINGER 2025 (Nature Genetics); gLike 2025 (Nature Genetics)
Strength:    strong — direct limitation-to-tool connection with published tools
Verified:    yes

---

## Insight 5
Strategy:    Cross-Domain Transfer
Observation: Paleoclimate scientists have reconstructed temperature over the past 800,000 years using MULTIPLE PROXIES simultaneously: δ18O isotopes, CO2, methane, dust, sea surface temperature, pollen records. Each proxy is affected differently by confounders (ice volume, ocean circulation, local effects). The power of paleoclimate reconstruction comes from combining correlated but independently confounded signals — inconsistencies between proxies reveal errors.
Implication: Deep-time demographic inference currently uses ONE proxy (SFS or coalescence rate). The "climate science" version of this would use multiple genomic statistics simultaneously: (1) SFS (affected by BGS, invariant to structure), (2) LD decay at different scales (affected differently by recombination, NOT by BGS), (3) cross-population coalescence rates (sensitive to population structure, timing of splits), (4) f3/f4 statistics (sensitive to admixture). Each statistic has different confounders. A multi-proxy method would exploit the complementarity and use inconsistencies to detect methodological artifacts.
Evidence:    Paleoclimate literature; genomic statistics literature; Wiuf et al. 2021 on LD vs SFS; MSMC2 cross-coalescence methods
Strength:    moderate — the analogy is strong but the genomic multi-proxy framework has not been formalized
Verified:    no — needs verification that the four statistics have sufficiently different confounder profiles

---

## Insight 6
Strategy:    Cross-Domain Transfer
Observation: In geophysics, seismic tomography reconstructs the 3D structure of Earth's interior from surface wave measurements. The key insight is: waves traveling through different paths are affected by different structures, so using waves from many angles simultaneously gives better resolution than any single wave. This is an inverse problem solved by combining many imperfect observations with complementary sensitivity.
Implication: Deep-time demographic inference is also an inverse problem (infer Ne(t), split times, admixture from present-day genomes). Different populations are like "waves traveling different paths" — they diverged at different times and carry different information about shared deep history. A tomographic approach: use all 26 1kGP populations as "probes" of the shared ancestral history, each contributing different information, combine via a joint likelihood. This is especially powerful for the pre-split (>100Ka) regime which ALL populations share.
Evidence:    MSMC2 cross-coalescence uses pairs of populations; SMC++ uses multiple populations but assumes panmixia before split; cobraa uses structure but only one structured model
Strength:    moderate — the analogy provides conceptual framing; the technical implementation requires work
Verified:    no — needs verification that the multi-population approach gives significantly more information than current 2-population methods

---

## Insight 7
Strategy:    Trend Extrapolation
Observation: The field has progressed: PSMC (2011, pairs of sequences) → MSMC (2014, 4-8 sequences) → SMC++ (2017, hundreds of unphased sequences) → FitCoal (2023, 3,154 sequences but using only SFS) → PHLASH (2025, Bayesian, GPU). The trend is toward: more samples, better UQ, and GPU acceleration. But NONE of these methods perform JOINT MODEL COMPARISON between fundamentally different model classes (panmictic vs. structured).
Implication: The natural next methodological step is a framework that (a) spans model space (not just parameter space), (b) performs Bayesian model comparison, (c) jointly infers panmixia vs. structure + demographic history, (d) uses modern deep learning for scalability. This is the "D step" that nobody has yet taken.
Evidence:    Publication timeline of PSMC→MSMC→SMC++→FitCoal→PHLASH; absence of model-comparison papers in deep-time inference
Strength:    strong — trend is clear, the gap is clearly identified
Verified:    yes (confirmed no model-comparison paper for panmictic vs. structured exists at this timescale)

---

## Insight 8
Strategy:    Counterfactual Reasoning
Observation: What if we asked "What does the DISTRIBUTION of coalescence times across the genome tell us, rather than just the mean (Ne)?". Current methods model Ne(t) as the mean coalescence rate at each time point. But population STRUCTURE creates a BIMODAL distribution of coalescence times: pairs of lineages within the same subpopulation coalesce faster than pairs from different subpopulations. This bimodality should leave a signature in the shape of the genome-wide distribution of TMRCA (time to most recent common ancestor).
Implication: If we compute the genome-wide distribution of TMRCA at each time window (using SINGER ARGs), and fit a mixture model (bimodal = structure, unimodal = panmixia), we can DIRECTLY TEST panmixia vs. structure without assuming either model. The mixture proportion over time = "degree of structure" as a function of time. This bypasses the FitCoal vs. Cobraa debate entirely by letting the data speak.
Evidence:    SINGER provides ARG posterior; TMRCA distributions computable from ARGs; mixture modeling is standard statistics
Strength:    strong — concrete, testable, directly addresses the controversy
Verified:    no — need to verify that TMRCA bimodality is detectable at deep time scales with current data

---

## Insight 9
Strategy:    Contradiction Mining
Observation: FitCoal uses the FULL SFS (all minor allele frequency categories). Cousins & Durvasula 2025 also uses SFS but finds no bottleneck. The FitCoal paper shows that the bottleneck signal comes primarily from rare variants and the specific shape of the SFS in the low-frequency bins. However, rare variants are (a) most affected by BGS, (b) most affected by recent population bottlenecks at all timescales (not just 930Ka), (c) most affected by sequencing errors and ascertainment bias.
Implication: The 930Ka bottleneck signal in FitCoal may be entirely driven by the low-frequency bins of the SFS, which are the LEAST reliable and MOST confounded bins. A sensitivity analysis that systematically masks different SFS bins and tests whether the bottleneck signal persists would be informative. If the signal disappears when rare variants are excluded, it's likely an artifact.
Evidence:    FitCoal 2023; Cousins & Durvasula 2025; general SFS theory
Strength:    moderate — speculative about which bins drive the signal, but testable
Verified:    no — need to run FitCoal with truncated SFS to verify

---

## Insight 10
Strategy:    Limitation-to-Opportunity Conversion
Observation: The PHLASH paper (Terhorst 2025) explicitly states it improves over FitCoal in uncertainty quantification and deep-time accuracy. However, PHLASH is still limited to single-population (panmictic) history. Its GitHub is public and it's GPU-accelerated.
Implication: PHLASH provides a well-tested, open-source Bayesian foundation. Extending PHLASH to handle two-population histories (with a split + optional admixture at deep times) would be a natural and impactful extension. The original FitCoal's SFS approach would be replaced by PHLASH's more principled Bayesian approach, but within a structured two-population framework. This is more tractable than building a new method from scratch.
Evidence:    PHLASH 2025; Cobraa 2025 (two-population idea); gap in literature
Strength:    moderate — builds on existing tools, feasible, but incremental
Verified:    yes — no multi-population Bayesian extension of PHLASH exists

---

## Insight 11
Strategy:    Assumption Challenging
Observation: ALL current methods assume a CONSTANT MUTATION RATE across time. But mutation rate in hominins has likely changed — generation time has increased, DNA repair mechanisms have evolved, and different studies estimate significantly different per-generation mutation rates. A 20% error in mutation rate translates to a 20% error in TIMING of all inferred events (though not their ordering). The 930Ka event timing has ~10% uncertainty from mutation rate alone.
Implication: A method that JOINTLY INFERS mutation rate variation and demographic history (or that is robust to mutation rate uncertainty) would produce more reliable timing for deep-time events. This is especially critical when comparing inferred genomic events against the fossil/archaeological record. Bayesian inference with a prior on mutation rate, or a sensitivity analysis framework, would be valuable.
Evidence:    Known mutation rate debates in human genetics; Bayesian approaches in molecular clock literature
Strength:    moderate — well-known problem, but combining with demographic inference specifically is underexplored
Verified:    no — need to check if existing methods model mutation rate uncertainty

---

## Discard List (speculative, no verification path)
- Insight about horizontal gene transfer effects: not applicable to humans
- Using environmental data (paleoclimate proxies directly) as a covariate: too many assumptions
