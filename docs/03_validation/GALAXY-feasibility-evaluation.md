# GALAXY Project: Technical Feasibility Evaluation

**Evaluator:** Senior postdoc in computational population genetics  
**Date:** 2025-02-25  
**Scope:** Research question feasibility, methodology gaps, data sufficiency, timeline, hidden traps, and verdict.

---

## 1. FEASIBILITY

### Can this be answered with available data, tools, and compute?

**Partially.** Core components exist, but several critical mismatches and gaps remain.

| Component | Status | Notes |
|-----------|--------|-------|
| **Data** | ✅ Adequate | Phased 1kGP VCF at `~/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/` — chr1–chr22, all 3202 samples, tabix-indexed. File naming: `1kGP_high_coverage_Illumina.chr{N}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz`. |
| **SINGER** | ⚠️ Mismatch | **SINGER is CPU-only.** No GPU support in the codebase or paper. The 8× L20 GPUs will sit idle for ARG inference. Compute bottleneck is 128-thread CPU. |
| **gLike** | ⚠️ Partial | Accepts `tskit.Tree` objects — compatible with SINGER's `.trees` output. But: (a) requires "selectively neutral and independent" trees; ARG marginal trees are **not** independent (linkage). (b) No built-in 4-model comparison; must construct demography via `Phase` objects. |
| **Bayes factors** | ⚠️ Not built-in | gLike returns log-likelihood. Bayes factors require marginal likelihood (path sampling, stepping-stone, or harmonic mean). gLike does not provide this; you must implement it. |
| **TMRCA mixture model** | ❌ Missing | "Structure fraction over time" from nonparametric mixture fitting to TMRCA distributions per 100Ka window — **no existing tool**. Requires new method development. |

### What's missing?

1. **Tree subsampling strategy** for gLike: How to select "independent" trees from an ARG when neighboring trees are highly correlated? Thinning by physical distance (e.g., every 100 kb) is heuristic; no principled likelihood correction for non-independence.
2. **Bayes factor computation pipeline**: Path sampling or stepping-stone over parameter space for each of 4 models; integration with gLike.
3. **TMRCA mixture model implementation**: Nonparametric mixture (e.g., Dirichlet process) fitted to genome-wide TMRCA distributions in 100Ka windows from 100Ka to 2Ma.
4. **Validation simulations** at 500Ka–2Ma: SINGER paper validated pairwise coalescence accuracy but did not explicitly report performance at >500Ka. Need msprime simulations under known demography to quantify TMRCA bias/variance at deep time.

---

## 2. METHODOLOGY GAP

### Challenges the question glosses over

| Gap | Severity | Description |
|-----|----------|-------------|
| **Tree non-independence** | **High** | gLike's `glike_trees()` documentation: *"it is the user's duty to manually pick out trees that are selective neutral and independent."* ARG marginal trees are correlated at the scale of the recombination rate (~1 cM ≈ 1 Mb). Using all trees inflates effective sample size and invalidates likelihood-based model comparison. Standard approach: thin to every 1–2 Mb. But then you discard most of the ARG; likelihood interpretation is approximate. |
| **Model identifiability** | **High** | At 100Ka–2Ma, the information in present-day genomes is limited. FitCoal's 930Ka bottleneck and Cobraa's 1.5Ma split + 300Ka admixture may be **nearly non-identifiable** from the same data. Cousins & Durvasula 2025 already showed that a panmictic no-bottleneck model fits SFS better than FitCoal's bottleneck — suggesting the data may not strongly favor either. |
| **Cobraa model ≠ gLike structured model** | **Medium** | Cobraa uses an HMM over ancestral segments (2 pops, split 1.5Ma, admixture 300Ka, 80:20). gLike requires a Phase-based demography (split times, admixture proportions, migration matrix). You must **map** Cobraa's inferred structure into a gLike-compatible model. The "structured+bottleneck" and "structured+no-bottleneck" models need explicit parameterization (e.g., split at 1.5Ma, optional bottleneck in one lineage, optional admixture at 300Ka). |
| **SINGER prior vs. demographic models** | **Medium** | SINGER uses a constant-Ne prior (or user-specified Ne) with ARG re-scaling. The inferred ARG is conditioned on this prior. gLike then evaluates likelihood under a *different* demographic model (e.g., bottleneck). There is a **prior mismatch**: the ARG was not sampled from the posterior under the model you are testing. Theoretically, the correct approach is to use the ARG as a sufficient statistic for the genealogy, but the re-scaling step may introduce bias when the true demography differs strongly from the prior. |
| **Deep-time TMRCA accuracy** | **Medium** | SINGER paper: best accuracy among methods for pairwise coalescence; captures bimodality under CEU demography. But validation was at human-relevant timescales (Out-of-Africa, etc.). No explicit validation at 500Ka–2Ma. SMC-based methods can have systematic bias at deep time due to (a) limited mutation information, (b) SMC approximation error. |

---

## 3. DATA SUFFICIENCY

### Is the data adequate?

**Yes for the stated design**, with caveats:

- **Sample size**: 200–400 samples across continental groups is feasible. SINGER paper ran 200 African (chrs 1–4) and 198 British. 400 samples is a ~2× increase; threading cost scales with the number of leaves, so expect roughly 2–4× longer runtime per window.
- **Phasing**: Phased VCF is required. You have it.
- **Population labels**: 1kGP has 26 populations. For "structure fraction over time" and model comparison, you need a clear sampling scheme: e.g., 40 per population × 5 continental groups = 200, or 80 × 5 = 400. Ensure balanced representation to avoid bias toward any one lineage.
- **Coverage**: High-coverage 1kGP. Adequate for ARG inference.

**Potential issue**: The VCF includes SVs. SINGER README: *"non-polymorphic, multi-allelic sites and structural variants are excluded from inference."* So SVs are filtered out automatically — fine.

---

## 4. TIMELINE REALISM

### How long would this take properly?

| Phase | Duration | Notes |
|-------|----------|-------|
| **SINGER ARG inference** | 4–8 weeks | 200–400 samples × 22 chromosomes × 1 Mb windows. SINGER is CPU-only; parallel_singer parallelizes over windows. With 128 threads, estimate ~1–2 weeks per 100 samples for full genome. 400 samples ≈ 4–8 weeks. **No GPU acceleration.** |
| **Validation simulations** | 2–4 weeks | msprime under 4 demographic models; run SINGER on simulated data; compare inferred vs. true TMRCA at 100Ka–2Ma. Essential before real-data analysis. |
| **gLike integration** | 2–3 weeks | Build 4 Phase-based demography models; implement tree subsampling (e.g., 1 Mb spacing); compute log-likelihood per model. |
| **Bayes factor pipeline** | 2–4 weeks | Path sampling or stepping-stone for marginal likelihood; handle optimization and convergence. |
| **TMRCA mixture model** | 4–6 weeks | New method: extract TMRCA per window from ARG; fit nonparametric mixture (e.g., DP mixture) in 100Ka bins; compute "structure fraction" curve. Requires method design, implementation, validation. |
| **Analysis + iteration** | 4–8 weeks | Run on real data; iterate on subsampling, model parameterization; sensitivity analyses. |

**Total: 6–8 months** for a careful execution, assuming no major blockers.

---

## 5. HIDDEN TRAPS

### Practical pitfalls a professor or editor might miss

1. **GPU assumption is wrong**  
   The proposal states "8× L20 GPUs" for SINGER. SINGER uses CPU threading only. **Correct**: 128-thread CPU is the critical resource. GPUs are useful for PHLASH (baseline) or future neural SBI, but not for SINGER.

2. **gLike runtime explosion**  
   gLike: *"number of ancestral populations increase the number of states in an exponential manner."* The structured models (2+ populations) have many more states than panmictic. For 400 samples and deep-time phases, a single gLike evaluation may take minutes to hours. Model comparison across 4 models × 100+ ARG samples × 22 chromosomes × 1 Mb tree spacing = **millions of evaluations**. Need to profile and possibly subsample chromosomes or ARG samples.

3. **Tree selection bias**  
   If you thin trees to every 1 Mb, you get ~3,000 trees per chromosome. But recombination rate varies by region (low in deserts, high near telomeres). Uniform thinning may over-represent high-recombination regions. Consider recombination-aware sampling.

4. **Mutation rate / generation time**  
   SINGER uses generations; you need years for 100Ka–2Ma. Standard: 28 years/generation. 930Ka ≈ 33,000 generations. Mutation rate (1.2e-8/bp/gen) is critical for scaling. Any misspecification propagates to TMRCA interpretation.

5. **BGS and selection**  
   Background selection reduces effective Ne in gene-rich regions. FitCoal's bottleneck signal may be partly BGS. ARG-inferred TMRCAs in gene regions may be biased. Consider masking or down-weighting genic regions.

6. **Cobraa's model is complex**  
   Cobraa infers: 2 pops, split 1.5Ma, bottleneck in major pop, admixture 300Ka (80:20). The "structured" models in your 4-way comparison need to be parameterized to match this. gLike's `twoway_admixture_demo` is a template but has different structure (admixture at t1, not split then admixture). You need a custom Phase sequence.

7. **"First formal resolution"**  
   Cousins & Durvasula 2025 already showed that a panmictic no-bottleneck model fits SFS better than FitCoal's bottleneck. Your contribution is moving from SFS to **full genealogical likelihood** and adding **structured** alternatives. The novelty is the method, not the first ever test.

---

## 6. SUGGESTED IMPROVEMENTS

### How to scope better for feasible execution

1. **Reduce scope for Phase 1**  
   - Start with **200 samples** (e.g., 40 each from 5 populations: YRI, GBR, CHB, PEL, GIH) instead of 400.  
   - Run SINGER on **2–3 chromosomes** (e.g., 1, 10, 22) for validation and pipeline development before full genome.

2. **Validate before scaling**  
   - Run msprime under 4 demographic models (panmictic+bottleneck, etc.).  
   - Infer ARGs with SINGER.  
   - Compute gLike likelihoods and Bayes factors on simulated data.  
   - Verify: (a) correct model is favored when data are simulated under it; (b) TMRCA distributions are recoverable at 500Ka–2Ma.

3. **Defer or simplify "structure fraction" curve**  
   - The novel "structure fraction over time" curve is high-risk (new method, no validation).  
   - Option A: Make it a **secondary** aim; primary = Bayes factor comparison.  
   - Option B: Replace with a simpler summary: e.g., genome-wide TMRCA distribution in 100Ka bins, without mixture fitting. Report bimodality or skew as a descriptive statistic.

4. **Clarify compute resource**  
   - Revise proposal: "128-thread CPU for SINGER; 8× L20 GPUs for PHLASH baseline and potential future SBI extensions."  
   - Avoid implying GPUs accelerate SINGER.

5. **Prioritize 2-way model comparison first**  
   - Before 4-way: panmictic+bottleneck vs. panmictic+no-bottleneck (replicating Cousins & Durvasula with gLike).  
   - Then add structured vs. panmictic.  
   - Reduces complexity and allows incremental validation.

6. **Contact gLike and SINGER developers**  
   - Ask: (a) recommended tree subsampling for ARG; (b) any plans for marginal likelihood / Bayes factors; (c) SINGER runtime for 400 samples genome-wide.

---

## 7. VERDICT

### CONDITIONAL PASS

**Conditions for PASS:**

1. **Validation simulation** (mandatory): Demonstrate on simulated data that (a) gLike + SINGER can recover the correct model among the 4, and (b) TMRCA accuracy at 500Ka–2Ma is acceptable before proceeding to real data.

2. **Correct compute narrative**: Acknowledge that SINGER is CPU-only; do not rely on GPUs for ARG inference.

3. **Scope reduction**: Either (a) start with 200 samples + 2–3 chromosomes, or (b) extend timeline to 8+ months for full genome at 400 samples.

4. **Tree subsampling**: Document and justify the strategy for selecting "independent" trees for gLike (e.g., 1 Mb spacing, recombination-aware sampling).

5. **"Structure fraction" curve**: Either validate on simulations first, or demote to exploratory/secondary aim.

**If conditions are met:** The project is feasible and addresses a significant scientific question. The combination of SINGER + gLike for formal model comparison is novel and could resolve the FitCoal vs. Cobraa controversy.

**If conditions are not met:** High risk of (a) inconclusive results (models non-identifiable), (b) methodological criticism (tree non-independence, prior mismatch), or (c) timeline overrun (SINGER runtime, gLike complexity).

---

## Summary Table

| Criterion | Rating | Key point |
|-----------|--------|-----------|
| Feasibility | Conditional | Data ✓; SINGER CPU-only; gLike partial; TMRCA mixture new |
| Methodology gap | High | Tree non-independence; model identifiability; prior mismatch |
| Data sufficiency | Adequate | 200–400 samples, phased, high quality |
| Timeline | 6–8 months | With validation; 4–8 weeks SINGER alone |
| Hidden traps | Several | GPU myth; gLike runtime; BGS; Cobraa model mapping |
| Verdict | **CONDITIONAL PASS** | Validation + scope + compute narrative required |
