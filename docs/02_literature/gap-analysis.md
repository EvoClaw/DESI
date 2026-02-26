# Gap Analysis
# Phase 1, Step 5

## Gap 1: No method jointly tests panmixia vs. structured ancestry AND infers demographics
**Confidence: HIGH**

FitCoal assumes panmixia. Cobraa assumes structure. PHLASH assumes panmixia.
No existing method formally tests which ancestral model (panmictic with bottleneck vs. structured with divergence/admixture) better fits the data, nor does any method simultaneously infer both.
Why unsolved: Computationally hard; requires likelihood over model space; structured models have many more parameters.
Approach: ARG-based likelihood (gLike + SINGER) can compute likelihood under BOTH models; Bayes factor comparison.
Status: Not solved. Cousins & Durvasula 2025 does model comparison but only between panmictic models (bottleneck vs. no bottleneck), NOT between panmictic and structured.

## Gap 2: SFS information-theoretic limit blocks deep-time resolution
**Confidence: HIGH**

Myers et al. 2015 proved SFS converges at O(1/log s) — adding more samples does NOT help past a certain point. This means FitCoal, PHLASH, SMC++ all hit the same theoretical ceiling, regardless of sample size.
Why unsolved: The limit is fundamental. Need DIFFERENT DATA TYPES (not more individuals), specifically data types with different information content.
Approach: LD decay encodes Ne at different time scales; cross-coalescence rates encode population separation; combination of statistics from different data types can bypass SFS limits.

## Gap 3: Multi-population information underutilized for pre-split deep history
**Confidence: HIGH**

26 populations in 1kGP all share a common deep history before population splits. All existing methods use 1-4 populations. The pre-split Ne(t) is constrained by ALL populations equally.
Why unsolved: Joint inference across many populations is computationally hard; model complexity explodes.
Approach: The shared deep history creates a natural hierarchical structure. Pre-split component can be inferred jointly; post-split can be inferred population-specifically. Approximate inference (e.g., neural posterior estimation) scales to this.

## Gap 4: Background selection correction not integrated into deep-time inference
**Confidence: HIGH**

LD-based Ne estimates are not affected by BGS, but coalescent/SFS-based estimates are strongly affected (Cousins et al. 2024, Wiuf et al. 2021). No method uses this differential sensitivity as a diagnostic/correction tool.
Why unsolved: Requires parallel computation of LD-based and SFS-based estimates + a model of how BGS affects the discrepancy.
Approach: Compare LD-based and SFS-based Ne trajectories; regions of discrepancy signal BGS contamination; correct SFS-based estimates accordingly.

## Gap 5: No rigorous Bayesian UQ for events beyond 500 Ka
**Confidence: MEDIUM**

PHLASH provides automatic UQ but is still a panmictic method and may not extend well to the 500Ka–2Ma regime. FitCoal provides point estimates. How uncertain is the 930Ka bottleneck timing and depth?
Why unsolved: Bayesian inference in the >500Ka regime requires either: a well-specified prior + likelihood (hard, because model uncertainty is large) or simulation-based inference with massive simulation budgets (now feasible with 8x L20).
Approach: SBI with GPU-accelerated coalescent simulations; train a neural posterior estimator over a broad model class covering both panmictic and structured deep histories.

## Gap 6: The 300Ka-1.5Ma window is poorly characterized across populations
**Confidence: HIGH**

This window encompasses: (a) the 930Ka proposed bottleneck, (b) the Cobraa-inferred admixture at 300Ka, (c) the emergence of Homo sapiens (~300Ka), (d) divergence from Neanderthals/Denisovans (~600-700Ka). Yet all these events are inferred from different analyses using different methods, data, and assumptions.
Why unsolved: Each method can only see part of the picture. No unified analysis simultaneously infers all events in this window.
Approach: A unified framework that infers Ne(t), population splits, admixture events, and their uncertainties simultaneously, across multiple populations.

## Already-Solved Gaps (do NOT pursue)

- Ne inference for recent time (<50Ka): HapNe, TTNE, ancIBD cover this
- Single-population Ne(t) with UQ: PHLASH published Sept 2025
- ARG inference for hundreds of genomes: SINGER published Sept 2025
- Full-likelihood inference from ARG: gLike published March 2025
