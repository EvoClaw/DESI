# Multi-Agent Brainstorming Results
# Phase 1, Step 5d — 2 Rounds

## Round 1 Summary

| Idea | Visionary | Pragmatic | Scout |
|------|-----------|-----------|-------|
| ARCHON (ARG model comparison) | STRONG | VIABLE | WEAK (timing risk: Durbin likely working on this) |
| PRISM (4-proxy SBI) | VIABLE | WEAK (data constraint) | VIABLE |
| FOSSIL (26-pop joint) | VIABLE | VIABLE | VIABLE |
| DEIMOS (TMRCA mixture) | STRONG | WEAK (TMRCA uncertainty) | STRONG (lowest competition) |
| COMET (BGS correction) | WEAK | KILL | KILL |

Round 1 actions:
- COMET: KILLED (unanimous)
- ARCHON + DEIMOS: MERGED into GALAXY (Visionary proposal)
- PRISM: REDUCED to dual-proxy (SFS + f4 only, works with existing data)
- FOSSIL: REDUCED to FOSSIL-5 (5 continental groups, more tractable)
- ARCHON alone: DEPRECATED (Durbin timing risk)

## Round 2 Summary

| Idea | Visionary | Pragmatic | Scout |
|------|-----------|-----------|-------|
| A: GALAXY (ARG + Bayes factors + TMRCA mixture) | STRONG | VIABLE | VIABLE |
| B: FOSSIL-5 (hierarchical joint NPE on 5 groups) | VIABLE | STRONG | VIABLE |
| C: PRISM-Dual (SFS + f4 dual proxy) | VIABLE | VIABLE | WEAK |

## Convergence

**Consensus reached on top 2 directions: GALAXY and FOSSIL-5.**

PRISM-Dual is insufficient alone for Nature Genetics but is viable as a component of either main direction.

---

## Final Ranked Ideas

═══════════════════════════════════════
DIRECTION 1 (Recommended): GALAXY
"Genealogy-Anchored Long-range ancestrY analysis"
═══════════════════════════════════════

Core question:
Use SINGER ARG inference on ~200-400 diverse 1kGP samples (from phased VCF, publicly available) to simultaneously:
(1) Formally test the FitCoal vs. Cobraa controversy via Bayes factor comparison of 4 competing deep-time demographic models (panmictic+bottleneck, panmictic+no-bottleneck, structured+no-bottleneck, structured+bottleneck) using gLike likelihood computation.
(2) Extract genome-wide TMRCA distributions per 100Ka time window and fit nonparametric mixture models to produce a "structure fraction over time" curve — a continuous profile of how much ancestral population structure existed from 100Ka to 2Ma.

Why this is the recommended direction:
- Directly resolves the hottest controversy in human evolutionary genomics
- Produces a NEW type of output (structure fraction curve) that no existing method provides
- Even if Bayes factors are inconclusive, the TMRCA mixture curve is independently publishable
- The combination of parametric (Bayes factors) + nonparametric (TMRCA mixture) gives both "which model wins" and "what did deep structure look like continuously"
- Differentiates from competing groups: Durbin would likely do Bayes factors (ARCHON-style), but the TMRCA mixture curve is conceptually orthogonal

Key novelty:
- First formal Bayes factor comparison between panmictic and structured deep-time demographic models
- The "structure fraction over time" curve is a novel, interpretable summary: 0 = fully panmictic, 1 = fully structured, over 100Ka–2Ma windows
- Can directly compare the 930Ka window: does structure fraction show a transition here? Is there evidence of two populations being "close to merging" in this window?

Expected findings (if successful):
- If structured model wins: 930Ka signal is an artifact of model misspecification (panmixia assumption), consistent with Cobraa
- If panmictic+bottleneck wins: confirms Hu et al. with better statistical rigor
- Either way: the structure fraction curve provides a novel continuous characterization

Compute plan:
- Download phased 1kGP VCF: ~50 samples per continental group × 4 groups = ~200 samples
- Run SINGER: ~2-4 weeks on 8x L20 GPUs
- Run gLike: fast (days)
- TMRCA extraction + mixture model: days
- Simulation validation: 4-6 weeks

Competition: Durbin lab may do Bayes factors, but TMRCA mixture is not their approach. 12-18 month window.
Risk: SINGER ARG accuracy at >1Ma; TMRCA uncertainty

═══════════════════════════════════════
DIRECTION 2 (Alternative): FOSSIL-5
"Hierarchical Neural Paleodemographic Inference using Five Continental Groups"
═══════════════════════════════════════

Core question:
Build a hierarchical coalescent model for 5 continental groups (African, European, East Asian, South Asian, American) with a SHARED deep history (pre-split, >100Ka) and group-specific recent history (0-100Ka). Use neural posterior estimation (NPE) on the joint folded SFS across all 5 groups + cross-group SFS. Infer: (a) shared Ne(t) from 100Ka-2Ma with calibrated uncertainty, (b) group-specific Ne(t) from 0-100Ka, (c) split times between groups.

Key design: directly exploit that the SFS information-theoretic limit (O(1/log s)) is per-population — but the SHARED parameter is constrained by ALL 5 groups simultaneously, dramatically increasing effective sample size for the pre-split window.

Why this is a strong alternative:
- Can use our EXISTING LD-pruned PLINK data (SFS computable directly)
- Lower competition risk than GALAXY
- Cleaner execution path
- NPE with 8x L20 GPUs is ideal for this
- momi3 as simulator is established and well-tested
- Naturally integrates the 930Ka question as a test: does the shared Ne(t) show a bottleneck at 930Ka?

Expected findings:
- A well-characterized shared Ne(t) for the 100Ka-2Ma window
- Either confirmation or refutation of the 930Ka bottleneck via hierarchical SFS
- The FIRST multi-continental joint characterization of deep human ancestral history

Compute plan:
- Joint SFS from existing PLINK data: days
- momi3 simulations for NPE training: 4-8 weeks on 8x L20 GPUs
- NPE training: 2-4 weeks
- Analysis and comparison with FitCoal/Cobraa: 2 weeks

Competition: Lower — Terhorst's PHLASH is single-population; momi3 is available but not used for this hierarchical deep-time focus.
Risk: Model complexity, parameter identifiability for 5-pop hierarchical model

═══════════════════════════════════════
OPTION C: COMBINED (Both as components of one paper)
═══════════════════════════════════════
If resources allow, GALAXY (ARG-based) + FOSSIL-5 (SFS-based) together constitute a SINGLE paper with two independent lines of evidence converging on the same question:
- GALAXY: ARG-based, formal model testing + TMRCA structure profile
- FOSSIL-5: SFS-based, hierarchical joint inference
- Convergent conclusion: whether the 930Ka signal is real is tested by two independent methodologies
This combined approach is the strongest possible publication for Nature Genetics.

---

## Deliberation History

Round 1 ideas that were killed: COMET (unanimous)
Round 1 ideas that were merged: ARCHON + DEIMOS → GALAXY
Round 1 ideas that were reduced: PRISM (4-proxy → dual-proxy); FOSSIL (26-pop → 5-group)
Key debates:
- ARCHON timing risk: resolved by merging with DEIMOS (TMRCA mixture is not what Durbin would do)
- DEIMOS TMRCA uncertainty: resolved by making Bayes factors the primary result
- PRISM data constraint: resolved by reducing to SFS+f4 or treating as supplement to GALAXY/FOSSIL-5

## Baseline Collection (Type H)

Methods to compare against (must include in evaluation):
- FitCoal (Hu et al. 2023): primary baseline for demographic models
- PSMC/MSMC2: standard coalescent methods
- PHLASH (Terhorst 2025): latest state-of-the-art for single-population Ne(t)
- Cobraa (Cousins et al. 2025): structured ancestry alternative
- Simulations with known ground truth: validation on controlled scenarios
