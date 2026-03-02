"""
DESI Sub-group TMRCA Analysis
Addresses reviewer concern E: does the AFR-EAS difference reflect
deep ancestral signal or modern within-AFR sub-population structure?

Uses existing /tmp/desi_pchmm_*.npy outputs — no VCF re-processing needed.

Computes per-pair genome-wide TMRCA for:
  - All 5 continental groups (existing) + their mixed-subpop breakdown
  - Single-subpopulation within: YRI, LWK, GWD, MSL, ESN, CHB, JPT, CHS, CEU, GBR
  - Cross-subpopulation within AFR: YRI×LWK, YRI×GWD, etc.

Key question: what fraction of the ~370 ka AFR-EAS difference persists
when controlling for modern within-AFR subpopulation structure?
"""

import numpy as np
import glob
import os
from itertools import combinations
import json

CHROMS = [f"chr{i}" for i in list(range(1, 23))]
NPY_PATTERN = "/tmp/desi_pchmm_{chrom}.npy"

# Subpopulations of interest
AFR_SUBS  = ['YRI', 'LWK', 'GWD', 'MSL', 'ESN']
EAS_SUBS  = ['CHB', 'JPT', 'CHS']
EUR_SUBS  = ['CEU', 'TSI', 'FIN', 'GBR', 'IBS']
SAS_SUBS  = ['GIH', 'PJL', 'BEB', 'STU', 'ITU']
AMR_SUBS  = ['MXL', 'PUR', 'CLM', 'PEL']

GEN_TIME = 28

def collect_pairs(chroms, pop_filter_fn):
    """Collect genome-wide TMRCA for pairs where pop_filter_fn(pop_i, pop_j) is True."""
    all_tmrca = []
    for chrom in chroms:
        path = NPY_PATTERN.format(chrom=chrom)
        if not os.path.exists(path):
            continue
        d = np.load(path, allow_pickle=True).item()
        hap_pop = d['hap_pop']         # shape (n_hap,) population label per haplotype
        pair_rows = d['pair_rows']     # shape (n_pairs,)
        pair_cols = d['pair_cols']     # shape (n_pairs,)
        pair_tmrca = d['pair_tmrca']   # shape (n_pairs,) genome-wide TMRCA in Ka
        
        # Filter to valid (non-NaN) pairs
        valid = ~np.isnan(pair_tmrca)
        rows = pair_rows[valid]
        cols = pair_cols[valid]
        tmrca = pair_tmrca[valid]
        
        # Get population of each haplotype in the pair
        # hap_pop indexed by haplotype index (2*sample_idx for hap0, 2*sample_idx+1 for hap1)
        pop_i = np.array([hap_pop[r] for r in rows])
        pop_j = np.array([hap_pop[c] for c in cols])
        
        mask = np.array([pop_filter_fn(pi, pj) for pi, pj in zip(pop_i, pop_j)])
        all_tmrca.extend(tmrca[mask].tolist())
    
    return np.array(all_tmrca)


def summarize(tmrca_arr, label):
    if len(tmrca_arr) == 0:
        return {"label": label, "n_pairs": 0, "mean_ka": None, "median_ka": None, "se_ka": None}
    mean = np.mean(tmrca_arr)
    se   = np.std(tmrca_arr) / np.sqrt(len(tmrca_arr))
    return {
        "label": label,
        "n_pairs": len(tmrca_arr),
        "mean_ka": round(float(mean), 1),
        "median_ka": round(float(np.median(tmrca_arr)), 1),
        "se_ka": round(float(se), 2),
    }


def within_sub(pop_list):
    """Return filter: both haplotypes from same subpopulation in pop_list."""
    def fn(pi, pj):
        return pi == pj and pi in pop_list
    return fn


def across_sub(pop_list):
    """Return filter: haplotypes from different subpopulations within pop_list."""
    pop_set = set(pop_list)
    def fn(pi, pj):
        return pi != pj and pi in pop_set and pj in pop_set
    return fn


def between_groups(pop_a, pop_b):
    set_a = set(pop_a)
    set_b = set(pop_b)
    def fn(pi, pj):
        return (pi in set_a and pj in set_b) or (pi in set_b and pj in set_a)
    return fn


def per_sub_within(sub):
    """Single subpop within-pairs."""
    def fn(pi, pj):
        return pi == sub and pj == sub
    return fn


results = {}
print("=" * 70)
print("DESI Sub-group TMRCA Analysis")
print("=" * 70)
print()

# ── 1. Original group-level estimates (reproduce published values) ─────────
print("1. Original continental groups (for reference)")
for label, pop_list in [
    ('within-AFR (all sub)',  list(set(AFR_SUBS))),
    ('within-EAS (all sub)',  list(set(EAS_SUBS))),
    ('within-EUR (all sub)',  list(set(EUR_SUBS))),
]:
    tmrca = collect_pairs(CHROMS, lambda pi, pj, pl=pop_list: pi in pl and pj in pl and pi == pj or (pi in pl and pj in pl))
    # Actually: within any of the group means same-or-different subpop
    pass

# Correct approach: within continental = both haplotypes from continental group
for label, pop_list in [
    ('within-AFR (all sub)',  AFR_SUBS),
    ('within-EAS (all sub)',  EAS_SUBS),
    ('within-EUR (all sub)',  EUR_SUBS),
    ('within-SAS (all sub)',  SAS_SUBS),
    ('within-AMR (all sub)',  AMR_SUBS),
]:
    pop_set = set(pop_list)
    tmrca = collect_pairs(CHROMS, lambda pi, pj, ps=pop_set: pi in ps and pj in ps)
    s = summarize(tmrca, label)
    results[label] = s
    print(f"  {label:<35} n={s['n_pairs']:>8,}  mean={s['mean_ka']:>8.1f} ± {s['se_ka']:.1f} Ka")

print()

# ── 2. Same-subpopulation only ────────────────────────────────────────────
print("2. Single-subpopulation pairs (same-sub only)")
for sub in AFR_SUBS + EAS_SUBS + EUR_SUBS[:3]:
    tmrca = collect_pairs(CHROMS, per_sub_within(sub))
    s = summarize(tmrca, f'within-{sub}')
    results[f'within-{sub}'] = s
    print(f"  {s['label']:<30} n={s['n_pairs']:>8,}  mean={s['mean_ka'] if s['mean_ka'] else 'N/A':>8} Ka")

print()

# ── 3. Cross-subpopulation within AFR ─────────────────────────────────────
print("3. Cross-subpopulation pairs within AFR")
afr_cross_tmrca = collect_pairs(CHROMS, across_sub(AFR_SUBS))
s = summarize(afr_cross_tmrca, 'cross-sub within AFR')
results['cross-sub-AFR'] = s
print(f"  {s['label']:<35} n={s['n_pairs']:>8,}  mean={s['mean_ka']:>8.1f} ± {s['se_ka']:.1f} Ka")

print()

# ── 4. Same-sub AFR vs same-sub EAS comparison ────────────────────────────
print("4. KEY COMPARISON: single-subpop AFR vs single-subpop EAS")
# Best matched: largest single AFR sub vs largest single EAS sub
for afr_sub in ['YRI', 'LWK']:
    for eas_sub in ['CHB', 'JPT']:
        ta = collect_pairs(CHROMS, per_sub_within(afr_sub))
        te = collect_pairs(CHROMS, per_sub_within(eas_sub))
        sa = summarize(ta, f'within-{afr_sub}')
        se = summarize(te, f'within-{eas_sub}')
        if sa['mean_ka'] and se['mean_ka']:
            diff = sa['mean_ka'] - se['mean_ka']
            print(f"  within-{afr_sub} ({sa['mean_ka']:.1f} Ka) − within-{eas_sub} ({se['mean_ka']:.1f} Ka) = {diff:.1f} Ka")

print()

# ── 5. Summary: mixed-sub vs single-sub AFR ───────────────────────────────
print("5. Inflation due to AFR within-continent structure")
afr_all_tmrca = collect_pairs(CHROMS, lambda pi, pj: pi in set(AFR_SUBS) and pj in set(AFR_SUBS))
afr_same_sub  = collect_pairs(CHROMS, lambda pi, pj: pi in set(AFR_SUBS) and pj in set(AFR_SUBS) and pi == pj)
afr_cross_sub = collect_pairs(CHROMS, lambda pi, pj: pi in set(AFR_SUBS) and pj in set(AFR_SUBS) and pi != pj)

s_all  = summarize(afr_all_tmrca, 'AFR all (mixed)')
s_same = summarize(afr_same_sub,  'AFR same-sub only')
s_cross= summarize(afr_cross_sub, 'AFR cross-sub only')

results['AFR_all']        = s_all
results['AFR_same_sub']   = s_same
results['AFR_cross_sub']  = s_cross

print(f"  AFR all-pairs (mixed sub):    mean = {s_all['mean_ka']:.1f} Ka")
print(f"  AFR same-sub pairs only:      mean = {s_same['mean_ka']:.1f} Ka")
print(f"  AFR cross-sub pairs only:     mean = {s_cross['mean_ka']:.1f} Ka")
print(f"  Inflation from cross-sub:     {s_all['mean_ka'] - s_same['mean_ka']:.1f} Ka")

eas_all = collect_pairs(CHROMS, lambda pi, pj: pi in set(EAS_SUBS) and pj in set(EAS_SUBS))
eas_same = collect_pairs(CHROMS, lambda pi, pj: pi in set(EAS_SUBS) and pj in set(EAS_SUBS) and pi == pj)
s_eas_all  = summarize(eas_all,  'EAS all (mixed)')
s_eas_same = summarize(eas_same, 'EAS same-sub only')

results['EAS_all']      = s_eas_all
results['EAS_same_sub'] = s_eas_same

print(f"  EAS all-pairs (mixed sub):    mean = {s_eas_all['mean_ka']:.1f} Ka")
print(f"  EAS same-sub pairs only:      mean = {s_eas_same['mean_ka']:.1f} Ka")
print(f"  Inflation from cross-sub:     {s_eas_all['mean_ka'] - s_eas_same['mean_ka']:.1f} Ka")

diff_mixed    = s_all['mean_ka'] - s_eas_all['mean_ka']
diff_same_sub = s_same['mean_ka'] - s_eas_same['mean_ka']

print()
print(f"  AFR-EAS gap (all pairs, mixed sub):   {diff_mixed:.1f} Ka  ← published value")
print(f"  AFR-EAS gap (same-sub pairs only):    {diff_same_sub:.1f} Ka  ← sub-structure controlled")
print(f"  Attribution to modern AFR structure:  {diff_mixed - diff_same_sub:.1f} Ka ({(diff_mixed-diff_same_sub)/diff_mixed*100:.1f}%)")

print()

# ── 6. Save results ────────────────────────────────────────────────────────
out_path = "/home/yanlin/popghistory/results/final/subgroup_tmrca.json"
with open(out_path, 'w') as f:
    json.dump(results, f, indent=2)
print(f"Results saved to {out_path}")
