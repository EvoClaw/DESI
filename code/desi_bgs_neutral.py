#!/usr/bin/env python3
"""
E_BGS: Neutral Region BGS Robustness Analysis
- Identifies neutral windows (>50kb from any gene)
- Re-runs per-chromosome AFR vs EAS TMRCA comparison on neutral windows only
- Compares neutral vs genome-wide results to test BGS confounding
"""
import numpy as np
import os, glob, sys

# ── Configuration ────────────────────────────────────────────────────────────
WSIZE      = 100_000          # 100 kb window (must match desi_run_chr.py v3)
CHR_FILES  = sorted(glob.glob("/tmp/desi_pchmm_chr*.npy"))
BED_FILE   = "/home/yanlin/popghistory/data/annotations/genes_50kb_buffer.bed"
OUT_DIR    = "/home/yanlin/popghistory/results/validation"
OUT_NPY    = os.path.join(OUT_DIR, "e_bgs_results.npy")

# Population groups (must match desi_run_chr.py GROUP_PAIRS ordering)
GROUP_LABELS = [
    "AFR-AFR", "EAS-EAS", "EUR-EUR", "SAS-SAS", "AMR-AMR",
    "AFR-EAS", "AFR-EUR", "AFR-SAS", "AFR-AMR",
    "EAS-EUR", "EAS-SAS", "EAS-AMR",
    "EUR-SAS", "EUR-AMR", "SAS-AMR",
]
N_GROUPS = len(GROUP_LABELS)
GEN_TIME = 30.0  # years per generation (Henn 2018)

os.makedirs(OUT_DIR, exist_ok=True)

# ── Load gene buffer regions ─────────────────────────────────────────────────
print("Loading gene buffer BED...", flush=True)
genic = {}  # chrom_int -> list of (start_win_idx, end_win_idx)  [half-open]
with open(BED_FILE) as f:
    for line in f:
        chrom, s, e = line.strip().split('\t')
        c = int(chrom[3:])
        s_win = int(s) // WSIZE
        e_win = (int(e) + WSIZE - 1) // WSIZE
        if c not in genic:
            genic[c] = []
        genic[c].append((s_win, e_win))

def is_neutral(chrom_int, win_idx, genic_regions):
    """True if window does NOT overlap any genic region."""
    if chrom_int not in genic_regions:
        return True
    for s, e in genic_regions[chrom_int]:
        if s <= win_idx < e:
            return False
    return True

# Pre-build per-chrom neutral window sets
print("Building neutral window sets...", flush=True)
neutral_wins = {}
# hg38 chrom sizes (approximate, in windows)
CHR_SIZES = {
    1:2490,2:2420,3:1980,4:1900,5:1810,6:1710,7:1590,8:1460,
    9:1380,10:1330,11:1350,12:1330,13:1140,14:1070,15:1020,
    16: 903,17: 833,18: 803,19: 589,20: 644,21: 468,22: 507
}
for c, size in CHR_SIZES.items():
    neutral_wins[c] = set()
    for w in range(size):
        if is_neutral(c, w, genic):
            neutral_wins[c].add(w)
    frac = len(neutral_wins[c]) / size
    print(f"  chr{c:2d}: {len(neutral_wins[c])}/{size} neutral windows ({frac*100:.1f}%)", flush=True)

# ── Load per-chromosome win_t data ───────────────────────────────────────────
print(f"\nFound {len(CHR_FILES)} chromosome files", flush=True)
if len(CHR_FILES) == 0:
    print("ERROR: No chr files found in /tmp/. Run desi_run_chr.py first.", flush=True)
    sys.exit(1)

all_t_all     = {g: [] for g in range(N_GROUPS)}   # all windows
all_t_neutral = {g: [] for g in range(N_GROUPS)}   # neutral windows only
all_blk_all     = {g: [] for g in range(N_GROUPS)}
all_blk_neutral = {g: [] for g in range(N_GROUPS)}

for fpath in CHR_FILES:
    cname = os.path.basename(fpath).replace("desi_pchmm_chr","").replace(".npy","")
    try:
        c_int = int(cname)
    except ValueError:
        print(f"Skipping {fpath}", flush=True)
        continue

    d = np.load(fpath, allow_pickle=True).item()
    win_t = d.get("win_t")
    if win_t is None:
        print(f"  chr{c_int}: no win_t, skipping", flush=True)
        continue

    # win_t shape: (n_windows, N_GROUPS)
    n_wins = win_t.shape[0]
    neut_set = neutral_wins.get(c_int, set())

    for win_idx in range(n_wins):
        for g in range(N_GROUPS):
            t_val = win_t[win_idx, g]
            if np.isnan(t_val) or t_val <= 0:
                continue
            # Block index (1 Mb = 10 windows of 100kb)
            blk_id = (c_int - 1) * 1000 + (win_idx // 10)

            all_t_all[g].append(t_val)
            all_blk_all[g].append(blk_id)

            if win_idx in neut_set:
                all_t_neutral[g].append(t_val)
                all_blk_neutral[g].append(blk_id)

    n_neut = sum(1 for w in range(n_wins) if w in neut_set)
    print(f"  chr{c_int}: {n_wins} total windows, {n_neut} neutral ({n_neut/n_wins*100:.1f}%)", flush=True)

# ── Block jackknife for group means ─────────────────────────────────────────
def block_jackknife_mean(t_vals, blks):
    t_vals = np.array(t_vals, dtype=np.float64)
    blks   = np.array(blks,   dtype=np.int64)
    ok = ~np.isnan(t_vals) & (t_vals > 0)
    t_vals, blks = t_vals[ok], blks[ok]
    if len(t_vals) < 10:
        return np.nan, np.nan, 0

    grand_mean = np.mean(t_vals)
    ublks = np.unique(blks)
    n = len(ublks)
    jk_means = []
    for bk in ublks:
        mask = blks != bk
        if mask.sum() < 2:
            continue
        jk_means.append(np.mean(t_vals[mask]))
    jk_means = np.array(jk_means)
    pseudo = n * grand_mean - (n - 1) * jk_means
    se = np.sqrt(np.nanvar(pseudo, ddof=1) / n)
    return grand_mean, se, len(t_vals)

print("\n=== TMRCA Group Means: All windows vs Neutral windows ===", flush=True)
print(f"{'Group':<15} {'All_Mean':>10} {'All_SE':>8} {'All_N':>8} | {'Neut_Mean':>10} {'Neut_SE':>8} {'Neut_N':>8}", flush=True)
print("-" * 75, flush=True)

results = {}
for g, label in enumerate(GROUP_LABELS):
    mean_all,  se_all,  n_all  = block_jackknife_mean(all_t_all[g],     all_blk_all[g])
    mean_neut, se_neut, n_neut = block_jackknife_mean(all_t_neutral[g], all_blk_neutral[g])
    results[label] = {
        "all_mean": mean_all,  "all_se": se_all,  "all_n": n_all,
        "neut_mean": mean_neut,"neut_se": se_neut,"neut_n": n_neut,
    }
    print(f"{label:<15} {mean_all:>10.1f} {se_all:>8.1f} {n_all:>8d} | {mean_neut:>10.1f} {se_neut:>8.1f} {n_neut:>8d}", flush=True)

# Focus on AFR-AFR vs EAS-EAS difference
afr = results["AFR-AFR"]
eas = results["EAS-EAS"]

diff_all  = afr["all_mean"]  - eas["all_mean"]
diff_neut = afr["neut_mean"] - eas["neut_mean"]
se_diff_all  = np.sqrt(afr["all_se"]**2  + eas["all_se"]**2)
se_diff_neut = np.sqrt(afr["neut_se"]**2 + eas["neut_se"]**2)

print("\n=== BGS Robustness Test ===", flush=True)
print(f"AFR-AFR - EAS-EAS (All windows):     {diff_all:.1f} ± {se_diff_all:.1f} Ka  (Z={diff_all/se_diff_all:.1f})", flush=True)
print(f"AFR-AFR - EAS-EAS (Neutral windows): {diff_neut:.1f} ± {se_diff_neut:.1f} Ka  (Z={diff_neut/se_diff_neut:.1f})", flush=True)
print(f"Change with neutral filter: {abs(diff_neut - diff_all)/diff_all*100:.1f}%", flush=True)

# ── Save results ─────────────────────────────────────────────────────────────
save_dict = {
    "results": results,
    "diff_all": diff_all, "se_diff_all": se_diff_all,
    "diff_neutral": diff_neut, "se_diff_neutral": se_diff_neut,
    "group_labels": GROUP_LABELS,
    "n_neutral_wins": {g: len(all_t_neutral[g]) for g in range(N_GROUPS)},
}
np.save(OUT_NPY, save_dict)
print(f"\nSaved to {OUT_NPY}", flush=True)

# ── Figure ────────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 3, figsize=(16, 6))
fig.suptitle("BGS Robustness: Neutral vs All Windows", fontsize=14, fontweight="bold")

# Panel A: Within-group means comparison
ax = axes[0]
within_groups = ["AFR-AFR", "EAS-EAS", "EUR-EUR", "SAS-SAS", "AMR-AMR"]
x = np.arange(len(within_groups))
w = 0.35

means_all  = [results[g]["all_mean"]  for g in within_groups]
ses_all    = [results[g]["all_se"]    for g in within_groups]
means_neut = [results[g]["neut_mean"] for g in within_groups]
ses_neut   = [results[g]["neut_se"]   for g in within_groups]

b1 = ax.bar(x - w/2, means_all,  w, yerr=ses_all,  label="All windows",     color="#4C72B0", capsize=4, alpha=0.85)
b2 = ax.bar(x + w/2, means_neut, w, yerr=ses_neut, label="Neutral windows", color="#DD8452", capsize=4, alpha=0.85)

ax.set_xticks(x)
ax.set_xticklabels([g.split("-")[0] for g in within_groups], fontsize=11)
ax.set_ylabel("Mean within-group TMRCA (Ka)", fontsize=11)
ax.set_title("A. Within-group TMRCA\n(All vs Neutral windows)", fontsize=11)
ax.legend(fontsize=9)
ax.set_ylim(0, max(means_all)*1.3)

# Panel B: AFR–EAS difference with/without neutral filter
ax = axes[1]
cats  = ["All windows", "Neutral windows"]
diffs = [diff_all, diff_neut]
ses   = [se_diff_all, se_diff_neut]
colors= ["#4C72B0", "#DD8452"]

bars = ax.bar(cats, diffs, yerr=ses, color=colors, capsize=6, width=0.5, alpha=0.85)
ax.axhline(0, color='k', lw=0.8, ls='--')
ax.set_ylabel("AFR − EAS mean TMRCA (Ka)", fontsize=11)
ax.set_title("B. AFR−EAS Difference\n(BGS Robustness)", fontsize=11)
for bar, d, se in zip(bars, diffs, ses):
    ax.text(bar.get_x() + bar.get_width()/2, d + se + 5,
            f"{d:.0f}±{se:.0f}", ha='center', va='bottom', fontsize=10, fontweight='bold')

pct_change = abs(diff_neut - diff_all) / diff_all * 100
ax.text(0.5, 0.05, f"Change: {pct_change:.1f}%\n(BGS effect negligible if <10%)",
        transform=ax.transAxes, ha='center', va='bottom', fontsize=9,
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

# Panel C: Window count breakdown
ax = axes[2]
chrom_data = {}
for fpath in CHR_FILES:
    cname = os.path.basename(fpath).replace("desi_pchmm_chr","").replace(".npy","")
    try:
        c_int = int(cname)
    except ValueError:
        continue
    d = np.load(fpath, allow_pickle=True).item()
    win_t = d.get("win_t")
    if win_t is None:
        continue
    n_wins = win_t.shape[0]
    n_neut = sum(1 for w in range(n_wins) if w in neutral_wins.get(c_int, set()))
    chrom_data[c_int] = (n_wins, n_neut)

chroms = sorted(chrom_data.keys())
frac_neutral = [chrom_data[c][1]/chrom_data[c][0] for c in chroms]
ax.bar([str(c) for c in chroms], [f*100 for f in frac_neutral], color="#55A868", alpha=0.8)
ax.axhline(30.7, color='red', lw=1.5, ls='--', label='Overall 30.7%')
ax.set_xlabel("Chromosome", fontsize=10)
ax.set_ylabel("% Neutral windows (>50kb from gene)", fontsize=10)
ax.set_title("C. Neutral Window Fraction\nper Chromosome", fontsize=11)
ax.tick_params(axis='x', labelsize=7, rotation=45)
ax.legend(fontsize=8)

plt.tight_layout()
out_fig = "/home/yanlin/popghistory/results/winhist/fig_bgs_robustness.png"
plt.savefig(out_fig, dpi=150, bbox_inches='tight')
plt.close()
print(f"Figure saved: {out_fig}", flush=True)
