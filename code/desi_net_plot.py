"""
DESI Ne(t) — smooth continuous curve via Gaussian-smoothed histogram.

Method: hazard-function inversion
  Ne(T) = S(T) / (2 * f(T) * GEN_TIME/1000)

Density f(T) is estimated using a fine log-T histogram smoothed with a
Gaussian kernel (scipy.ndimage.gaussian_filter1d).  This gives 300+
continuous time points like PSMC/MSMC output.

Fast block jackknife: pre-bin per block; each jackknife iteration is
O(N_BINS) → all 2700 blocks take < 3 seconds total.
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
import glob, os, warnings
warnings.filterwarnings("ignore")

GEN_TIME = 28
OUT_DIR  = "/home/yanlin/popghistory/results/winhist"
os.makedirs(OUT_DIR, exist_ok=True)

GROUP_NAMES = [
    "within-AFR","within-EAS","within-EUR","within-AMR","within-SAS",
    "AFR-EAS","AFR-EUR","AFR-AMR","AFR-SAS",
    "EAS-EUR","EAS-AMR","EAS-SAS","EUR-AMR","EUR-SAS","AMR-SAS",
]
WITHIN = list(range(5))
COLORS = ["#d62728","#1f77b4","#2ca02c","#ff7f0e","#9467bd"]
LABELS = ["AFR","EAS","EUR","AMR","SAS"]

# Fine log-T grid: 300 bins spanning 100–2500 Ka
N_BINS    = 300
LOG_LO    = np.log(100)
LOG_HI    = np.log(2500)
LOG_EDGES = np.linspace(LOG_LO, LOG_HI, N_BINS + 1)
T_EDGES   = np.exp(LOG_EDGES)
T_MIDS    = np.exp((LOG_EDGES[:-1] + LOG_EDGES[1:]) / 2)
DT        = T_EDGES[1:] - T_EDGES[:-1]   # bin widths in Ka

# Gaussian smoothing sigma in log-bin units
# sigma=12 → smoothing scale ≈ 0.19 log-units ≈ resolves 2-fold changes in time
SMOOTH_SIGMA = 12


def ne_from_logcounts(counts_raw, sigma=SMOOTH_SIGMA):
    """
    Compute Ne(t) from raw histogram counts in log-T space.
    Returns Ne array on T_MIDS grid.
    """
    total = counts_raw.sum()
    if total < 10:
        return np.full(N_BINS, np.nan)

    # Gaussian smooth in log-T space
    counts = gaussian_filter1d(counts_raw.astype(np.float64), sigma=sigma,
                               mode="nearest")
    counts = np.maximum(counts, 0)

    # Density in linear-T space: f(T) = counts / (total * dT)
    f_T = counts / (total * DT)

    # Survival function using smoothed counts
    S = 1.0 - np.cumsum(counts) / counts.sum()
    S = np.clip(S, 0, 1)

    Ne = S / (2.0 * f_T * GEN_TIME / 1000.0)
    # Mask where density is negligible (< 0.5% of peak)
    Ne[f_T < f_T.max() * 0.005] = np.nan
    return Ne


def block_jackknife_ne(t_vals, blks):
    """
    Fast block jackknife for smooth Ne(t).
    Pre-bins all data; each JK iteration subtracts one block's histogram.
    Returns (Ne_all, se_Ne).  Both shape (N_BINS,).
    """
    ok = ~np.isnan(t_vals) & (t_vals > T_EDGES[0]) & (t_vals < T_EDGES[-1])
    v, b = t_vals[ok], blks[ok]

    # Bin all windows once
    bin_idx = np.clip(np.digitize(v, T_EDGES) - 1, 0, N_BINS - 1)
    hist_all = np.bincount(bin_idx, minlength=N_BINS).astype(np.float64)
    Ne_all   = ne_from_logcounts(hist_all)

    ublks = np.unique(b)
    n = len(ublks)

    # Pre-compute per-block histograms
    blk_hists = {bk: np.bincount(bin_idx[b == bk], minlength=N_BINS).astype(np.float64)
                 for bk in ublks}

    # Jackknife
    jk_ne = np.array([ne_from_logcounts(hist_all - blk_hists[bk]) for bk in ublks])

    pseudo = n * Ne_all[None, :] - (n - 1) * jk_ne
    se = np.sqrt(np.nanvar(pseudo, axis=0, ddof=1) / n)
    return Ne_all, se


def load_data():
    files = sorted(glob.glob("/tmp/desi_pchmm_chr*.npy"))
    all_t, all_blk = [], []
    blk_off = 0
    for f in files:
        d = np.load(f, allow_pickle=True).item()
        if "win_t" not in d: continue
        wt = d["win_t"]
        wb = d["win_blk"].astype(np.int64) + blk_off
        blk_off = int(wb.max()) + 1000
        all_t.append(wt); all_blk.append(wb)
    win_t   = np.vstack(all_t)
    win_blk = np.concatenate(all_blk)
    print(f"Loaded: {win_t.shape[0]} windows, {len(np.unique(win_blk))} blocks")
    return win_t, win_blk


def main():
    win_t, win_blk = load_data()

    # ── compute Ne(t) for all within-groups ──────────────────────────────────
    results = {}
    for gi in WITHIN:
        gn = GROUP_NAMES[gi]
        print(f"  Ne(t) jackknife: {gn}...", flush=True)
        Ne, se = block_jackknife_ne(win_t[:, gi], win_blk)
        results[gn] = (Ne, se)

    mids = T_MIDS   # shared x-axis

    # ── Figure 1: all 5 groups — continuous PSMC-style Ne(t) ────────────────
    fig, ax = plt.subplots(figsize=(11, 6))

    for gi, col in zip(WITHIN, COLORS):
        gn = GROUP_NAMES[gi]; lbl = LABELS[gi]
        Ne, se = results[gn]
        ax.fill_between(mids,
                        np.maximum(Ne - 2*se, 300),
                        Ne + 2*se,
                        color=col, alpha=0.15)
        ax.plot(mids, Ne, color=col, lw=2.0, label=lbl)

    ax.axvspan(813, 930, alpha=0.10, color="gray")
    ax.axhline(1280, color="gray", lw=1.2, ls=":", label="FitCoal Ne=1,280")
    ax.axvline(930, color="r",         lw=1.2, ls="--", alpha=0.8)
    ax.axvline(65,  color="steelblue", lw=1.0, ls=":",  alpha=0.7)
    ax.text(930 * 1.04, 400, "930 Ka",   color="r",         fontsize=8, rotation=90, va="bottom")
    ax.text(65  * 1.08, 400, "OOA ~65 Ka", color="steelblue", fontsize=8, rotation=90, va="bottom")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(100, 2500); ax.set_ylim(300, 400_000)
    ax.set_xlabel("Time before present (Ka)", fontsize=12)
    ax.set_ylabel("Effective population size  $N_e$", fontsize=12)
    ax.set_title(f"DESI Ne(t)  —  {N_BINS} time points, Gaussian-smoothed hazard function\n"
                 "(100 kb windows, all pairs aggregated, shaded = ±2 jackknife SE)",
                 fontsize=11)
    ax.legend(fontsize=10, loc="upper left")
    ax.grid(True, which="both", alpha=0.2)

    # Secondary x-axis in years
    ax2 = ax.twiny()
    ax2.set_xscale("log")
    ax2.set_xlim(100e3, 2500e3)
    ax2.set_xlabel("Years before present", fontsize=10, color="gray")
    ax2.tick_params(labelcolor="gray", labelsize=8)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{OUT_DIR}/fig5_net_all_groups.{ext}", dpi=180, bbox_inches="tight")
    plt.close()
    print("  fig5_net_all_groups saved")

    # ── Figure 2: AFR vs EAS side-by-side ────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, gi, col, lbl in zip(axes, [0, 1], ["#d62728", "#1f77b4"], ["AFR", "EAS"]):
        Ne, se = results[GROUP_NAMES[gi]]
        ax.fill_between(mids, np.maximum(Ne - 2*se, 300), Ne + 2*se,
                        color=col, alpha=0.20)
        ax.plot(mids, Ne, color=col, lw=2.5, label=lbl)
        ax.axvspan(813, 930, alpha=0.10, color="gray")
        ax.axhline(1280, color="gray", ls=":", lw=1.2, label="FitCoal Ne=1,280")
        ax.axvline(930, color="r",         ls="--", lw=1.2, alpha=0.8)
        ax.axvline(65,  color="steelblue", ls=":",  lw=1.0, alpha=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(100, 2500); ax.set_ylim(300, 400_000)
        ax.set_xlabel("Time before present (Ka)", fontsize=11)
        ax.set_ylabel("$N_e$", fontsize=11)
        ax.set_title(f"{lbl}:  Ne(t)", fontsize=12, fontweight="bold")
        ax.legend(fontsize=9); ax.grid(True, which="both", alpha=0.2)
    plt.suptitle(f"Ne(t) — AFR vs EAS  (DESI v3, {N_BINS} time points)",
                 fontsize=11, y=1.02)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{OUT_DIR}/fig6_net_afr_eas.{ext}", dpi=180, bbox_inches="tight")
    plt.close()
    print("  fig6_net_afr_eas saved")

    # ── Print table at key time points ───────────────────────────────────────
    checkpoints = [150, 300, 500, 700, 930, 1200, 1500, 2000]
    print(f"\n{'Group':6s}  " + "  ".join(f"{t:>7} Ka" for t in checkpoints))
    for gi in WITHIN:
        Ne, se = results[GROUP_NAMES[gi]]
        vals = []
        for t in checkpoints:
            idx = np.argmin(np.abs(mids - t))
            vals.append(f"{Ne[idx]:>7.0f}" if not np.isnan(Ne[idx]) else "    nan")
        print(f"  {LABELS[gi]:4s}  " + "  ".join(vals))

    print(f"\nDone -> {OUT_DIR}/")


if __name__ == "__main__":
    main()
