"""DESI v3 analysis - aggregate window T estimates, 100 kb windows."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, sys, glob, warnings
warnings.filterwarnings("ignore")

GEN_TIME = 28
MU = 1.2e-8
OUT_DIR = "/home/yanlin/popghistory/results/winhist"
os.makedirs(OUT_DIR, exist_ok=True)

GROUP_NAMES = [
    "within-AFR","within-EAS","within-EUR","within-AMR","within-SAS",
    "AFR-EAS","AFR-EUR","AFR-AMR","AFR-SAS",
    "EAS-EUR","EAS-AMR","EAS-SAS","EUR-AMR","EUR-SAS","AMR-SAS",
]
N_GROUPS = len(GROUP_NAMES)
WITHIN_GROUPS = list(range(5))
GROUP_COLOR = {
    "within-AFR":"#d62728","within-EAS":"#1f77b4","within-EUR":"#2ca02c",
    "within-AMR":"#ff7f0e","within-SAS":"#9467bd",
    "AFR-EAS":"#8c564b","AFR-EUR":"#e377c2","AFR-AMR":"#7f7f7f",
    "AFR-SAS":"#bcbd22","EAS-EUR":"#17becf","EAS-AMR":"#aec7e8",
    "EAS-SAS":"#ffbb78","EUR-AMR":"#98df8a","EUR-SAS":"#ff9896","AMR-SAS":"#c5b0d5",
}


def load_data(chroms=None):
    files = ([f"/tmp/desi_pchmm_{c}.npy" for c in chroms]
             if chroms else sorted(glob.glob("/tmp/desi_pchmm_chr*.npy")))
    if not files:
        print("No result files found."); sys.exit(1)
    all_t, all_n, all_blk = [], [], []
    blk_off = 0
    for f in files:
        try:
            d = np.load(f, allow_pickle=True).item()
        except Exception as e:
            print(f"  Skip {f}: {e}"); continue
        if "win_t" not in d:
            print(f"  Skip {f}: old format (no win_t key)"); continue
        chrom = d.get("chrom", os.path.basename(f))
        wt = d["win_t"]; wn = d["win_n"]
        wb = d["win_blk"].astype(np.int64) + blk_off
        blk_off = int(wb.max()) + 1000
        all_t.append(wt); all_n.append(wn); all_blk.append(wb)
        print(f"  {chrom}: {wt.shape[0]} windows  blk {int(wb.min())}-{int(wb.max())}")
    win_t = np.vstack(all_t)
    win_n = np.vstack(all_n)
    win_blk = np.concatenate(all_blk)
    print(f"Total: {win_t.shape[0]} windows")
    return win_t, win_n, win_blk


def block_jk(vals, blks):
    ok = ~np.isnan(vals)
    v, b = vals[ok], blks[ok]
    if len(v) == 0:
        return np.nan, np.nan
    ublks = np.unique(b)
    n = len(ublks)
    if n < 2:
        return float(np.nanmean(v)), np.nan
    theta = float(np.nanmean(v))
    jk = np.array([np.nanmean(v[b != bk]) for bk in ublks])
    pseudo = n * theta - (n - 1) * jk
    se = float(np.sqrt(np.var(pseudo, ddof=1) / n))
    return theta, se


def compute_stats(win_t, win_blk):
    stats = {}
    for gi, gn in enumerate(GROUP_NAMES):
        m, se = block_jk(win_t[:, gi], win_blk)
        nok = int((~np.isnan(win_t[:, gi])).sum())
        stats[gn] = (m, se, nok)
    return stats


def expected_bottle(n=300000, seed=42):
    """Simulate per-window mean T under panmictic 930 Ka bottleneck."""
    rng = np.random.default_rng(seed)
    Ne_bot = 1280; T_start = 930000 / GEN_TIME; T_end = (930000 - 117000) / GEN_TIME
    Ne_anc = 98000; rate = 1.0 / (2 * Ne_bot); dur = T_start - T_end
    p = 1 - np.exp(-rate * dur)
    u = rng.random(n); ib = u < p
    T = np.empty(n)
    T[ib] = (rng.exponential(1 / rate, ib.sum()).clip(0, dur) + T_end) * GEN_TIME / 1000
    T[~ib] = (T_start + rng.exponential(1.0 / (2 * Ne_anc), (~ib).sum())) * GEN_TIME / 1000
    return T


def fig1(win_t, stats):
    fig, axes = plt.subplots(1, 5, figsize=(18, 4))
    bins = np.linspace(0, 3000, 60)
    for gi in WITHIN_GROUPS:
        ax = axes[gi]; gn = GROUP_NAMES[gi]; col = GROUP_COLOR[gn]
        v = win_t[:, gi]; v = v[~np.isnan(v)]
        ax.hist(v, bins=bins, color=col, alpha=0.75, density=True, edgecolor="none")
        m, se, nw = stats[gn]
        ax.axvline(m, color="k", lw=1.5, label=f"mean {m:.0f} Ka")
        ax.axvline(930, color="r", lw=1.2, ls="--", label="930 Ka")
        if not np.isnan(se):
            ax.axvspan(m - 2 * se, m + 2 * se, alpha=0.15, color="k")
        ax.set_title(gn.replace("within-", ""), fontsize=12, fontweight="bold")
        ax.set_xlabel("Mean TMRCA (Ka)", fontsize=9)
        if gi == 0:
            ax.set_ylabel("Density", fontsize=9)
        ax.legend(fontsize=7)
        ax.text(0.97, 0.97, f"n={nw}", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    plt.suptitle("Per-window mean TMRCA (100 kb windows, all pairs aggregated, shaded=2SE)", fontsize=11, y=1.02)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{OUT_DIR}/fig1_winhist_within.{ext}", dpi=180, bbox_inches="tight")
    plt.close()
    print("  fig1 saved")


def fig2(win_t, win_blk, stats):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    bins = np.linspace(0, 2500, 50)
    # Panel A: within-AFR obs vs expected
    ax = axes[0]
    v = win_t[:, 0]; v = v[~np.isnan(v)]
    ax.hist(v, bins=bins, density=True, color="#d62728", alpha=0.6,
            label=f"within-AFR obs (n={len(v)})")
    exp = expected_bottle()
    ax.hist(exp, bins=bins, density=True, color="gray", alpha=0.5,
            histtype="step", lw=2, ls="--", label="expected: panmictic 930 Ka")
    ax.axvline(930, color="r", lw=1.5, ls=":", label="930 Ka")
    m, se, _ = stats["within-AFR"]
    lbl = f"obs mean {m:.0f}+/-{se:.0f} Ka" if not np.isnan(se) else f"obs mean {m:.0f} Ka"
    ax.axvline(m, color="#d62728", lw=1.8, label=lbl)
    ax.set_title("within-AFR: observed vs expected", fontweight="bold")
    ax.set_xlabel("TMRCA (Ka)"); ax.set_ylabel("Density"); ax.legend(fontsize=7)
    # Panel B: bar chart all within groups
    ax = axes[1]
    means, ses, labels, colors = [], [], [], []
    for gi in WITHIN_GROUPS:
        gn = GROUP_NAMES[gi]; m2, s2, _ = stats[gn]
        if not np.isnan(m2):
            means.append(m2); ses.append(s2 if not np.isnan(s2) else 0)
            labels.append(gn.replace("within-", "")); colors.append(GROUP_COLOR[gn])
    y = np.arange(len(means))
    ax.barh(y, means, xerr=[2 * s for s in ses], color=colors, alpha=0.75,
            error_kw=dict(ecolor="k", capsize=4, lw=1.2))
    ax.axvline(930, color="r", lw=1.5, ls="--", label="930 Ka")
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel("Mean TMRCA (Ka)"); ax.set_title("Mean +/- 2 jackknife SE", fontweight="bold")
    ax.legend(fontsize=8)
    # Panel C: CDF AFR vs EAS
    ax = axes[2]
    for gi, col, lbl in zip([0, 1], ["#d62728", "#1f77b4"], ["within-AFR", "within-EAS"]):
        v = win_t[:, gi]; v = np.sort(v[~np.isnan(v)])
        ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=col, lw=1.8, label=lbl)
    ax.axvline(930, color="r", lw=1.2, ls="--", label="930 Ka")
    ax.set_xlabel("TMRCA (Ka)"); ax.set_ylabel("CDF")
    ax.set_title("AFR vs EAS CDF", fontweight="bold")
    ax.legend(fontsize=8); ax.set_xlim(0, 2500)
    plt.suptitle("930 Ka bottleneck test (100 kb windows, per-group aggregate TMRCA)", fontsize=12, y=1.02)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{OUT_DIR}/fig2_bottleneck_test.{ext}", dpi=180, bbox_inches="tight")
    plt.close()
    print("  fig2 saved")


def fig3(stats):
    fig, ax = plt.subplots(figsize=(9, 5))
    for gi in WITHIN_GROUPS:
        gn = GROUP_NAMES[gi]; m, se, nw = stats[gn]
        if np.isnan(m): continue
        ax.errorbar(m, gi, xerr=2 * se if not np.isnan(se) else 0,
                    fmt="o", color=GROUP_COLOR[gn], ms=8, capsize=5,
                    label=f"{gn.replace('within-', '')}: {m:.0f}+/-{se:.0f} Ka")
    ax.axvline(930, color="r", ls="--", lw=1.5, label="FitCoal 930 Ka")
    ax.axvline(65, color="steelblue", ls=":", lw=1.2, label="OOA ~65 Ka")
    ax.set_xlabel("Mean per-window TMRCA (Ka)", fontsize=11)
    ax.set_yticks(WITHIN_GROUPS)
    ax.set_yticklabels([GROUP_NAMES[i].replace("within-", "") for i in WITHIN_GROUPS])
    ax.set_title("Population mean TMRCA with block-jackknife 95% CI", fontsize=11)
    ax.legend(fontsize=8); ax.set_xlim(0, 2000)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{OUT_DIR}/fig3_net.{ext}", dpi=180, bbox_inches="tight")
    plt.close()
    print("  fig3 saved")


def fig4(stats):
    mat = np.full((5, 5), np.nan)
    gmap = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (3, 3): 3, (4, 4): 4,
            (0, 1): 5, (0, 2): 6, (0, 3): 7, (0, 4): 8,
            (1, 2): 9, (1, 3): 10, (1, 4): 11,
            (2, 3): 12, (2, 4): 13, (3, 4): 14}
    for (r, c), gi in gmap.items():
        m, _, _ = stats[GROUP_NAMES[gi]]
        mat[r, c] = mat[c, r] = m
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(mat, cmap="YlOrRd", vmin=np.nanmin(mat) * 0.85)
    plt.colorbar(im, ax=ax, label="Mean TMRCA (Ka)")
    labs = ["AFR", "EAS", "EUR", "AMR", "SAS"]
    ax.set_xticks(range(5)); ax.set_xticklabels(labs)
    ax.set_yticks(range(5)); ax.set_yticklabels(labs)
    ax.set_title("TMRCA matrix (100 kb windows, all pairs)", fontsize=10)
    for r in range(5):
        for c in range(5):
            if not np.isnan(mat[r, c]):
                ax.text(c, r, f"{mat[r, c]:.0f}", ha="center", va="center",
                        fontsize=8.5, fontweight="bold",
                        color="white" if mat[r, c] > np.nanmedian(mat) * 1.3 else "k")
    plt.tight_layout()
    for ext in ("png", "pdf"):
        plt.savefig(f"{OUT_DIR}/fig4_tmrca_matrix.{ext}", dpi=180, bbox_inches="tight")
    plt.close()
    print("  fig4 saved")


def print_summary(stats):
    print()
    print(f"{'Group':22s}  {'Mean(Ka)':>10}  {'SE(Ka)':>8}  {'Windows':>8}")
    for gn in GROUP_NAMES:
        m, se, n = stats[gn]
        if np.isnan(m): continue
        s_str = f"{se:.1f}" if not np.isnan(se) else "n/a"
        print(f"  {gn:20s}  {m:>10.1f}  {s_str:>8}  {n:>8}")
    m, se, _ = stats["within-AFR"]
    if not np.isnan(se) and se > 0:
        z = abs(m - 930) / se
        print(f"\nwithin-AFR: {m:.1f} +/- {se:.1f} Ka  =>  {z:.1f} SE from 930 Ka")
        print("  =>", "SIGNIFICANT departure from panmictic 930 Ka" if z > 3 else "consistent with 930 Ka")


def main():
    chroms = sys.argv[1:] if len(sys.argv) > 1 else None
    print("Loading data...")
    win_t, win_n, win_blk = load_data(chroms)
    print("Computing block-jackknife statistics...")
    stats = compute_stats(win_t, win_blk)
    print_summary(stats)
    print("\nGenerating figures...")
    fig1(win_t, stats)
    fig2(win_t, win_blk, stats)
    fig3(stats)
    fig4(stats)
    print(f"\nDone -> {OUT_DIR}/")


if __name__ == "__main__":
    main()
