"""
plot_bottleneck_scan.py
───────────────────────
绘制瓶颈参数扫描结果：
  Fig A: P(AFR > 930 ka) 热图（anc=20k），标注 FitCoal 点和观测线
  Fig B: P(AFR > 930 ka) 热图（anc=50k），同上
  Fig C: AFR 均值 TMRCA 热图（两个 anc_ne）
  Bottom panel: FitCoal 参数下预测 vs 观测的直接对比条形图
"""

import pickle, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import os

OUTDIR   = "/home/yanlin/popghistory/results/bottleneck_scan"
SCAN_PKL = os.path.join(OUTDIR, "scan_results.pkl")

with open(SCAN_PKL, "rb") as f:
    data = pickle.load(f)

df  = pd.DataFrame(data["results"])
OBS_P   = data["obs_p_afr_gt930"]
OBS_AFR = data["obs_afr_mean"]
OBS_DIFF = data["obs_diff"]
NE_GRID  = data["ne_bn_grid"]
DUR_GRID = data["dur_ka_grid"]

print(f"Loaded {len(df)} parameter combinations")
print(f"Observed: P(>930)={OBS_P:.3f}, AFR mean={OBS_AFR:.0f} ka, AFR-EAS={OBS_DIFF:.0f} ka")

# FitCoal exact point
FC_NE = 1280; FC_DUR = 117
fitcoal_rows = df[(df.ne_bn == FC_NE) & (df.dur_ka == FC_DUR)]

# ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle("Bottleneck Parameter Scan: P(AFR TMRCA > 930 ka)\nFitCoal Parameters vs. Observed Multi-population TMRCA Distribution",
             fontsize=14, fontweight='bold', y=0.98)

ne_labels  = [f"{n//1000}k" for n in NE_GRID]
dur_labels = [str(d) for d in DUR_GRID]

def make_heatmap(ax, anc_ne, metric, title, obs_val, cmap, vmin, vmax, fmt=".2f"):
    sub = df[df.anc_ne == anc_ne]
    mat = np.zeros((len(NE_GRID), len(DUR_GRID)))
    for i, ne in enumerate(NE_GRID):
        for j, dur in enumerate(DUR_GRID):
            row = sub[(sub.ne_bn == ne) & (sub.dur_ka == dur)]
            mat[i, j] = row[metric].values[0] if len(row) else np.nan
    im = ax.imshow(mat, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, origin='upper')
    ax.set_xticks(range(len(DUR_GRID)));  ax.set_xticklabels(dur_labels, fontsize=9)
    ax.set_yticks(range(len(NE_GRID)));   ax.set_yticklabels(ne_labels, fontsize=9)
    ax.set_xlabel("Bottleneck duration (ka)", fontsize=10)
    ax.set_ylabel("Bottleneck Ne", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    # Mark FitCoal point
    fc_j = DUR_GRID.index(FC_DUR)
    fc_i = NE_GRID.index(FC_NE)
    ax.add_patch(Rectangle((fc_j-0.5, fc_i-0.5), 1, 1,
                            fill=False, edgecolor='red', linewidth=3))
    ax.text(fc_j, fc_i, f"FitCoal\n{fmt.format(mat[fc_i, fc_j]) if not np.isnan(mat[fc_i,fc_j]) else 'N/A'}",
            ha='center', va='center', fontsize=8, color='red', fontweight='bold')
    # Add value annotations
    for i in range(len(NE_GRID)):
        for j in range(len(DUR_GRID)):
            v = mat[i, j]
            if not np.isnan(v):
                color = 'white' if v < (vmin + (vmax-vmin)*0.5) else 'black'
                ax.text(j, i, format(v, fmt), ha='center', va='center', fontsize=7.5, color=color)
    plt.colorbar(im, ax=ax, shrink=0.7)
    return mat, im

# Row 1: P(AFR > 930 ka)
for col, anc_ne in enumerate([20000, 50000]):
    mat, _ = make_heatmap(
        axes[0, col], anc_ne, "p_afr_gt930",
        f"P(AFR win TMRCA > 930 ka)\nanc Ne = {anc_ne//1000}k",
        OBS_P, "RdYlGn", 0, 1, ".2f"
    )

# Row 1, col 2: direct comparison bar chart
ax_bar = axes[0, 2]
categories = ["FitCoal\n(Ne=1280, 117ka)\nanc=20k",
              "FitCoal\n(Ne=1280, 117ka)\nanc=50k",
              "Observed\n(1kGP DESI)"]
fitcoal_20k = df[(df.anc_ne==20000)&(df.ne_bn==FC_NE)&(df.dur_ka==FC_DUR)]["p_afr_gt930"].values
fitcoal_50k = df[(df.anc_ne==50000)&(df.ne_bn==FC_NE)&(df.dur_ka==FC_DUR)]["p_afr_gt930"].values

vals = [fitcoal_20k[0] if len(fitcoal_20k) else 0,
        fitcoal_50k[0] if len(fitcoal_50k) else 0,
        OBS_P]
colors_bar = ['#d62728', '#d62728', '#2ca02c']
bars = ax_bar.bar(categories, vals, color=colors_bar, edgecolor='black', linewidth=1.5, width=0.5)
ax_bar.axhline(OBS_P, color='green', linestyle='--', linewidth=2, label=f'Observed ({OBS_P:.2f})')
ax_bar.set_ylim(0, 1.05)
ax_bar.set_ylabel("P(AFR win TMRCA > 930 ka)", fontsize=10)
ax_bar.set_title("FitCoal Prediction vs. Observed", fontsize=11, fontweight='bold')
ax_bar.legend(fontsize=9)
for bar, val in zip(bars, vals):
    ax_bar.text(bar.get_x() + bar.get_width()/2, val + 0.02,
                f'{val:.2f}', ha='center', fontsize=11, fontweight='bold')
ax_bar.tick_params(axis='x', labelsize=9)

# Row 2: AFR mean TMRCA
for col, anc_ne in enumerate([20000, 50000]):
    make_heatmap(
        axes[1, col], anc_ne, "afr_mean",
        f"AFR Mean TMRCA (ka)\nanc Ne = {anc_ne//1000}k",
        OBS_AFR, "YlOrRd", 500, 2500, ".0f"
    )
    axes[1, col].set_title(axes[1, col].get_title() + f"\n(observed = {OBS_AFR:.0f} ka)", fontsize=10)

# Row 2, col 2: scatter plot of P vs AFR mean for all params
ax_sc = axes[1, 2]
for anc_ne, marker, color in [(20000, 'o', '#1f77b4'), (50000, 's', '#ff7f0e')]:
    sub = df[df.anc_ne == anc_ne]
    sc = ax_sc.scatter(sub["afr_mean"], sub["p_afr_gt930"],
                       c=np.log10(sub["ne_bn"]), cmap='viridis',
                       marker=marker, s=60, alpha=0.7,
                       label=f'anc={anc_ne//1000}k')
# FitCoal points
for anc_ne in [20000, 50000]:
    row = df[(df.anc_ne==anc_ne)&(df.ne_bn==FC_NE)&(df.dur_ka==FC_DUR)]
    if len(row):
        ax_sc.scatter(row["afr_mean"], row["p_afr_gt930"],
                      color='red', marker='*', s=300, zorder=5,
                      label=f'FitCoal ({anc_ne//1000}k)' if anc_ne==20000 else '_')
# Observed point
ax_sc.scatter([OBS_AFR], [OBS_P], color='green', marker='D', s=200, zorder=6, label='Observed')
ax_sc.axhline(OBS_P, color='green', linestyle='--', alpha=0.6)
ax_sc.axvline(OBS_AFR, color='green', linestyle='--', alpha=0.6)
ax_sc.set_xlabel("AFR Mean TMRCA (ka)", fontsize=10)
ax_sc.set_ylabel("P(AFR TMRCA > 930 ka)", fontsize=10)
ax_sc.set_title("Parameter Space vs. Observed", fontsize=11, fontweight='bold')
ax_sc.legend(fontsize=8, loc='upper left')

plt.tight_layout(rect=[0, 0, 1, 0.97])
out_fig = os.path.join(OUTDIR, "bottleneck_scan_heatmap.png")
plt.savefig(out_fig, dpi=150, bbox_inches='tight')
print(f"\nFigure saved: {out_fig}")

# ─── Summary statistics ───────────────────────────────────────────────
print("\n" + "=" * 75)
print("Key result: FitCoal (Ne=1280, 117 ka) vs. Observed")
print("=" * 75)
for anc_ne in [20000, 50000]:
    row = df[(df.anc_ne==anc_ne)&(df.ne_bn==FC_NE)&(df.dur_ka==FC_DUR)]
    if len(row):
        print(f"  anc={anc_ne//1000}k: AFR mean={row.afr_mean.values[0]:.0f} ka, "
              f"P(>930)={row.p_afr_gt930.values[0]:.3f}")
print(f"  Observed:  AFR mean={OBS_AFR:.0f} ka, P(>930)={OBS_P:.3f}")

# Compatible region (both P and mean within tolerance)
print("\nCompatible region (|ΔP|<0.10 AND |ΔAFR|<200 ka):")
compat = df[(df.delta_p_gt930 < 0.10) & (df.delta_afr_mean < 200)]
if len(compat):
    print(compat[["anc_ne","ne_bn","dur_ka","afr_mean","p_afr_gt930"]].to_string(index=False))
else:
    print("  No combination passes both criteria simultaneously.")
    # Relax to find closest
    df["combined_score"] = df["delta_p_gt930"] / OBS_P + df["delta_afr_mean"] / OBS_AFR
    top5 = df.nsmallest(5, "combined_score")
    print("  Top 5 closest combinations:")
    print(top5[["anc_ne","ne_bn","dur_ka","afr_mean","p_afr_gt930","delta_p_gt930","delta_afr_mean"]].to_string(index=False))
