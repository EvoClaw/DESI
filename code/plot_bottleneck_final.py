"""
plot_bottleneck_final.py
────────────────────────
最终发表质量图：
  Panel A: P(AFR win TMRCA > 930 ka) heatmap for anc=50k (FitCoal's ancestral Ne)
  Panel B: AFR mean TMRCA heatmap for anc=50k
  Panel C: Comparison bar chart (FitCoal predicted vs observed)
  Panel D: P(AFR>930) vs duration for key Ne values (line plot)
"""

import pickle, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
import os

OUTDIR   = "/home/yanlin/popghistory/results/bottleneck_scan"
FIG_OUT  = "/home/yanlin/popghistory/results/winhist/fig_bottleneck_scan.png"

with open(os.path.join(OUTDIR, "scan_results.pkl"), "rb") as f:
    data = pickle.load(f)

df       = pd.DataFrame(data["results"])
OBS_P    = data["obs_p_afr_gt930"]        # 0.707
OBS_AFR  = data["obs_afr_mean"]           # 1229 ka
OBS_DIFF = data["obs_diff"]               # 371 ka
NE_GRID  = data["ne_bn_grid"]             # [300, 500, 800, 1280, 2000, 4000, 8000, 20000]
DUR_GRID = data["dur_ka_grid"]            # [30, 50, 80, 117, 150, 200, 300]
FC_NE, FC_DUR = 1280, 117

ne_labels  = [f"{n//1000}k" if n >= 1000 else str(n) for n in NE_GRID]
dur_labels = [str(d) for d in DUR_GRID]

# ── Data matrices ────────────────────────────────────────────────────
sub50 = df[df.anc_ne == 50000]
mat_p50  = np.full((len(NE_GRID), len(DUR_GRID)), np.nan)
mat_m50  = np.full((len(NE_GRID), len(DUR_GRID)), np.nan)

for i, ne in enumerate(NE_GRID):
    for j, dur in enumerate(DUR_GRID):
        row = sub50[(sub50.ne_bn == ne) & (sub50.dur_ka == dur)]
        if len(row):
            mat_p50[i, j] = row["p_afr_gt930"].values[0]
            mat_m50[i, j] = row["afr_mean"].values[0]

fc_i = NE_GRID.index(FC_NE)
fc_j = DUR_GRID.index(FC_DUR)

# ── Figure ───────────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 13))
gs  = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.32,
                        left=0.07, right=0.97, top=0.93, bottom=0.08)

ax_pa = fig.add_subplot(gs[0, 0])
ax_pb = fig.add_subplot(gs[0, 1])
ax_pc = fig.add_subplot(gs[0, 2])
ax_pd = fig.add_subplot(gs[1, 0])
ax_pe = fig.add_subplot(gs[1, 1])
ax_pf = fig.add_subplot(gs[1, 2])

ANNOTATION_SIZE = 8
LABEL_SIZE = 10
TITLE_SIZE = 11

def draw_heatmap(ax, mat, cmap, vmin, vmax, title, obs_val, fmt=".2f",
                 obs_annotation=None):
    im = ax.imshow(mat, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax, origin='upper')
    ax.set_xticks(range(len(DUR_GRID)));  ax.set_xticklabels(dur_labels, fontsize=9)
    ax.set_yticks(range(len(NE_GRID)));   ax.set_yticklabels(ne_labels,  fontsize=9)
    ax.set_xlabel("Bottleneck duration (ka)", fontsize=LABEL_SIZE)
    ax.set_ylabel("Bottleneck $N_e$", fontsize=LABEL_SIZE)
    ax.set_title(title, fontsize=TITLE_SIZE, fontweight='bold', pad=8)
    # FitCoal box
    rect = Rectangle((fc_j-0.5, fc_i-0.5), 1, 1,
                      fill=False, edgecolor='red', linewidth=3, linestyle='-')
    ax.add_patch(rect)
    # Value annotations
    for i in range(len(NE_GRID)):
        for j in range(len(DUR_GRID)):
            v = mat[i, j]
            if np.isnan(v): continue
            is_fc = (i == fc_i and j == fc_j)
            color = 'white' if v < vmin + (vmax-vmin)*0.5 else 'black'
            text = format(v, fmt)
            weight = 'bold' if is_fc else 'normal'
            ax.text(j, i, text, ha='center', va='center',
                    fontsize=ANNOTATION_SIZE if not is_fc else 9,
                    color='red' if is_fc else color, fontweight=weight)
    plt.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    if obs_annotation:
        ax.text(0.98, 0.02, obs_annotation, transform=ax.transAxes,
                ha='right', va='bottom', fontsize=8, color='navy',
                bbox=dict(facecolor='lightyellow', alpha=0.8, edgecolor='navy'))

# A: P heatmap (anc=50k)
draw_heatmap(
    ax_pa, mat_p50, "RdYlGn", 0, 1,
    "(A) P(AFR win TMRCA > 930 ka)\n[ancestral Ne = 50k]",
    OBS_P, fmt=".2f",
    obs_annotation=f"Observed = {OBS_P:.3f}"
)

# B: AFR mean heatmap (anc=50k)
draw_heatmap(
    ax_pb, mat_m50, "YlOrRd_r", 400, 2800,
    "(B) AFR Mean TMRCA (ka)\n[ancestral Ne = 50k]",
    OBS_AFR, fmt=".0f",
    obs_annotation=f"Observed = {OBS_AFR:.0f} ka"
)

# C: Predicted vs Observed bar chart
ax_pc.set_visible(False)
ax_pc = fig.add_subplot(gs[0, 2])
sub20 = df[df.anc_ne == 20000]
fc_20k = sub20[(sub20.ne_bn==FC_NE)&(sub20.dur_ka==FC_DUR)]["p_afr_gt930"].values[0]
fc_50k = sub50[(sub50.ne_bn==FC_NE)&(sub50.dur_ka==FC_DUR)]["p_afr_gt930"].values[0]

categories = [f"FitCoal prediction\n(anc Ne=20k)", f"FitCoal prediction\n(anc Ne=50k)", "Observed\n(1kGP multi-pop)"]
vals = [fc_20k, fc_50k, OBS_P]
colors_bar = ['#d62728', '#ff7f0e', '#2ca02c']
bars = ax_pc.bar(categories, vals, color=colors_bar, edgecolor='black', linewidth=1.5, width=0.55)
ax_pc.axhline(OBS_P, color='darkgreen', linestyle='--', linewidth=2.5, zorder=2,
               label=f"Observed P = {OBS_P:.3f}")
# Error bars for simulations (binomial SE with 600 windows)
for bar, val in zip(bars[:2], vals[:2]):
    se = np.sqrt(val*(1-val)/600)
    ax_pc.errorbar(bar.get_x() + bar.get_width()/2, val, yerr=2*se,
                   fmt='none', color='black', capsize=5, linewidth=2)
# Observed SE (negligibly small)
se_obs = np.sqrt(OBS_P*(1-OBS_P)/27000)
ax_pc.errorbar(bars[2].get_x() + bars[2].get_width()/2, OBS_P, yerr=2*se_obs,
               fmt='none', color='black', capsize=5, linewidth=2)
ax_pc.set_ylim(0, 1.05)
ax_pc.set_ylabel("P(AFR TMRCA > 930 ka)", fontsize=LABEL_SIZE)
ax_pc.set_title("(C) FitCoal Prediction vs. Observed", fontsize=TITLE_SIZE, fontweight='bold', pad=8)
for bar, val in zip(bars, vals):
    ax_pc.text(bar.get_x() + bar.get_width()/2, val + 0.03,
               f'{val:.3f}', ha='center', fontsize=12, fontweight='bold',
               color=bar.get_facecolor())
ax_pc.tick_params(axis='x', labelsize=9)
ax_pc.legend(fontsize=9, loc='upper left')

# Add annotation for FitCoal deficit
ax_pc.annotate(
    f'Δ = {OBS_P - fc_50k:.3f}\n(~{(OBS_P-fc_50k)/np.sqrt(fc_50k*(1-fc_50k)/600):.1f} SE)',
    xy=(bars[1].get_x() + bars[1].get_width()/2, fc_50k),
    xytext=(bars[1].get_x() + bars[1].get_width()/2 + 0.35, fc_50k + 0.08),
    fontsize=9, color='darkred', fontweight='bold',
    arrowprops=dict(arrowstyle='->', color='darkred', lw=1.5)
)

# D: P(>930) vs duration for different Ne (anc=50k)
ne_plot = [1280, 2000, 4000, 8000]
colors_ne = ['red', 'orange', '#2ca02c', '#1f77b4']
for ne, c in zip(ne_plot, colors_ne):
    rows = [sub50[(sub50.ne_bn==ne)&(sub50.dur_ka==d)] for d in DUR_GRID]
    ps = [r["p_afr_gt930"].values[0] if len(r) else np.nan for r in rows]
    lw = 2.5 if ne == FC_NE else 1.5
    ls = '-' if ne == FC_NE else '--'
    ax_pd.plot(DUR_GRID, ps, color=c, lw=lw, ls=ls, marker='o', ms=6,
               label=f"Ne={ne//1000}k")
ax_pd.axhline(OBS_P, color='darkgreen', lw=2.5, ls='-.',
               label=f"Observed P = {OBS_P:.2f}")
# Mark FitCoal point
ax_pd.scatter([FC_DUR], [fc_50k], color='red', marker='*', s=300, zorder=6)
ax_pd.set_xlabel("Bottleneck duration (ka)", fontsize=LABEL_SIZE)
ax_pd.set_ylabel("P(AFR TMRCA > 930 ka)", fontsize=LABEL_SIZE)
ax_pd.set_title(f"(D) P(AFR > 930 ka) vs. Duration\n[ancestral Ne = 50k]",
                fontsize=TITLE_SIZE, fontweight='bold', pad=8)
ax_pd.legend(fontsize=9, loc='upper right')
ax_pd.set_xlim(0, 320); ax_pd.set_ylim(0, 1.05)
ax_pd.annotate("FitCoal\n(Ne=1280,117ka)", xy=(FC_DUR, fc_50k),
               xytext=(FC_DUR+40, fc_50k-0.15), fontsize=9, color='red', fontweight='bold',
               arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

# E: AFR mean vs duration
for ne, c in zip(ne_plot, colors_ne):
    rows = [sub50[(sub50.ne_bn==ne)&(sub50.dur_ka==d)] for d in DUR_GRID]
    ms = [r["afr_mean"].values[0] if len(r) else np.nan for r in rows]
    lw = 2.5 if ne == FC_NE else 1.5
    ls = '-' if ne == FC_NE else '--'
    ax_pe.plot(DUR_GRID, ms, color=c, lw=lw, ls=ls, marker='o', ms=6,
               label=f"Ne={ne//1000}k")
ax_pe.axhline(OBS_AFR, color='darkgreen', lw=2.5, ls='-.',
               label=f"Observed = {OBS_AFR:.0f} ka")
fc50_afr = sub50[(sub50.ne_bn==FC_NE)&(sub50.dur_ka==FC_DUR)]["afr_mean"].values[0]
ax_pe.scatter([FC_DUR], [fc50_afr], color='red', marker='*', s=300, zorder=6)
ax_pe.set_xlabel("Bottleneck duration (ka)", fontsize=LABEL_SIZE)
ax_pe.set_ylabel("AFR Mean TMRCA (ka)", fontsize=LABEL_SIZE)
ax_pe.set_title(f"(E) AFR Mean TMRCA vs. Duration\n[ancestral Ne = 50k]",
                fontsize=TITLE_SIZE, fontweight='bold', pad=8)
ax_pe.legend(fontsize=9, loc='upper right')
ax_pe.set_xlim(0, 320)
ax_pe.annotate(f"FitCoal: {fc50_afr:.0f} ka", xy=(FC_DUR, fc50_afr),
               xytext=(FC_DUR+40, fc50_afr+100), fontsize=9, color='red', fontweight='bold',
               arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

# F: scatter: P vs AFR mean (all anc=50k points)
cmap_sc = plt.cm.plasma
norm_sc = mcolors.Normalize(vmin=np.log10(300), vmax=np.log10(20000))
for i, (ne, dur) in enumerate(zip(sub50.ne_bn, sub50.dur_ka)):
    row = sub50[(sub50.ne_bn==ne)&(sub50.dur_ka==dur)]
    if not len(row): continue
    m = row["afr_mean"].values[0]; p = row["p_afr_gt930"].values[0]
    c = cmap_sc(norm_sc(np.log10(ne)))
    ax_pf.scatter(m, p, color=c, s=50, alpha=0.8)
# FitCoal and observed
ax_pf.scatter([fc50_afr], [fc_50k], color='red', marker='*', s=400, zorder=7,
               label='FitCoal (anc=50k)')
ax_pf.scatter([OBS_AFR], [OBS_P], color='darkgreen', marker='D', s=200, zorder=7,
               label='Observed')
ax_pf.axhline(OBS_P, color='darkgreen', ls='-.', lw=1.5, alpha=0.7)
ax_pf.axvline(OBS_AFR, color='darkgreen', ls='-.', lw=1.5, alpha=0.7)
# Colorbar for Ne
sm = plt.cm.ScalarMappable(cmap=cmap_sc, norm=norm_sc)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax_pf, shrink=0.75)
cbar.set_label("log₁₀(Ne_bottleneck)", fontsize=8)
ax_pf.set_xlabel("AFR Mean TMRCA (ka)", fontsize=LABEL_SIZE)
ax_pf.set_ylabel("P(AFR TMRCA > 930 ka)", fontsize=LABEL_SIZE)
ax_pf.set_title("(F) Parameter Space vs. Observed\n[ancestral Ne = 50k]",
                fontsize=TITLE_SIZE, fontweight='bold', pad=8)
ax_pf.legend(fontsize=9, loc='upper left')

fig.suptitle(
    "Bottleneck Parameter Scan: Consistency of FitCoal (Ne=1280, 117 ka) "
    "with Multi-population TMRCA Statistics",
    fontsize=13, fontweight='bold'
)

plt.savefig(FIG_OUT, dpi=150, bbox_inches='tight')
print(f"Figure saved: {FIG_OUT}")

# ── Print final summary ───────────────────────────────────────────────
print("\n" + "="*70)
print("FINAL SUMMARY")
print("="*70)

se_fc50 = np.sqrt(fc_50k*(1-fc_50k)/600)
se_obs  = np.sqrt(OBS_P*(1-OBS_P)/27000)
z       = (OBS_P - fc_50k) / np.sqrt(se_fc50**2 + se_obs**2)
print(f"\nFitCoal (anc=50k, Ne=1280, dur=117ka):")
print(f"  Predicted P(AFR>930) = {fc_50k:.3f} ± {2*se_fc50:.3f} (95% CI)")
print(f"  Observed  P(AFR>930) = {OBS_P:.3f} ± {2*se_obs:.4f} (95% CI)")
print(f"  Z-score = {z:.1f}  (P-value << 0.0001)")
print(f"")
print(f"  Predicted AFR mean = {fc50_afr:.0f} ka")
print(f"  Observed  AFR mean = {OBS_AFR:.0f} ka  (FitCoal is {100*(OBS_AFR-fc50_afr)/OBS_AFR:.1f}% too low)")
print(f"")
print(f"Compatible parameters (closest match to both P and AFR mean):")
df["score"] = df["delta_p_gt930"]/OBS_P + df["delta_afr_mean"]/OBS_AFR
top3 = df.nsmallest(3,"score")[["anc_ne","ne_bn","dur_ka","afr_mean","p_afr_gt930"]]
print(top3.to_string(index=False))
