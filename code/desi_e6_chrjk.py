"""
E6: Per-chromosome consistency check (BGS robustness proxy)

Rationale: If the AFR > EAS TMRCA signal is a genome-wide artifact (BGS,
centromere effects, mapping bias), it should disappear on some chromosomes.
True ancestry signal should be consistent across all 22 autosomes.

Also computes:
- Leave-one-chromosome-out jackknife (true BGS robustness test)
- Rank-correlation between chromosome length and TMRCA (sanity check)
- Acrocentric chromosomes (13,14,15,21,22) vs metacentric comparison
"""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

os.makedirs("/home/yanlin/popghistory/results/validation", exist_ok=True)

POP_MAP = {
    'YRI':'AFR','LWK':'AFR','GWD':'AFR','MSL':'AFR','ESN':'AFR',
    'ASW':'AMR','ACB':'AMR','MXL':'AMR','PUR':'AMR','CLM':'AMR','PEL':'AMR',
    'CHB':'EAS','JPT':'EAS','CHS':'EAS','CDX':'EAS','KHV':'EAS',
    'CEU':'EUR','TSI':'EUR','FIN':'EUR','GBR':'EUR','IBS':'EUR',
    'GIH':'SAS','PJL':'SAS','BEB':'SAS','STU':'SAS','ITU':'SAS'
}
CHR_LEN = {
    "chr1":248956422,"chr2":242193529,"chr3":198295559,"chr4":190214555,
    "chr5":181538259,"chr6":170805979,"chr7":159345973,"chr8":145138636,
    "chr9":138394717,"chr10":133797422,"chr11":135086622,"chr12":133275309,
    "chr13":114364328,"chr14":107043718,"chr15":101991189,"chr16":90338345,
    "chr17":83257441,"chr18":80373285,"chr19":58617616,"chr20":64444167,
    "chr21":46709983,"chr22":50818468
}
ACROCENTRIC = {"chr13","chr14","chr15","chr21","chr22"}

import glob
files = sorted(glob.glob("/tmp/desi_pchmm_chr*.npy"))
print(f"Found {len(files)} chromosome files")

# Load first file to get metadata
d0 = np.load(files[0], allow_pickle=True).item()
hap_pop = d0["hap_pop"] if "hap_pop" in d0 else None
pair_rows = d0.get("pair_rows", None)
pair_cols = d0.get("pair_cols", None)

if hap_pop is None:
    print("ERROR: hap_pop not found in npy file. Keys:", list(d0.keys()))
    raise SystemExit(1)

sup = np.array([POP_MAP.get(p,'UNK') for p in hap_pop])
rs  = sup[pair_rows]; cs = sup[pair_cols]
afr_mask = (rs=="AFR") & (cs=="AFR")
eas_mask = (rs=="EAS") & (cs=="EAS")
eur_mask = (rs=="EUR") & (cs=="EUR")

print(f"AFR pairs: {afr_mask.sum()}, EAS pairs: {eas_mask.sum()}, EUR pairs: {eur_mask.sum()}")

# Per-chromosome analysis
chr_res = []
for fn in files:
    chrom = fn.split("_chr")[1].replace(".npy","")
    d = np.load(fn, allow_pickle=True).item()
    pt = d["pair_tmrca"]
    afr_m = pt[afr_mask].mean() if afr_mask.sum()>0 else np.nan
    eas_m = pt[eas_mask].mean() if eas_mask.sum()>0 else np.nan
    eur_m = pt[eur_mask].mean() if eur_mask.sum()>0 else np.nan
    n_snps = d.get("n_snps", np.nan)
    chlen  = CHR_LEN.get(f"chr{chrom}", CHR_LEN.get(chrom, np.nan))
    chr_res.append({"chrom":chrom, "afr":afr_m, "eas":eas_m, "eur":eur_m,
                    "diff":afr_m-eas_m, "length":chlen, "n_snps":n_snps,
                    "acrocentric":f"chr{chrom}" in ACROCENTRIC})
    print(f"  chr{chrom:>2}: AFR={afr_m:.1f}, EAS={eas_m:.1f}, diff={afr_m-eas_m:.1f} Ka  "
          f"({'acrocentric' if f'chr{chrom}' in ACROCENTRIC else 'metacentric'})")

# Genome-wide averages
afr_all = np.mean([r["afr"] for r in chr_res if not np.isnan(r["afr"])])
eas_all = np.mean([r["eas"] for r in chr_res if not np.isnan(r["eas"])])
diff_all = afr_all - eas_all
diffs = np.array([r["diff"] for r in chr_res if not np.isnan(r["diff"])])
print(f"\nGenome-wide mean: AFR={afr_all:.1f}, EAS={eas_all:.1f}, diff={diff_all:.1f} Ka")
print(f"Per-chromosome diff: mean={diffs.mean():.1f}, std={diffs.std():.1f}, min={diffs.min():.1f}, max={diffs.max():.1f}")
print(f"All chromosomes show AFR>EAS: {(diffs>0).all()}")

# Leave-one-out jackknife
chr_tmrca = {}
for fn in files:
    chrom = fn.split("_chr")[1].replace(".npy","")
    d = np.load(fn, allow_pickle=True).item()
    chr_tmrca[chrom] = d["pair_tmrca"]

all_tmrca = np.array([chr_tmrca[fn.split("_chr")[1].replace(".npy","")] for fn in files])
# shape: (22, n_pairs)
afr_jk = []; eas_jk = []; diff_jk = []
for i, fn in enumerate(files):
    chrom = fn.split("_chr")[1].replace(".npy","")
    loo = np.delete(all_tmrca, i, axis=0).mean(axis=0)
    afr_jk.append(loo[afr_mask].mean())
    eas_jk.append(loo[eas_mask].mean())
    diff_jk.append(loo[afr_mask].mean() - loo[eas_mask].mean())

diff_jk = np.array(diff_jk)
pseudo = 22 * diff_all - 21 * diff_jk
se_jk  = np.sqrt(np.var(pseudo, ddof=1) / 22)
print(f"\nLeave-one-out jackknife: diff={diff_all:.1f} ± {se_jk:.1f} Ka (SE)")
print(f"All LOO diffs: min={diff_jk.min():.1f}, max={diff_jk.max():.1f} Ka")

# Spearman correlation: chromosome length vs TMRCA
lens = [r["length"] for r in chr_res if not np.isnan(r["afr"])]
afrs = [r["afr"] for r in chr_res if not np.isnan(r["afr"])]
eass = [r["eas"] for r in chr_res if not np.isnan(r["eas"])]
rho_afr, p_afr = spearmanr(lens, afrs)
rho_eas, p_eas = spearmanr(lens, eass)
rho_diff, p_diff = spearmanr(lens, diffs)
print(f"\nSpearman r (chr length vs TMRCA): AFR r={rho_afr:.3f} p={p_afr:.3f}, "
      f"EAS r={rho_eas:.3f} p={p_eas:.3f}, diff r={rho_diff:.3f} p={p_diff:.3f}")

# Acrocentric vs metacentric
acro_diffs = [r["diff"] for r in chr_res if r["acrocentric"] and not np.isnan(r["diff"])]
meta_diffs = [r["diff"] for r in chr_res if not r["acrocentric"] and not np.isnan(r["diff"])]
print(f"\nAcrocentric chr: mean diff={np.mean(acro_diffs):.1f} Ka")
print(f"Metacentric chr: mean diff={np.mean(meta_diffs):.1f} Ka")

# Save
np.save("/home/yanlin/popghistory/results/validation/e6_results.npy",
        {"chr_res":chr_res, "diff_jk":diff_jk, "se_jk":se_jk,
         "diff_all":diff_all, "rho_afr":rho_afr, "rho_eas":rho_eas, "rho_diff":rho_diff})

# Figure
chroms_sorted = sorted(chr_res, key=lambda r: int(r["chrom"]) if r["chrom"].isdigit() else 99)
labels = [r["chrom"] for r in chroms_sorted]
afr_v  = [r["afr"]  for r in chroms_sorted]
eas_v  = [r["eas"]  for r in chroms_sorted]
diff_v = [r["diff"] for r in chroms_sorted]
colors = ["#FF6B6B" if r["acrocentric"] else "#4ECDC4" for r in chroms_sorted]

fig, axes = plt.subplots(2, 1, figsize=(14, 8))
fig.suptitle("E6: Per-chromosome TMRCA consistency (BGS robustness)", fontsize=12)

ax = axes[0]
x = np.arange(len(labels))
ax.plot(x, afr_v, "rs-", ms=7, label=f"within-AFR (mean={afr_all:.1f} Ka)")
ax.plot(x, eas_v, "bs-", ms=7, label=f"within-EAS (mean={eas_all:.1f} Ka)")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, fontsize=8)
ax.set_ylabel("Mean pairwise TMRCA (Ka)")
ax.set_title("(a) Per-chromosome within-group TMRCA")
ax.legend(fontsize=9)
for xi, r in enumerate(chroms_sorted):
    if r["acrocentric"]:
        ax.axvspan(xi-0.4, xi+0.4, alpha=0.12, color="orange", label="_nolegend_")
ax.axhline(afr_all, color="red", ls="--", lw=1, alpha=0.5)
ax.axhline(eas_all, color="blue", ls="--", lw=1, alpha=0.5)

ax = axes[1]
bars = ax.bar(x, diff_v, color=colors, alpha=0.85, edgecolor="white", lw=0.5)
ax.axhline(diff_all, color="black", ls="--", lw=2, label=f"Genome-wide mean ({diff_all:.1f} Ka)")
ax.axhline(0, color="grey", lw=0.5)
ax.fill_between(x[[0,-1]], diff_all-se_jk, diff_all+se_jk, alpha=0.15, color="black",
                label=f"Jackknife SE=±{se_jk:.1f} Ka")
ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, fontsize=8)
ax.set_ylabel("AFR − EAS TMRCA (Ka)")
ax.set_title(f"(b) AFR-EAS difference per chromosome  (teal=metacentric, coral=acrocentric)")
ax.legend(fontsize=9)

plt.tight_layout()
plt.savefig("/home/yanlin/popghistory/results/validation/fig_e6.png", dpi=150, bbox_inches="tight")
print("\nSaved: results/validation/fig_e6.png")
print("Saved: results/validation/e6_results.npy")
