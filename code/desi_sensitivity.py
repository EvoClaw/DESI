"""E7+E8: Mutation rate sensitivity + AMR subpopulation stratification"""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("/home/yanlin/popghistory/results/validation", exist_ok=True)

d = np.load("/home/yanlin/popghistory/results/final/desi_genome_tmrca.npy", allow_pickle=True).item()
pair_tmrca = d["pair_tmrca"]
hap_pop    = d["hap_pop"]
pair_rows  = d["pair_rows"]
pair_cols  = d["pair_cols"]

# Fine-grained population map
pop_to_super = {
    'YRI':'AFR','LWK':'AFR','GWD':'AFR','MSL':'AFR','ESN':'AFR',
    'ASW':'AMR','ACB':'AMR','MXL':'AMR','PUR':'AMR','CLM':'AMR','PEL':'AMR',
    'CHB':'EAS','JPT':'EAS','CHS':'EAS','CDX':'EAS','KHV':'EAS',
    'CEU':'EUR','TSI':'EUR','FIN':'EUR','GBR':'EUR','IBS':'EUR',
    'GIH':'SAS','PJL':'SAS','BEB':'SAS','STU':'SAS','ITU':'SAS'
}
sup  = np.array([pop_to_super.get(p,'UNK') for p in hap_pop])
fine = np.array(hap_pop)  # keep fine-grained pop label
rs   = sup[pair_rows];  cs  = sup[pair_cols]
rf   = fine[pair_rows]; cf  = fine[pair_cols]

MU_BASE = 1.2e-8

# ==== E7: Mutation rate sensitivity ====
print("=== E7: Mutation rate sensitivity ===")
print(f"{'mu':<15} {'AFR (Ka)':>10} {'EAS (Ka)':>10} {'Diff (Ka)':>12} {'Ratio':>8}")
results_mu = []
for mu in [1.0e-8, 1.2e-8, 1.4e-8, 1.6e-8]:
    scale = MU_BASE / mu  # T ∝ 1/mu
    t_scaled = pair_tmrca * scale
    afr_m = t_scaled[(rs=="AFR")&(cs=="AFR")].mean()
    eas_m = t_scaled[(rs=="EAS")&(cs=="EAS")].mean()
    eur_m = t_scaled[(rs=="EUR")&(cs=="EUR")].mean()
    diff  = afr_m - eas_m
    results_mu.append({"mu":mu,"afr":afr_m,"eas":eas_m,"eur":eur_m,"diff":diff})
    print(f"  mu={mu:.1e}  AFR={afr_m:>9.1f}  EAS={eas_m:>9.1f}  diff={diff:>11.1f}  ratio={afr_m/eas_m:>7.3f}")

# ==== E8: AMR subpopulation analysis ====
print("\n=== E8: AMR subpopulation stratification ===")
AMR_SUBS = ['ASW','ACB','MXL','PUR','CLM','PEL']
sub_res = {}
for sub in AMR_SUBS:
    mask = (rf==sub)&(cf==sub)
    vals = pair_tmrca[mask]
    if len(vals) >= 5:
        sub_res[sub] = {"n":len(vals),"mean":vals.mean(),"std":vals.std(),"vals":vals}
        print(f"  {sub}: N={len(vals):>4} pairs, mean={vals.mean():.1f} Ka, std={vals.std():.1f} Ka")
    else:
        print(f"  {sub}: too few pairs ({len(vals)})")

# AMR overall
amr_all = pair_tmrca[(rs=="AMR")&(cs=="AMR")]
print(f"  AMR overall: N={len(amr_all)} pairs, mean={amr_all.mean():.1f} Ka, std={amr_all.std():.1f} Ka")
print("  Note: high AMR std likely driven by admixed populations (e.g. ASW/ACB with AFR ancestry)")

# Cross-group for AMR subpops
print("\nCross-group TMRCA for AMR subs:")
for sub in AMR_SUBS:
    for sp in ["AFR","EAS","EUR"]:
        mask = ((rf==sub)&(cs==sp))|((cf==sub)&(rs==sp))
        vals = pair_tmrca[mask]
        if len(vals)>=5:
            print(f"  {sub}-{sp}: N={len(vals)}, mean={vals.mean():.0f} Ka")

# Save
np.save("/home/yanlin/popghistory/results/validation/e7e8_results.npy",
        {"mu_sensitivity":results_mu, "amr_subs":
         {k:{kk:vv for kk,vv in v.items() if kk!="vals"} for k,v in sub_res.items()}})

# ---- Figures ----
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("E7+E8: Sensitivity + AMR Stratification", fontsize=12)

ax = axes[0]
mus  = [r["mu"] for r in results_mu]
afrs = [r["afr"] for r in results_mu]; eass = [r["eas"] for r in results_mu]
diffs= [r["diff"] for r in results_mu]
ax.plot(mus, afrs, "rs-", ms=7, label="within-AFR")
ax.plot(mus, eass, "bs-", ms=7, label="within-EAS")
ax2 = ax.twinx()
ax2.plot(mus, diffs, "k^--", ms=8, label="AFR-EAS diff")
ax.set_xlabel("Mutation rate (μ)"); ax.set_ylabel("Mean TMRCA (Ka)", color="black")
ax2.set_ylabel("AFR-EAS difference (Ka)", color="black")
ax.set_title("(a) E7: Mutation rate sensitivity")
ax.legend(fontsize=9, loc="upper left"); ax2.legend(fontsize=9, loc="upper right")
ax.set_xscale("log")

ax = axes[1]
subs = [s for s in AMR_SUBS if s in sub_res]
means= [sub_res[s]["mean"] for s in subs]
stds = [sub_res[s]["std"] for s in subs]
ns   = [sub_res[s]["n"] for s in subs]
y    = range(len(subs))
bars = ax.barh(y, means, xerr=[s/np.sqrt(n) for s,n in zip(stds,ns)],
               color=["#E91E63" if s in ("ASW","ACB") else "#FF9800" for s in subs],
               alpha=0.8, capsize=4)
ax.set_yticks(y); ax.set_yticklabels([f"{s}\n(N={n})" for s,n in zip(subs,ns)])
ax.axvline(afr_m, color="red", ls="--", lw=1.5, label=f"within-AFR ({results_mu[1]['afr']:.0f} Ka)")
ax.axvline(eas_m, color="blue", ls="--", lw=1.5, label=f"within-EAS ({results_mu[1]['eas']:.0f} Ka)")
ax.axvline(amr_all.mean(), color="orange", ls=":", lw=2, label=f"within-AMR all ({amr_all.mean():.0f} Ka)")
ax.set_xlabel("Within-subpopulation mean TMRCA (Ka)")
ax.set_title("(b) E8: AMR subpopulation breakdown\n(pink=African-admixed, orange=Latin American)")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("/home/yanlin/popghistory/results/validation/fig_e7e8.png", dpi=150, bbox_inches="tight")
print("\nSaved: results/validation/fig_e7e8.png")
