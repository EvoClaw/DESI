"""E3+E5: Structure fraction over time + per-population K-mixture"""
import numpy as np, os
from scipy.special import gammaln
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("/home/yanlin/popghistory/results/validation", exist_ok=True)

d = np.load("/home/yanlin/popghistory/results/final/desi_genome_tmrca.npy", allow_pickle=True).item()
pair_tmrca = d["pair_tmrca"]
hap_pop    = d["hap_pop"]
pair_rows  = d["pair_rows"]
pair_cols  = d["pair_cols"]

pop3 = {'YRI':'AFR','LWK':'AFR','GWD':'AFR','MSL':'AFR','ESN':'AFR',
        'ASW':'AMR','ACB':'AMR','MXL':'AMR','PUR':'AMR','CLM':'AMR','PEL':'AMR',
        'CHB':'EAS','JPT':'EAS','CHS':'EAS','CDX':'EAS','KHV':'EAS',
        'CEU':'EUR','TSI':'EUR','FIN':'EUR','GBR':'EUR','IBS':'EUR',
        'GIH':'SAS','PJL':'SAS','BEB':'SAS','STU':'SAS','ITU':'SAS'}
sup = np.array([pop3.get(p,'UNK') for p in hap_pop])
rs  = sup[pair_rows];  cs = sup[pair_cols]

def em_gamma(x, K, n_iter=400, n_restarts=10):
    n = len(x); best_ll=-1e18; best=None
    for seed in range(n_restarts):
        rng = np.random.default_rng(seed)
        mu0 = np.quantile(x, np.sort(rng.uniform(0.05,0.95,K)))
        var0 = np.var(x)*rng.uniform(0.3,3.0,K)
        sh = np.clip(mu0**2/(var0+1e-8),0.5,300.); sc = np.clip(var0/(mu0+1e-8),0.01,2000.)
        pi = np.ones(K)/K; ll_p=-1e18
        for _ in range(n_iter):
            lc = np.zeros((n,K))
            for k in range(K):
                lc[:,k] = ((sh[k]-1)*np.log(x+1e-300)-x/(sc[k]+1e-10)
                           -gammaln(sh[k])-sh[k]*np.log(sc[k]+1e-10)+np.log(pi[k]+1e-300))
            lZ = np.logaddexp.reduce(lc,axis=1,keepdims=True)
            resp = np.exp(lc-lZ); ll=float(lZ.sum())
            if ll-ll_p<1e-5: break
            ll_p=ll; Nk=resp.sum(0)+1e-10; pi=Nk/Nk.sum()
            muv=(resp*x[:,None]).sum(0)/Nk; ex2=(resp*x[:,None]**2).sum(0)/Nk
            v=np.clip(ex2-muv**2,1e-3,None)
            sh=np.clip(muv**2/(v+1e-8),0.3,500.); sc=np.clip(v/(muv+1e-8),0.01,5000.)
        if ll>best_ll: best_ll=ll; best=(pi.copy(),sh.copy(),sc.copy())
    pi,sh,sc=best; order=np.argsort(sh*sc); pi,sh,sc=pi[order],sh[order],sc[order]
    bic=-2*best_ll+(3*K-1)*np.log(n)
    return dict(ll=best_ll,bic=bic,means=sh*sc,pi=pi,sh=sh,sc=sc)

# ==== E5: Per-population K-mixture ====
print("=== E5: Per-population K-mixture ===")
POPS = ["AFR","EAS","EUR","SAS"]
km_res = {}
for sp in POPS:
    mask = (rs==sp) & (cs==sp)
    vals = pair_tmrca[mask]
    print(f"within-{sp}: N={len(vals):,} pairs, mean={vals.mean():.1f} Ka")
    res = {}
    for K in [1,2,3]:
        r = em_gamma(vals, K)
        res[K] = r
        print(f"  K={K}: BIC={r['bic']:.1f}, means={r['means'].round(1)}, pi={r['pi'].round(3)}")
    best_K = min(res, key=lambda k: res[k]['bic'])
    print(f"  Best K={best_K} (delta_BIC vs K=1: {res[best_K]['bic']-res[1]['bic']:.1f})")
    km_res[sp] = {"results":res, "best_K":best_K, "vals":vals}

# Also cross-group AFR-EAS
mask_ae = ((rs=="AFR")&(cs=="EAS"))|((rs=="EAS")&(cs=="AFR"))
vals_ae = pair_tmrca[mask_ae]
print(f"AFR-EAS cross: N={len(vals_ae):,}, mean={vals_ae.mean():.1f} Ka")
res_ae = {}
for K in [1,2]:
    r = em_gamma(vals_ae, K)
    res_ae[K] = r
    print(f"  K={K}: BIC={r['bic']:.1f}, means={r['means'].round(1)}")
km_res["AFR-EAS"] = {"results":res_ae, "best_K":min(res_ae,key=lambda k:res_ae[k]['bic']), "vals":vals_ae}

# ==== E3: Structure fraction over time ====
print("\n=== E3: Structure fraction over time ===")
# Bin pairs by TMRCA into 100 Ka windows, track population composition
bins = np.arange(100, 2601, 100)  # 100, 200, ... 2600 Ka
bin_idx = np.digitize(pair_tmrca, bins) - 1  # 0-indexed

struct_frac = []
for bi in range(len(bins)-1):
    mask = bin_idx == bi
    if mask.sum() < 20: struct_frac.append(np.nan); continue
    afr_in_bin  = ((rs[mask]=="AFR") | (cs[mask]=="AFR")).mean()  # fraction with at least one AFR
    pure_within = ((rs[mask]==cs[mask])).mean()  # within-group pairs
    # Fraction of pairs that are "non-AFR within-group" at this TMRCA depth
    eas_eur_within = (((rs[mask]!="AFR")&(cs[mask]!="AFR")&(rs[mask]==cs[mask]))).mean()
    struct_frac.append({"t_mid":(bins[bi]+bins[bi+1])/2, "n":mask.sum(),
                        "afr_frac":afr_in_bin, "within_frac":pure_within,
                        "eas_eur_within":eas_eur_within})

# AFR-only within pairs: what fraction of within-AFR pairs have deep coalescence?
afr_only = pair_tmrca[(rs=="AFR")&(cs=="AFR")]
print(f"within-AFR TMRCA: mean={afr_only.mean():.1f}, "
      f"frac>1500Ka: {(afr_only>1500).mean():.3f}, "
      f"frac>2000Ka: {(afr_only>2000).mean():.3f}")
eas_only = pair_tmrca[(rs=="EAS")&(cs=="EAS")]
print(f"within-EAS TMRCA: mean={eas_only.mean():.1f}, "
      f"frac>1500Ka: {(eas_only>1500).mean():.3f}, "
      f"frac>2000Ka: {(eas_only>2000).mean():.3f}")

# Save
np.save("/home/yanlin/popghistory/results/validation/e3e5_results.npy",
        {"km_res":{k:{kk:vv for kk,vv in v.items() if kk!="vals"} for k,v in km_res.items()},
         "struct_frac":struct_frac})

# ---- Figures ----
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("E3+E5: Per-population K-mixture & Structure Fraction", fontsize=12)

PCOLS = {"AFR":"#E91E63","EAS":"#2196F3","EUR":"#4CAF50","SAS":"#FF9800"}
ax = axes[0,0]
for sp in POPS:
    vals = km_res[sp]["vals"]
    ax.hist(vals, bins=100, range=(0,3000), density=True, alpha=0.45,
            color=PCOLS[sp], label=f"within-{sp} (mean={vals.mean():.0f}Ka)")
ax.set_xlabel("TMRCA (Ka)"); ax.set_ylabel("Density"); ax.set_title("(a) Within-group TMRCA distributions")
ax.legend(fontsize=8); ax.axvline(930, color="black", ls="--", lw=1, label="930 Ka")

ax = axes[0,1]
for si, sp in enumerate(POPS):
    res = km_res[sp]["results"]
    bics = [res[k]["bic"] for k in [1,2,3]]
    rel  = [b - bics[0] for b in bics]
    ax.plot([1,2,3], rel, "o-", color=PCOLS[sp], label=f"within-{sp}")
ax.axhline(0, color="grey", ls="--", lw=1)
ax.set_xlabel("K"); ax.set_ylabel("ΔBIC (relative to K=1)")
ax.set_title("(b) K-mixture BIC per population\n(negative=better than K=1)")
ax.legend(fontsize=8)

ax = axes[1,0]
ts_struct = [s["t_mid"] for s in struct_frac if isinstance(s,dict)]
afr_frac  = [s["afr_frac"] for s in struct_frac if isinstance(s,dict)]
eas_eur_w = [s["eas_eur_within"] for s in struct_frac if isinstance(s,dict)]
ax.plot(ts_struct, afr_frac,  "r-o", ms=4, label="Fraction of pairs with ≥1 AFR hap")
ax.plot(ts_struct, eas_eur_w, "b-s", ms=4, label="Fraction within-EAS/EUR (non-AFR)")
ax.axvline(930, color="black", ls="--", lw=1, label="930 Ka")
ax.axvline(1500, color="purple", ls="--", lw=1, label="1500 Ka (cobraa split)")
ax.set_xlabel("Mean TMRCA bin (Ka)"); ax.set_ylabel("Fraction of pairs")
ax.set_title("(c) Population composition vs TMRCA depth\n(E3: structure fraction over time)")
ax.legend(fontsize=8)

ax = axes[1,1]
bins100 = np.linspace(0, 3000, 61)
ax.hist(afr_only, bins=bins100, density=True, alpha=0.6, color="#E91E63", label=f"within-AFR (K={km_res['AFR']['best_K']})")
ax.hist(eas_only, bins=bins100, density=True, alpha=0.6, color="#2196F3", label=f"within-EAS (K={km_res['EAS']['best_K']})")
ax.hist(vals_ae,  bins=bins100, density=True, alpha=0.35, color="#9E9E9E", label=f"AFR-EAS cross (K={km_res['AFR-EAS']['best_K']})")
# Overlay K=2 fit for AFR
for sp, col in [("AFR","#E91E63"),("EAS","#2196F3")]:
    r2 = km_res[sp]["results"][2]
    x = np.linspace(10, 3000, 500)
    from scipy.stats import gamma as gamma_dist
    pdf = sum(r2["pi"][k]*gamma_dist.pdf(x, a=r2["sh"][k], scale=r2["sc"][k]) for k in range(2))
    ax.plot(x, pdf, color=col, lw=2, ls="--", label=f"{sp} K=2 fit")
ax.axvline(930, color="black", ls="--", lw=1)
ax.set_xlabel("TMRCA (Ka)"); ax.set_ylabel("Density")
ax.set_title("(d) AFR vs EAS distributions + K=2 fits")
ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig("/home/yanlin/popghistory/results/validation/fig_e3e5.png", dpi=150, bbox_inches="tight")
print("\nSaved: results/validation/fig_e3e5.png")
print("Saved: results/validation/e3e5_results.npy")
