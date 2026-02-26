"""
E5 (corrected): Per-population K-mixture on per-WINDOW TMRCA values.

Uses the window-level TMRCA estimates (win_t, shape n_windows×15) from all
22 chromosomes. These have ~600 Ka std, suitable for Gamma mixture fitting.

Key question: Does within-AFR window distribution show K=2 bimodality
(indicating ancient sub-groups within Africa), while within-EAS does not?
"""
import numpy as np, glob, os
from scipy.special import gammaln
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("/home/yanlin/popghistory/results/validation", exist_ok=True)

# Load all 22 chromosomes
files = sorted(glob.glob("/tmp/desi_pchmm_chr*.npy"))
print(f"Loading {len(files)} chromosomes ...")

WIN_T_ALL = []  # will be (total_windows, 15)
for fn in files:
    d = np.load(fn, allow_pickle=True).item()
    WIN_T_ALL.append(d["win_t"])
WIN_T = np.concatenate(WIN_T_ALL, axis=0)  # (total_windows, 15)
GROUP_NAMES = np.load(files[0], allow_pickle=True).item()["group_names"]
N_WIN = WIN_T.shape[0]
print(f"Total windows: {N_WIN}  Groups: {GROUP_NAMES[:5]}")

def em_gamma(x, K, n_iter=500, n_restarts=15):
    """Gamma mixture EM — no shape cap; handles broad distributions."""
    x = x[~np.isnan(x) & (x > 1)]  # valid values only
    n = len(x)
    if n < 20: return None
    best_ll = -1e18; best = None
    for seed in range(n_restarts):
        rng = np.random.default_rng(seed)
        quantiles = np.sort(rng.uniform(0.05, 0.95, K))
        mu0  = np.quantile(x, quantiles)
        var0 = np.var(x) * rng.uniform(0.3, 3.0, K)
        sh   = mu0**2 / (var0 + 1e-8)
        sc   = var0 / (mu0 + 1e-8)
        pi   = np.ones(K) / K
        ll_p = -1e18
        for _ in range(n_iter):
            lc = np.zeros((n, K))
            for k in range(K):
                lc[:,k] = ((sh[k]-1)*np.log(x+1e-300)
                           - x/(sc[k]+1e-10)
                           - gammaln(sh[k])
                           - sh[k]*np.log(sc[k]+1e-10)
                           + np.log(pi[k]+1e-300))
            lZ = np.logaddexp.reduce(lc, axis=1, keepdims=True)
            resp = np.exp(lc - lZ)
            ll = float(lZ.sum())
            if ll - ll_p < 1e-4: break
            ll_p = ll
            Nk = resp.sum(0) + 1e-10
            pi = Nk / Nk.sum()
            mu_v = (resp * x[:,None]).sum(0) / Nk
            ex2  = (resp * x[:,None]**2).sum(0) / Nk
            v    = np.clip(ex2 - mu_v**2, 1e-1, None)
            sh   = mu_v**2 / (v + 1e-8)
            sc   = v / (mu_v + 1e-8)
        if ll > best_ll:
            best_ll = ll
            best = (pi.copy(), sh.copy(), sc.copy())
    if best is None: return None
    pi, sh, sc = best
    order = np.argsort(sh * sc)
    pi, sh, sc = pi[order], sh[order], sc[order]
    n_par = 3*K - 1
    bic = -2*best_ll + n_par * np.log(n)
    return dict(ll=best_ll, bic=bic, means=sh*sc, pi=pi, sh=sh, sc=sc, n=n)

GROUPS_TO_TEST = [0,1,2,3,4,5]  # within-AFR,EAS,EUR,AMR,SAS, AFR-EAS cross
km_res = {}
print("\n=== E5: Per-population K-mixture on per-window TMRCA ===")
for gi in GROUPS_TO_TEST:
    gname = GROUP_NAMES[gi]
    vals = WIN_T[:, gi]
    ok_vals = vals[~np.isnan(vals) & (vals > 10)]
    print(f"\nGroup: {gname}  N_windows={len(ok_vals):,}, "
          f"mean={ok_vals.mean():.1f}, std={ok_vals.std():.1f} Ka")
    res = {}
    for K in [1, 2, 3]:
        r = em_gamma(ok_vals, K)
        if r:
            res[K] = r
            print(f"  K={K}: BIC={r['bic']:.1f}, means={r['means'].round(1)}, "
                  f"pi={r['pi'].round(3)}, n={r['n']}")
        else:
            print(f"  K={K}: fit failed")
    if res:
        best_K = min(res, key=lambda k: res[k]['bic'])
        delta = {k: res[k]['bic'] - res[1]['bic'] for k in res if k > 1}
        print(f"  Best K={best_K}  ΔBIC vs K=1: {delta}")
    km_res[gname] = {"results": res, "vals": ok_vals,
                     "best_K": best_K if res else 1}

# Save (no large arrays)
save_km = {}
for g, v in km_res.items():
    save_km[g] = {"results": v["results"],
                  "best_K": v["best_K"],
                  "n_wins": len(v["vals"]),
                  "mean": v["vals"].mean(),
                  "std":  v["vals"].std()}
np.save("/home/yanlin/popghistory/results/validation/e5_results.npy", save_km)

# ---- Figures ----
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle("E5: Per-population K-mixture on per-window TMRCA (all 22 chr)", fontsize=12)
PCOLS = {"within-AFR":"#E91E63","within-EAS":"#2196F3","within-EUR":"#4CAF50",
         "within-AMR":"#FF9800","within-SAS":"#9C27B0","AFR-EAS":"#607D8B"}

# (a) Distributions
ax = axes[0,0]
for gname, col in list(PCOLS.items())[:5]:
    if gname in km_res:
        ax.hist(km_res[gname]["vals"], bins=80, range=(50,3000),
                density=True, alpha=0.45, color=col,
                label=f"{gname} (n={len(km_res[gname]['vals']):,})")
ax.axvline(930, color="black", ls="--", lw=1.5, label="930 Ka (FitCoal)")
ax.axvline(1500, color="purple", ls="--", lw=1, label="1500 Ka (cobraa)")
ax.set_xlabel("Per-window TMRCA (Ka)"); ax.set_ylabel("Density")
ax.set_title("(a) Per-window TMRCA distributions\n(all 22 autosomes, 100kb windows)")
ax.legend(fontsize=7)

# (b) BIC comparison
ax = axes[0,1]
for gname, col in PCOLS.items():
    if gname in km_res and km_res[gname]["results"]:
        bics = [km_res[gname]["results"].get(k, {}).get("bic", np.nan) for k in [1,2,3]]
        rel  = [b - bics[0] for b in bics]
        ax.plot([1,2,3], rel, "o-", color=col, label=gname)
ax.axhline(0, color="grey", ls="--", lw=1)
ax.set_xlabel("K"); ax.set_ylabel("ΔBIC (vs K=1)")
ax.set_title("(b) BIC by K  (negative = better than K=1)")
ax.legend(fontsize=7)

# (c) AFR K=2 decomposition
ax = axes[1,0]
from scipy.stats import gamma as gamma_dist
gname = "within-AFR"
if gname in km_res and 2 in km_res[gname]["results"]:
    r2 = km_res[gname]["results"][2]
    x  = np.linspace(10, 3500, 500)
    total_pdf = sum(r2["pi"][k]*gamma_dist.pdf(x, a=r2["sh"][k], scale=r2["sc"][k]) for k in range(2))
    ax.hist(km_res[gname]["vals"], bins=80, range=(50,3500), density=True, alpha=0.5,
            color="#E91E63", label="within-AFR data")
    ax.plot(x, total_pdf, "r-", lw=2, label=f"K=2 fit (total)")
    for k in range(2):
        comp_pdf = r2["pi"][k]*gamma_dist.pdf(x, a=r2["sh"][k], scale=r2["sc"][k])
        ax.plot(x, comp_pdf, "--", lw=1.5, label=f"Comp {k+1}: μ={r2['means'][k]:.0f}Ka, π={r2['pi'][k]:.2f}")
    ax.axvline(930, color="black", ls="--", lw=1)
    ax.set_title(f"(c) within-AFR K=2 decomposition")
else:
    ax.text(0.5,0.5,"K=2 fit not available", ha="center", va="center", transform=ax.transAxes)
ax.set_xlabel("Per-window TMRCA (Ka)"); ax.set_ylabel("Density"); ax.legend(fontsize=7)

# (d) EAS K=2 decomposition
ax = axes[1,1]
gname = "within-EAS"
if gname in km_res and 2 in km_res[gname]["results"]:
    r2 = km_res[gname]["results"][2]
    x  = np.linspace(10, 3500, 500)
    total_pdf = sum(r2["pi"][k]*gamma_dist.pdf(x, a=r2["sh"][k], scale=r2["sc"][k]) for k in range(2))
    ax.hist(km_res[gname]["vals"], bins=80, range=(50,3500), density=True, alpha=0.5,
            color="#2196F3", label="within-EAS data")
    ax.plot(x, total_pdf, "b-", lw=2, label="K=2 fit (total)")
    for k in range(2):
        comp_pdf = r2["pi"][k]*gamma_dist.pdf(x, a=r2["sh"][k], scale=r2["sc"][k])
        ax.plot(x, comp_pdf, "--", lw=1.5, label=f"Comp {k+1}: μ={r2['means'][k]:.0f}Ka, π={r2['pi'][k]:.2f}")
    ax.axvline(930, color="black", ls="--", lw=1)
    ax.set_title(f"(d) within-EAS K=2 decomposition")
else:
    ax.text(0.5,0.5,"K=2 fit not available", ha="center", va="center", transform=ax.transAxes)
ax.set_xlabel("Per-window TMRCA (Ka)"); ax.set_ylabel("Density"); ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig("/home/yanlin/popghistory/results/validation/fig_e5.png", dpi=150, bbox_inches="tight")
print("\nSaved: results/validation/fig_e5.png")
print("Saved: results/validation/e5_results.npy")
