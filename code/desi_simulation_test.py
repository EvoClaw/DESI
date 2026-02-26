"""
Phase 4a: Validate DESI on msprime simulations.
Test if DESI can DETECT structure when it truly exists (power test).
Run 2 models x 5 reps each:
  Model A: panmictic + bottleneck at 930 Ka (FitCoal params)
  Model D: two-pop structured + bottleneck + admix 300 Ka (Cobraa params)
"""
import numpy as np
import msprime
from scipy.special import gammaln
import time

GEN_TIME = 28
MU = 1.2e-8
RHO = 1.0e-8
WINDOW_BP = 1_000_000
SEQ_LEN = 50_000_000  # 50 Mb simulation (much faster than full chr22)
N_SAMPLES = 30  # 30 diploid individuals
SEED = 42
np.random.seed(SEED)

def fit_gamma_mixture_em(data, K, n_iter=200, tol=1e-7):
    n = len(data)
    if n < K * 20:
        return None
    quantiles = np.linspace(0.1, 0.9, K)
    means = np.quantile(data, quantiles)
    alphas = np.full(K, 2.0)
    betas = alphas / means
    weights = np.ones(K) / K
    log_data = np.log(data + 1e-10)
    prev_ll = -np.inf
    for it in range(n_iter):
        log_resp = np.zeros((n, K))
        for k in range(K):
            log_resp[:, k] = (np.log(weights[k]+1e-300)
                              + (alphas[k]-1)*log_data
                              - betas[k]*data
                              + alphas[k]*np.log(betas[k])
                              - gammaln(alphas[k]))
        lrm = log_resp.max(axis=1, keepdims=True)
        log_resp -= lrm
        resp = np.exp(log_resp)
        resp /= resp.sum(axis=1, keepdims=True) + 1e-300
        Nk = resp.sum(axis=0)
        weights = np.maximum(Nk / n, 1e-10)
        weights /= weights.sum()
        for k in range(K):
            if Nk[k] < 1:
                continue
            wk = resp[:, k]
            mean_k = np.dot(wk, data) / Nk[k]
            mean_log_k = np.dot(wk, log_data) / Nk[k]
            s = np.log(max(mean_k, 1e-10)) - mean_log_k
            if s <= 0: s = 1e-4
            alphas[k] = max(0.1, (3-s+np.sqrt((s-3)**2+24*s))/(12*s))
            betas[k] = alphas[k] / max(mean_k, 1e-10)
        ll = sum(np.dot(resp[:,k], (alphas[k]-1)*log_data - betas[k]*data)
                 + Nk[k]*(np.log(weights[k]+1e-300)+alphas[k]*np.log(betas[k])-gammaln(alphas[k]))
                 for k in range(K))
        if abs(ll - prev_ll) < tol: break
        prev_ll = ll
    K_params = 3*K - 1
    bic = -2*ll + K_params*np.log(n)
    return {'weights': weights, 'Ne_k': alphas/betas/2,
            'loglik': ll, 'bic': bic, 'n_obs': n}

def simulate_and_analyze(demography, model_name, seed, n_samples=N_SAMPLES):
    """Simulate, compute pairwise TMRCA, run DESI, return K_hat at 930 Ka."""
    ts = msprime.sim_ancestry(
        samples=n_samples,
        demography=demography,
        sequence_length=SEQ_LEN,
        recombination_rate=RHO,
        random_seed=seed
    )
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)

    # Extract genotype matrix per 1Mb window
    n_haplo = n_samples * 2
    windows = {}
    for var in ts.variants():
        pos = int(var.site.position)
        win_id = pos // WINDOW_BP
        gt = var.genotypes  # shape: (n_haplo,)
        if win_id not in windows:
            windows[win_id] = []
        windows[win_id].append(gt.astype(np.float32))

    # Pairwise T per window
    pT = []
    for win_id, geno_list in windows.items():
        if len(geno_list) < 5:
            continue
        geno = np.array(geno_list)  # (n_snps, n_haplo)
        col_sum = geno.sum(axis=0)
        dot_mat = geno.T @ geno
        het_mat = col_sum[:, None] + col_sum[None, :] - 2*dot_mat
        T_mat = het_mat / (2 * MU * WINDOW_BP)
        for i in range(n_haplo):
            for j in range(i+1, n_haplo):
                if i // 2 == j // 2:
                    continue
                if T_mat[i,j] > 0:
                    pT.append(float(T_mat[i,j]))

    if not pT:
        return None
    pT = np.array(pT)

    # Time bin at 930 Ka
    center = 930 * 1000 / GEN_TIME
    lo, hi = center * 0.70, center * 1.30
    data = pT[(pT >= lo) & (pT <= hi)]

    if len(data) < 50:
        return {'k_hat': None, 'n_obs': len(data), 'w_minor': None, 'lbf': None}

    f1 = fit_gamma_mixture_em(data, K=1)
    f2 = fit_gamma_mixture_em(data, K=2)
    if f1 is None or f2 is None:
        return None

    k_hat = 1 if f1['bic'] <= f2['bic'] else 2
    lbf = (f2['loglik'] - f1['loglik']) / np.log(10)
    if k_hat == 2:
        idx = np.argsort(f2['Ne_k'])
        w_min = f2['weights'][idx[0]]
    else:
        w_min = 0.0
    return {'k_hat': k_hat, 'w_minor': w_min, 'lbf': lbf, 'n_obs': len(data)}


# ── Model A: panmictic + bottleneck 930 Ka (FitCoal parameters) ──────────────
def model_A_demography():
    d = msprime.Demography()
    d.add_population(name="ANC", initial_size=10000)
    # Bottleneck: Ne collapses to 1280 from 930 Ka to 813 Ka
    # Times in generations (generation time = 28 yr)
    t_bottle_start = 930000 / GEN_TIME  # ~33214 gen
    t_bottle_end   = 813000 / GEN_TIME  # ~29036 gen
    d.add_population_parameters_change(time=t_bottle_end,  population="ANC", initial_size=1280)
    d.add_population_parameters_change(time=t_bottle_start, population="ANC", initial_size=10000)
    return d

# ── Model D: 2-pop structured + major lineage bottleneck + admix 300 Ka ──────
def model_D_demography():
    d = msprime.Demography()
    t_split  = 1500000 / GEN_TIME  # ~53571 gen
    t_admix  = 300000  / GEN_TIME  # ~10714 gen
    t_bottle = 930000  / GEN_TIME  # ~33214 gen
    # Two ancestral populations
    d.add_population(name="pop1", initial_size=8000)   # major (80%)
    d.add_population(name="pop2", initial_size=2000)   # minor (20%)
    # Admixture event at 300 Ka: pop2 merges 20% into pop1
    d.add_admixture(time=t_admix, derived="pop1",
                    ancestral=["pop1", "pop2"], proportions=[0.80, 0.20])
    # Bottleneck in pop1 at 930 Ka
    d.add_population_parameters_change(time=t_bottle * 0.87, population="pop1", initial_size=1280)
    d.add_population_parameters_change(time=t_bottle,        population="pop1", initial_size=8000)
    # Split at 1.5 Ma
    d.add_population_split(time=t_split, derived=["pop1", "pop2"], ancestral="pop1")
    return d


print("=" * 60)
print("DESI Simulation Power Test — Phase 4a")
print(f"N_samples={N_SAMPLES}, SeqLen={SEQ_LEN/1e6:.0f}Mb, n_reps=5 per model")
print("=" * 60)

models = [("Model_A_panmictic_bottleneck", model_A_demography),
          ("Model_D_structured_bottleneck", model_D_demography)]

all_results = {}
for model_name, dem_fn in models:
    print(f"\n--- {model_name} ---")
    results = []
    for rep in range(5):
        seed = SEED + rep * 100
        t0 = time.time()
        try:
            d = dem_fn()
            r = simulate_and_analyze(d, model_name, seed)
            dt = time.time()-t0
            if r:
                print(f"  Rep {rep+1}: K_hat={r['k_hat']}, w_minor={r.get('w_minor',0):.3f}, logBF={r.get('lbf',0):.2f}, n={r['n_obs']}, t={dt:.1f}s")
            else:
                print(f"  Rep {rep+1}: no result, t={dt:.1f}s")
        except Exception as e:
            print(f"  Rep {rep+1}: ERROR: {e}")
            r = None
        results.append(r)
    all_results[model_name] = results

    # Recovery stats
    k_hats = [r['k_hat'] for r in results if r and r['k_hat'] is not None]
    correct_k = 1 if "panmictic" in model_name else 2
    if k_hats:
        recovery = sum(1 for k in k_hats if k == correct_k) / len(k_hats)
        print(f"  K recovery rate (K={correct_k}): {recovery*100:.0f}% ({sum(1 for k in k_hats if k==correct_k)}/{len(k_hats)})")

print("\nSummary:")
print(f"  Model A (panmictic): DESI should output K=1 (correct = panmictic)")
print(f"  Model D (structured): DESI should output K=2 (correct = structured)")
