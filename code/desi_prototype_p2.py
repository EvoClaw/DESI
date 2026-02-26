"""DESI Prototype Part 2: K-component Gamma mixture inference per time bin."""
import numpy as np
from scipy.special import gammaln
import time

GEN_TIME = 28
SEED = 42
np.random.seed(SEED)

pairwise_T = np.load("/tmp/desi_pairwise_T.npy")
print(f"Loaded {len(pairwise_T):,} pairwise T estimates")
print(f"T range (Ka): {pairwise_T.min()*GEN_TIME/1000:.1f} -- {pairwise_T.max()*GEN_TIME/1000:.0f}")

bins_ka = np.array([100, 200, 300, 500, 700, 930, 1250, 1500, 2000])
bins_gen = bins_ka * 1000 / GEN_TIME

def fit_gamma_mixture_em(data, K, n_iter=300, tol=1e-7):
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
                 + Nk[k]*(np.log(weights[k]+1e-300) + alphas[k]*np.log(betas[k]) - gammaln(alphas[k]))
                 for k in range(K))
        if abs(ll - prev_ll) < tol: break
        prev_ll = ll
    K_params = 3*K - 1
    bic = -2*ll + K_params*np.log(n)
    return {'weights': weights, 'alphas': alphas, 'betas': betas,
            'Ne_k': alphas/betas/2, 'means_ka': alphas/betas*GEN_TIME/1000,
            'loglik': ll, 'bic': bic, 'n_obs': n}

print(f"\n{'Ka':>6} {'n_obs':>8} {'K_hat':>6} {'w_minor':>8} {'Ne_major':>10} {'Ne_minor':>10} {'logBF(2v1)':>11}")
print("-" * 70)

summary = {}
for ka, center in zip(bins_ka, bins_gen):
    lo, hi = center * 0.70, center * 1.30
    data = pairwise_T[(pairwise_T >= lo) & (pairwise_T <= hi)]
    if len(data) < 50:
        print(f"{ka:6d} {len(data):8d}   (too few)")
        continue
    f1 = fit_gamma_mixture_em(data, K=1)
    f2 = fit_gamma_mixture_em(data, K=2)
    if f1 is None or f2 is None:
        continue
    k_hat = 1 if f1['bic'] <= f2['bic'] else 2
    lbf = (f2['loglik'] - f1['loglik']) / np.log(10)
    if k_hat == 2:
        idx = np.argsort(f2['Ne_k'])
        w_min = f2['weights'][idx[0]]
        ne_maj = f2['Ne_k'][idx[1]]
        ne_min = f2['Ne_k'][idx[0]]
    else:
        w_min = 0.0; ne_maj = f1['Ne_k'][0]; ne_min = float('nan')
    summary[ka] = {'k_hat': k_hat, 'w_minor': w_min,
                   'Ne_major': ne_maj, 'Ne_minor': ne_min, 'lbf': lbf, 'n': len(data)}
    print(f"{ka:6d} {len(data):8,} {k_hat:6d} {w_min:8.3f} {ne_maj:10.0f} {ne_min:10.0f} {lbf:11.2f}")

import json
with open("/tmp/desi_summary.json", "w") as f:
    json.dump({int(k): {kk: (float(v) if not isinstance(v, str) else v)
                        for kk, v in vv.items()} for k,vv in summary.items()}, f, indent=2)

print("\nKey result at 930 Ka:")
if 930 in summary:
    r = summary[930]
    print(f"  K_hat={r['k_hat']}, w_minor={r['w_minor']:.3f}, log10_BF(K=2 vs K=1)={r['lbf']:.2f}")
    if r['k_hat'] == 2:
        print(f"  Ne_major={r['Ne_major']:.0f}, Ne_minor={r['Ne_minor']:.0f}")
        print(f"  --> STRUCTURE SIGNAL DETECTED at 930 Ka (consistent with Cobraa)")
    else:
        print(f"  Ne={r['Ne_major']:.0f}")
        print(f"  --> Panmictic signal at 930 Ka (consistent with FitCoal bottleneck model)")
