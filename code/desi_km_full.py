"""
Full K-mixture analysis (K=1..5) + window-level Ne(t) for representative pairs.
Two outputs:
  /tmp/km_full_results.npy   - BIC/means/pi for K=1..5
  /tmp/ne_t_results.npy      - per-window TMRCA for Ne(t) analysis
"""
import numpy as np, glob, subprocess, time
from scipy.special import gammaln

MU        = 1.2e-8
WINDOW    = 10_000
GEN_TIME  = 28
T_BINS_KA  = np.exp(np.linspace(np.log(10), np.log(5000), 40))
T_BINS_GEN = T_BINS_KA * 1000 / GEN_TIME
K_BINS     = len(T_BINS_KA)

POP_MAP = {
    'YRI':'AFR','GWD':'AFR','ESN':'AFR','MSL':'AFR','ACB':'AFR',
    'CEU':'EUR','TSI':'EUR','CHB':'EAS','JPT':'EAS','CHS':'EAS',
    'PJL':'SAS','GIH':'SAS','PUR':'AMR','PEL':'AMR',
}
VCF_BASE = ("/home/yanlin/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/"
            "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

# ── K-mixture EM ─────────────────────────────────────────────────────────────
def em_gamma(x, K, n_iter=500, n_restarts=15):
    n=len(x); best_ll=-1e18; best=None
    for seed in range(n_restarts):
        rng=np.random.default_rng(seed)
        mu0=np.quantile(x, np.sort(rng.uniform(0.05,0.95,K)))
        var0=np.var(x)*rng.uniform(0.3,3.0,K)
        sh=np.clip(mu0**2/(var0+1e-8),0.5,300.); sc=np.clip(var0/(mu0+1e-8),0.01,2000.)
        pi=np.ones(K)/K; ll_p=-1e18
        for _ in range(n_iter):
            lc=np.zeros((n,K))
            for k in range(K):
                lc[:,k]=((sh[k]-1)*np.log(x+1e-300)-x/(sc[k]+1e-10)
                         -gammaln(sh[k])-sh[k]*np.log(sc[k]+1e-10)+np.log(pi[k]+1e-300))
            lZ=np.logaddexp.reduce(lc,axis=1,keepdims=True); resp=np.exp(lc-lZ); ll=float(lZ.sum())
            if ll-ll_p<1e-5: break
            ll_p=ll; Nk=resp.sum(0)+1e-10; pi=Nk/Nk.sum()
            mu_v=(resp*x[:,None]).sum(0)/Nk; ex2=(resp*x[:,None]**2).sum(0)/Nk
            v=np.clip(ex2-mu_v**2,1e-3,None)
            sh=np.clip(mu_v**2/(v+1e-8),0.3,500.); sc=np.clip(v/(mu_v+1e-8),0.01,5000.)
        if ll>best_ll: best_ll=ll; best=(pi.copy(),sh.copy(),sc.copy())
    pi,sh,sc=best; order=np.argsort(sh*sc); pi,sh,sc=pi[order],sh[order],sc[order]
    return dict(ll=best_ll, bic=-2*best_ll+(3*K-1)*np.log(n), means=sh*sc, pi=pi, sh=sh, sc=sc)

# ── Part 1: K-mixture BIC scan ────────────────────────────────────────────────
print("=== Part 1: K-mixture BIC scan ===", flush=True)
files=sorted(glob.glob('/tmp/desi_pchmm_chr*.npy'))
per_chr=[]; hap_pop=None; rows0=cols0=None; samples=None
for fn in files:
    d=np.load(fn,allow_pickle=True).item()
    per_chr.append(d['pair_tmrca'])
    if hap_pop is None:
        hap_pop=list(d['hap_pop']); rows0=d['pair_rows']; cols0=d['pair_cols']
        samples=list(d['samples'])

tmrca_pool = np.concatenate(per_chr)
np.random.seed(99); idx=np.random.choice(len(tmrca_pool), 500000, replace=False)
x = tmrca_pool[idx]
print(f"  Fitting K=1..5 on {len(x):,} pooled data points...", flush=True)
km_results={}
for K in range(1,6):
    km_results[K]=em_gamma(x,K)
    m=km_results[K]
    print(f"  K={K}: BIC={m['bic']:.0f}  means={np.round(m['means']).astype(int)}  pi={np.round(m['pi'],3)}", flush=True)
best_K=min(km_results,key=lambda k:km_results[k]['bic'])
print(f"  Best K={best_K}  dBIC(2vs1)={km_results[2]['bic']-km_results[1]['bic']:.0f}  dBIC(3vs2)={km_results[3]['bic']-km_results[2]['bic']:.0f}", flush=True)

# ── Part 2: Window-level Ne(t) for representative pairs ───────────────────────
print("\n=== Part 2: Window-level Ne(t) analysis ===", flush=True)

# Build emission table for window-level
def make_emission(t_bins_gen, window_bp, mu):
    max_het = min(int(max(t_bins_gen) * 2 * mu * window_bp * 2), 2000)
    rates   = 2 * mu * window_bp * t_bins_gen
    n       = np.arange(max_het + 1)
    lf      = np.zeros(max_het + 1)
    for i in range(1, max_het + 1): lf[i] = lf[i-1] + np.log(i)
    return n[None,:]*np.log(rates[:,None]+1e-300) - rates[:,None] - lf[None,:], max_het

log_E, max_het_e = make_emission(T_BINS_GEN, WINDOW, MU)

# Select 4 representative pairs per continental group
n_hap = len(hap_pop)
rep_pairs = {}  # group -> list of (hap_i, hap_j, sample_i, sample_j)
for group in ['EAS','EUR','AFR','SAS']:
    haps = [i for i,p in enumerate(hap_pop) if POP_MAP.get(p,'OTH')==group]
    if len(haps) >= 4:
        chosen = [(haps[i], haps[i+2]) for i in range(0, min(8, len(haps)-2), 4)][:4]
        rep_pairs[group] = [(hi, hj) for hi,hj in chosen if hi//2 != hj//2][:3]
        print(f"  {group}: {len(rep_pairs[group])} pairs, e.g. hap{rep_pairs[group][0][0]}({hap_pop[rep_pairs[group][0][0]]}) vs hap{rep_pairs[group][0][1]}({hap_pop[rep_pairs[group][0][1]]})", flush=True)

# AFR-EAS cross-population pairs
afr_haps = [i for i,p in enumerate(hap_pop) if POP_MAP.get(p,'OTH')=='AFR'][:4]
eas_haps = [i for i,p in enumerate(hap_pop) if POP_MAP.get(p,'OTH')=='EAS'][:4]
rep_pairs['AFR_EAS'] = [(afr_haps[i], eas_haps[i]) for i in range(min(3, len(afr_haps), len(eas_haps)))]

# For each group, get per-window TMRCA on all autosomes
ne_t_data = {}  # group -> array of per-window T values (Ka)
chroms = [f'chr{i}' for i in range(1,23)]

for group, pairs in rep_pairs.items():
    print(f"\n  Processing {group} ({len(pairs)} pairs, all autosomes)...", flush=True)
    all_win_t = []
    for hap_i, hap_j in pairs:
        samp_i = samples[hap_i//2]
        samp_j = samples[hap_j//2]
        hap_idx_i = hap_i % 2  # 0 or 1
        hap_idx_j = hap_j % 2
        win_ts = []
        for chrom in chroms:
            vcf = VCF_BASE.format(chrom=chrom)
            if samp_i == samp_j:
                sample_str = samp_i
            else:
                sample_str = f"{samp_i},{samp_j}"
            cmd = (f"bcftools view -r {chrom} -s {sample_str} -v snps -m 2 -M 2 {vcf} | "
                   f"bcftools query -f '[%GT\\t]%POS\\n'")
            n_samp_local = 1 if samp_i == samp_j else 2
            GT_BYTES = n_samp_local * 4
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1<<22)
            cur_win = -1; win_cnt = 0; win_lp = np.zeros(K_BINS)
            for raw in proc.stdout:
                raw = raw.rstrip(b'\n')
                if len(raw) < GT_BYTES + 1: continue
                try: pos = int(raw[GT_BYTES:])
                except ValueError: continue
                w = pos // WINDOW
                arr = np.frombuffer(raw, dtype=np.uint8, count=GT_BYTES)
                # get the two specific alleles
                if samp_i == samp_j:  # same sample, two haplotypes
                    a_i = int(arr[hap_idx_i * 2]) - 48
                    a_j = int(arr[hap_idx_j * 2]) - 48
                else:
                    a_i = int(arr[hap_idx_i * 2]) - 48       # allele from sample i
                    a_j = int(arr[4 + hap_idx_j * 2]) - 48   # allele from sample j
                if a_i < 0 or a_i > 1 or a_j < 0 or a_j > 1: continue
                het = int(a_i != a_j)
                if w != cur_win:
                    if cur_win >= 0:
                        lZ = np.logaddexp.reduce(win_lp)
                        post = np.exp(win_lp - lZ)
                        win_ts.append(float(post @ T_BINS_KA))
                    win_lp = np.zeros(K_BINS); cur_win = w; win_cnt = 0
                h = min(het, max_het_e)  # single-site het (0 or 1)
                win_lp += log_E[:, h]
                win_cnt += 1
            if cur_win >= 0:
                lZ = np.logaddexp.reduce(win_lp)
                post = np.exp(win_lp - lZ)
                win_ts.append(float(post @ T_BINS_KA))
            proc.wait()
            print(f"    {chrom}: {len(win_ts)} windows so far", flush=True)
        all_win_t.extend(win_ts)
    ne_t_data[group] = np.array(all_win_t)
    print(f"  {group}: {len(ne_t_data[group])} window estimates, mean={np.mean(ne_t_data[group]):.0f} Ka", flush=True)

np.save('/tmp/ne_t_results.npy', ne_t_data, allow_pickle=True)
print("\nSaved /tmp/ne_t_results.npy", flush=True)

# ── Save K-mixture results ────────────────────────────────────────────────────
np.save('/tmp/km_full_results.npy', km_results, allow_pickle=True)
print("Saved /tmp/km_full_results.npy", flush=True)
print("DONE ALL", flush=True)
