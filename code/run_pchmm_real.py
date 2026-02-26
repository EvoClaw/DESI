"""
DESI PCHMM on real chr22 data (30 pilot samples: 10 YRI + 10 CEU + 10 CHB).
Reads VCF, computes per-pair TMRCA estimates using PCHMM, then runs K-test.
"""

import numpy as np
import subprocess, time, sys
from scipy.special import gammaln

GEN_TIME = 28
MU       = 1.2e-8
WINDOW   = 10_000      # 10 kb windows
RHO      = 1e-8

# T bins: 40 log-scale bins from 10 Ka to 5000 Ka
T_BINS_KA  = np.exp(np.linspace(np.log(10), np.log(5000), 40))
T_BINS_GEN = T_BINS_KA * 1000 / GEN_TIME

CHR22_LEN = 50_818_468  # bp
VCF_PATH  = "/home/yanlin/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/chr22.vcf.gz"

# 30 pilot samples (same as Phase 4a pilot)
SAMPLES_FILE = "/tmp/pilot_samples.txt"

def make_emission(t_bins_gen, window_bp, mu):
    max_het = int(max(t_bins_gen) * 2 * mu * window_bp * 2)
    max_het = min(max_het, 3000)
    rates = 2 * mu * window_bp * t_bins_gen
    n = np.arange(max_het + 1)
    log_factorials = np.zeros(len(n))
    for i in range(1, len(n)):
        log_factorials[i] = log_factorials[i-1] + np.log(i)
    log_E = (n[None,:] * np.log(rates[:,None] + 1e-300)
             - rates[:,None]
             - log_factorials[None,:])
    return log_E, max_het

def get_pair_tmrca_estimate(het_counts, log_E, t_bins_ka):
    K = len(t_bins_ka)
    log_pi = np.zeros(K)
    total_log_post = np.zeros(K)
    for n_het in het_counts:
        n_clipped = min(n_het, log_E.shape[1] - 1)
        log_post_w = log_pi + log_E[:, n_clipped]
        log_Z = np.logaddexp.reduce(log_post_w)
        total_log_post += log_post_w - log_Z
    avg_lp = total_log_post - np.logaddexp.reduce(total_log_post)
    post = np.exp(avg_lp)
    return float(np.dot(post, t_bins_ka))

def em_gamma_mixture(data_ka, K, n_iter=200, tol=1e-7):
    n = len(data_ka)
    a = np.full(K, 2.0)
    b = a / np.quantile(data_ka, np.linspace(0.1, 0.9, K))
    w = np.ones(K) / K
    ld = np.log(data_ka + 1e-10)
    prev_ll = -1e18
    for _ in range(n_iter):
        lr = np.zeros((n, K))
        for k in range(K):
            lr[:,k] = (np.log(w[k]+1e-300)+(a[k]-1)*ld-b[k]*data_ka
                       +a[k]*np.log(b[k])-gammaln(a[k]))
        lr -= lr.max(1, keepdims=True)
        r = np.exp(lr); r /= r.sum(1, keepdims=True)+1e-300
        Nk = r.sum(0); w = np.maximum(Nk/n, 1e-10); w /= w.sum()
        for k in range(K):
            if Nk[k] < 1: continue
            mk = np.dot(r[:,k], data_ka)/Nk[k]; ml = np.dot(r[:,k], ld)/Nk[k]
            s = np.log(max(mk,1e-10))-ml
            if s <= 0: s = 1e-4
            a[k] = max(0.1,(3-s+np.sqrt((s-3)**2+24*s))/(12*s)); b[k]=a[k]/max(mk,1e-10)
        ll = sum(np.dot(r[:,k],(a[k]-1)*ld-b[k]*data_ka)+Nk[k]*(np.log(w[k]+1e-300)+a[k]*np.log(b[k])-gammaln(a[k])) for k in range(K))
        if abs(ll-prev_ll)<tol: break
        prev_ll = ll
    bic = -2*ll+(3*K-1)*np.log(n)
    return {'w':w,'means_ka':a/b,'ll':ll,'bic':bic}

def main():
    print("DESI PCHMM on real chr22 data")
    print(f"Window={WINDOW//1000}kb, T_bins={len(T_BINS_KA)}, MU={MU}, RHO={RHO}")

    # Get sample list
    try:
        with open(SAMPLES_FILE) as f:
            samples = [l.strip() for l in f if l.strip()]
        print(f"Loaded {len(samples)} samples from {SAMPLES_FILE}")
    except FileNotFoundError:
        print(f"Sample file not found. Reading from VCF header...")
        cmd = f"bcftools query -l {VCF_PATH} | head -30"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        samples = result.stdout.strip().split("\n")[:30]
        with open(SAMPLES_FILE, 'w') as f:
            f.write("\n".join(samples) + "\n")
        print(f"Using {len(samples)} samples: {samples[:5]}...")

    n_samples = len(samples)
    sample_str = ",".join(samples)

    # Read VCF and extract genotypes
    print(f"\nReading VCF for chr22, {n_samples} samples...")
    t0 = time.time()

    cmd = (f"bcftools view -r chr22 -s {sample_str} -v snps "
           f"--min-af 0.0001 {VCF_PATH} | "
           f"bcftools query -f '[%GT\\t]%POS\\n'")

    print(f"Running: {cmd[:80]}...")
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1<<20)

    n_hap = n_samples * 2
    n_windows = CHR22_LEN // WINDOW + 1
    het_matrix = np.zeros((n_hap, n_hap, n_windows), dtype=np.int16)

    n_snps = 0
    for line in proc.stdout:
        line = line.decode().strip()
        parts = line.split('\t')
        if len(parts) < n_samples + 1:
            continue
        pos_str = parts[-1]
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        w = pos // WINDOW
        if w >= n_windows:
            continue

        # Parse haplotypes
        gts = []
        for gt in parts[:n_samples]:
            if '|' in gt:
                a, b = gt.split('|')
            elif '/' in gt:
                a, b = gt.split('/')
            else:
                continue
            try:
                gts.extend([int(a), int(b)])
            except ValueError:
                gts.extend([0, 0])

        if len(gts) != n_hap:
            continue

        gts_arr = np.array(gts, dtype=np.int8)
        # Vectorized het update
        for i in range(n_hap):
            if gts_arr[i] == 0:
                continue
            for j in range(i+1, n_hap):
                if gts_arr[j] != gts_arr[i]:
                    het_matrix[i, j, w] += 1
                    het_matrix[j, i, w] += 1

        n_snps += 1
        if n_snps % 100000 == 0:
            print(f"  {n_snps} SNPs processed ({time.time()-t0:.0f}s)...")

    proc.wait()
    print(f"Done reading: {n_snps} SNPs in {time.time()-t0:.0f}s")

    # Build emission matrix
    print("\nBuilding emission matrix...")
    log_E, max_het = make_emission(T_BINS_GEN, WINDOW, MU)
    print(f"  Emission matrix shape: {log_E.shape}")

    # Compute per-pair TMRCA estimates
    print("\nComputing per-pair TMRCA estimates using PCHMM...")
    rows, cols = np.triu_indices(n_hap, k=1)
    mask = ~(rows // 2 == cols // 2)
    i_pairs, j_pairs = rows[mask], cols[mask]
    n_pairs = len(i_pairs)
    print(f"  {n_pairs} cross-individual haplotype pairs")

    pair_tmrca = np.zeros(n_pairs)
    t1 = time.time()
    for p_idx, (i, j) in enumerate(zip(i_pairs, j_pairs)):
        het_counts = het_matrix[i, j, :]
        pair_tmrca[p_idx] = get_pair_tmrca_estimate(het_counts, log_E, T_BINS_KA)
        if (p_idx+1) % 200 == 0:
            elapsed = time.time() - t1
            rate = (p_idx+1) / elapsed
            remaining = (n_pairs - p_idx - 1) / rate
            print(f"  [{p_idx+1}/{n_pairs}] {elapsed:.0f}s, ~{remaining:.0f}s remaining")

    print(f"\nPer-pair TMRCA estimation complete: {n_pairs} pairs in {time.time()-t1:.0f}s")

    # Summary statistics
    print(f"\nPer-pair TMRCA distribution (Ka):")
    print(f"  Mean: {np.mean(pair_tmrca):.0f} Ka")
    print(f"  Median: {np.median(pair_tmrca):.0f} Ka")
    print(f"  10th pct: {np.percentile(pair_tmrca,10):.0f} Ka")
    print(f"  90th pct: {np.percentile(pair_tmrca,90):.0f} Ka")
    print(f"  Min: {np.min(pair_tmrca):.0f} Ka, Max: {np.max(pair_tmrca):.0f} Ka")

    # K-test
    data = pair_tmrca[pair_tmrca > 10]  # filter near-zero values
    print(f"\nK-test on {len(data)} pairs...")
    f1 = em_gamma_mixture(data, 1)
    f2 = em_gamma_mixture(data, 2)
    k_hat = 1 if f1['bic'] <= f2['bic'] else 2
    lbf = (f2['ll'] - f1['ll']) / np.log(10)
    idx = np.argsort(f2['means_ka'])

    print(f"\n=== DESI K-TEST RESULT (real chr22 data) ===")
    print(f"  n_pairs: {len(data)}")
    print(f"  K_hat: {k_hat}, log10_BF: {lbf:+.1f}")
    print(f"  K=1 fit: mean = {f1['means_ka'][0]:.0f} Ka")
    print(f"  K=2 fit: component 1 = {f2['means_ka'][idx[0]]:.0f} Ka (w={f2['w'][idx[0]]:.3f})")
    print(f"           component 2 = {f2['means_ka'][idx[1]]:.0f} Ka (w={f2['w'][idx[1]]:.3f})")
    print(f"  Separation ratio: {f2['means_ka'][idx[1]]/max(1,f2['means_ka'][idx[0]]):.1f}x")

    if k_hat == 2:
        print(f"\nINTERPRETATION: K=2 detected → ANCESTRAL STRUCTURE signal")
        print(f"  Split time indicator: ~{f2['means_ka'][idx[1]]:.0f} Ka (deeper component)")
        print(f"  Structure fraction: {f2['w'][idx[1]]:.3f}")
    else:
        print(f"\nINTERPRETATION: K=1 → consistent with panmixia at deep time")
        print(f"  (or insufficient signal with only 30 pilot samples)")

    # Save results
    np.save("/tmp/pchmm_pair_tmrca_chr22.npy", pair_tmrca)
    print(f"\nResults saved to /tmp/pchmm_pair_tmrca_chr22.npy")

if __name__ == "__main__":
    main()
