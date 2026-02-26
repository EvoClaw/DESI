"""
DESI Phase 4a Pilot (revised) — 10 kb windows for PCHMM
Reads chr22 VCF for 30 pilot samples, computes vectorized het-count matrix
in 10 kb windows, then runs PCHMM to get per-pair TMRCA estimates.
"""

import numpy as np
import subprocess, time
from scipy.special import gammaln

GEN_TIME = 28
MU       = 1.2e-8
WINDOW   = 10_000   # 10 kb

T_BINS_KA  = np.exp(np.linspace(np.log(10), np.log(5000), 40))
T_BINS_GEN = T_BINS_KA * 1000 / GEN_TIME

CHR22_LEN = 50_818_468
VCF_PATH  = ("/home/yanlin/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/"
             "1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

# 30 pilot samples (10 YRI + 10 CEU + 10 CHB) from 1kGP
SAMPLE_GROUPS = {
    "YRI": ["NA18484","NA18485","NA18486","NA18487","NA18488",
            "NA18489","NA18497","NA18498","NA18499","NA18500"],
    "CEU": ["NA06984","NA06985","NA06986","NA06989","NA06991",
            "NA06993","NA06994","NA06995","NA06997","NA07000"],
    "CHB": ["NA18525","NA18526","NA18528","NA18530","NA18531",
            "NA18532","NA18533","NA18534","NA18535","NA18536"],
}

def make_emission_log(t_bins_gen, window_bp, mu):
    max_het = int(max(t_bins_gen) * 2 * mu * window_bp * 2)
    max_het = min(max_het, 2000)
    rates = 2 * mu * window_bp * t_bins_gen
    n = np.arange(max_het + 1)
    lf = np.zeros(max_het + 1)
    for i in range(1, max_het + 1):
        lf[i] = lf[i-1] + np.log(i)
    log_E = (n[None,:] * np.log(rates[:,None] + 1e-300)
             - rates[:,None] - lf[None,:])
    return log_E, max_het

def pchmm_pair(het_counts_1d, log_E, t_bins_ka):
    """Independent-window PCHMM for one pair. Returns posterior mean TMRCA (Ka)."""
    K = len(t_bins_ka)
    log_pi = np.zeros(K)
    total = np.zeros(K)
    for n_het in het_counts_1d:
        nc = min(int(n_het), log_E.shape[1] - 1)
        lp = log_pi + log_E[:, nc]
        lZ = np.logaddexp.reduce(lp)
        total += lp - lZ
    avg = total - np.logaddexp.reduce(total)
    return float(np.dot(np.exp(avg), t_bins_ka))

def main():
    print("=== DESI Pilot 2: 10kb windows + PCHMM ===")
    t0 = time.time()

    # Select samples
    samples = []
    for pop, ids in SAMPLE_GROUPS.items():
        samples.extend(ids)
    n_samp = len(samples)
    n_hap  = n_samp * 2
    sample_str = ",".join(samples)
    print(f"Samples: {n_samp} ({','.join(SAMPLE_GROUPS.keys())}), {n_hap} haplotypes")

    # Compute number of windows
    n_windows = CHR22_LEN // WINDOW + 1
    print(f"chr22: {CHR22_LEN//1_000_000} Mb, {n_windows} windows of {WINDOW//1000}kb")

    # Read VCF using bcftools
    print(f"\nReading chr22 VCF...")
    cmd = (f"bcftools view -r chr22 -s {sample_str} -v snps "
           f"-m 2 -M 2 {VCF_PATH} | "
           f"bcftools query -f '[%GT\\t]%POS\\n'")

    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1<<22)

    # Per-pair het count matrix (stored as sparse: only accumulate per window)
    # For 1740 pairs × 5000 windows × 2 bytes = 17.4 MB
    n_pairs_max = n_hap * (n_hap - 1) // 2
    # Build pair index map
    rows, cols = np.triu_indices(n_hap, k=1)
    same_indiv = (rows // 2 == cols // 2)
    valid_mask = ~same_indiv
    valid_rows = rows[valid_mask]
    valid_cols = cols[valid_mask]
    n_pairs = len(valid_rows)

    # Map (i,j) to pair index
    pair_idx = {}
    for p, (r, c) in enumerate(zip(valid_rows, valid_cols)):
        pair_idx[(r, c)] = p

    het_matrix = np.zeros((n_pairs, n_windows), dtype=np.int16)

    n_snps = 0
    n_skipped = 0
    last_update = time.time()

    for line in proc.stdout:
        try:
            line = line.decode().strip()
            parts = line.split('\t')
            if len(parts) < n_samp + 1:
                n_skipped += 1
                continue
            pos = int(parts[-1])
            w = pos // WINDOW
            if w >= n_windows:
                continue

            # Parse haplotypes (all n_samp × 2 haplotypes)
            gts = np.empty(n_hap, dtype=np.int8)
            ok = True
            for si in range(n_samp):
                gt = parts[si]
                if '|' in gt:
                    a, b = gt.split('|')
                elif '/' in gt:
                    a, b = gt.split('/')
                else:
                    ok = False; break
                try:
                    gts[si*2] = int(a); gts[si*2+1] = int(b)
                except ValueError:
                    ok = False; break
            if not ok:
                n_skipped += 1; continue

            # Vectorized het update for this SNP and window
            # For each valid pair, check if they differ
            g_i = gts[valid_rows]  # (n_pairs,)
            g_j = gts[valid_cols]  # (n_pairs,)
            diff = (g_i != g_j).astype(np.int8)  # (n_pairs,)
            het_matrix[:, w] += diff

            n_snps += 1
            if time.time() - last_update > 30:
                print(f"  {n_snps//1000}k SNPs processed ({time.time()-t0:.0f}s)")
                last_update = time.time()
        except Exception as e:
            n_skipped += 1
            continue

    proc.wait()
    print(f"Done: {n_snps} SNPs, {n_skipped} skipped ({time.time()-t0:.0f}s)")

    # Save het matrix
    np.save("/tmp/desi_het_10kb_chr22.npy", het_matrix)
    print(f"Het matrix saved: {het_matrix.shape}, "
          f"mean het per window per pair = {het_matrix.mean():.2f}")

    # Distribution of window het counts
    flat = het_matrix.flatten()
    print(f"Het count distribution: 0:{(flat==0).mean()*100:.1f}% "
          f"1-5:{((flat>0)&(flat<=5)).mean()*100:.1f}% "
          f">5:{(flat>5).mean()*100:.1f}%")

    # Run PCHMM on all pairs
    print(f"\nRunning PCHMM for {n_pairs} pairs...")
    log_E, max_het = make_emission_log(T_BINS_GEN, WINDOW, MU)
    print(f"Emission matrix: {log_E.shape}, max_het={max_het}")

    t1 = time.time()
    pair_tmrca = np.zeros(n_pairs)
    for p_idx in range(n_pairs):
        pair_tmrca[p_idx] = pchmm_pair(het_matrix[p_idx, :], log_E, T_BINS_KA)
        if (p_idx+1) % 300 == 0:
            el = time.time()-t1
            rt = (n_pairs-p_idx-1)/(p_idx+1)*el
            print(f"  [{p_idx+1}/{n_pairs}] {el:.0f}s, ~{rt:.0f}s left")

    print(f"PCHMM complete: {n_pairs} pairs in {time.time()-t1:.0f}s")

    # Summary
    print(f"\nPer-pair TMRCA (Ka):")
    print(f"  Mean={np.mean(pair_tmrca):.0f}, Median={np.median(pair_tmrca):.0f}")
    print(f"  10pct={np.percentile(pair_tmrca,10):.0f}, 90pct={np.percentile(pair_tmrca,90):.0f}")
    print(f"  Min={np.min(pair_tmrca):.0f}, Max={np.max(pair_tmrca):.0f}")

    # Population stratification
    pop_labels = []
    for pop in SAMPLE_GROUPS:
        n_per_pop = len(SAMPLE_GROUPS[pop])
        pop_labels.extend([pop] * (n_per_pop * 2))
    pop_arr = np.array(pop_labels)
    hap_pop_r = pop_arr[valid_rows]
    hap_pop_c = pop_arr[valid_cols]
    for p1, p2 in [("YRI","YRI"),("CEU","CEU"),("CHB","CHB"),
                   ("YRI","CEU"),("YRI","CHB"),("CEU","CHB")]:
        pm = ((hap_pop_r==p1)&(hap_pop_c==p2)) | ((hap_pop_r==p2)&(hap_pop_c==p1))
        if pm.sum() > 0:
            print(f"  {p1}-{p2}: n={pm.sum()}, "
                  f"mean={np.mean(pair_tmrca[pm]):.0f} Ka, "
                  f"median={np.median(pair_tmrca[pm]):.0f} Ka")

    # K-test
    data = pair_tmrca[pair_tmrca > 10]

    def em_gm(d, K, ni=200):
        n=len(d); a=np.full(K,2.0)
        b=a/np.quantile(d,np.linspace(0.1,0.9,K)); w=np.ones(K)/K
        ld=np.log(d+1e-10); pl=-1e18
        for _ in range(ni):
            lr=np.zeros((n,K))
            for k in range(K): lr[:,k]=np.log(w[k]+1e-300)+(a[k]-1)*ld-b[k]*d+a[k]*np.log(b[k])-gammaln(a[k])
            lr-=lr.max(1,keepdims=True); r=np.exp(lr); r/=r.sum(1,keepdims=True)+1e-300
            Nk=r.sum(0); w=np.maximum(Nk/n,1e-10); w/=w.sum()
            for k in range(K):
                if Nk[k]<1: continue
                mk=np.dot(r[:,k],d)/Nk[k]; ml=np.dot(r[:,k],ld)/Nk[k]
                s=np.log(max(mk,1e-10))-ml
                if s<=0: s=1e-4
                a[k]=max(0.1,(3-s+np.sqrt((s-3)**2+24*s))/(12*s)); b[k]=a[k]/max(mk,1e-10)
            ll=sum(np.dot(r[:,k],(a[k]-1)*ld-b[k]*d)+Nk[k]*(np.log(w[k]+1e-300)+a[k]*np.log(b[k])-gammaln(a[k])) for k in range(K))
            if abs(ll-pl)<1e-7: break; pl=ll
        bic=-2*ll+(3*K-1)*np.log(n); return {'w':w,'means_ka':a/b,'ll':ll,'bic':bic}

    f1 = em_gm(data, 1); f2 = em_gm(data, 2)
    kh = 1 if f1['bic'] <= f2['bic'] else 2
    lbf = (f2['ll'] - f1['ll']) / np.log(10)
    idx = np.argsort(f2['means_ka'])

    print(f"\n=== DESI K-TEST (real chr22, PCHMM T estimates) ===")
    print(f"  K_hat={kh}, log10_BF={lbf:+.1f}, n={len(data)}")
    print(f"  K=1 mean: {f1['means_ka'][0]:.0f} Ka")
    print(f"  K=2 comp1: {f2['means_ka'][idx[0]]:.0f} Ka (w={f2['w'][idx[0]]:.3f})")
    print(f"  K=2 comp2: {f2['means_ka'][idx[1]]:.0f} Ka (w={f2['w'][idx[1]]:.3f})")
    print(f"  Separation: {f2['means_ka'][idx[1]]/max(1,f2['means_ka'][idx[0]]):.1f}x")

    np.save("/tmp/pchmm_pair_tmrca_chr22.npy", pair_tmrca)
    print(f"\nSaved /tmp/pchmm_pair_tmrca_chr22.npy")
    print(f"Total time: {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
