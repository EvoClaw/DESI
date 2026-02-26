"""
Baseline: SFS-based demographic inference using dadi.
Computes the folded SFS from 1kGP data for representative populations,
then fits a series of size-change models to infer Ne(t).

Used to compare with DESI's TMRCA-based demographic inference.
Usage: python3 baseline_sfs.py --pop YRI --chrom chr22
"""

import numpy as np
import argparse, os, subprocess, sys, time
import dadi

# ── Config ────────────────────────────────────────────────────────────────────

VCF_BASE = ("/home/yanlin/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/"
            "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
SAMPLES_FILE = "/home/yanlin/popghistory/docs/04_data_resource/desi_240_samples.txt"
POP_INFO_FILE = "/home/yanlin/public/1000GP/samples.info"

MU         = 1.2e-8     # mutation rate per bp per generation
GEN_TIME   = 28          # generation time in years
OUT_DIR    = "/home/yanlin/popghistory/results"


# ── Demographic models for dadi ────────────────────────────────────────────────

def model_panmictic_3epoch(params, ns, pts):
    """Three-epoch size change model: anc → intermediate → current."""
    nu1, nu2, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T1, nu1)
    phi = dadi.Integration.one_pop(phi, xx, T2, nu2)
    fs  = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs


def model_bottleneck(params, ns, pts):
    """Bottleneck model: ancestral pop → bottleneck → expansion → current."""
    nu_anc, nu_bot, nu_cur, T_bot_dur, T_rec = params
    xx  = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T_rec, nu_cur)   # expansion
    # Note: dadi integrates backwards, so this is actually:
    # ancient → bottleneck → recovery → current (if T_rec large)
    phi = dadi.Integration.one_pop(phi, xx, T_bot_dur, nu_bot)  # bottleneck
    phi = dadi.Integration.one_pop(phi, xx, 100.0, nu_anc)      # deep past
    fs  = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs


def model_bottleneck_simple(params, ns, pts):
    """Simple exponential-then-flat model to capture deep bottleneck."""
    nu1, nu2, T1, T2, T3 = params
    xx  = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.Integration.one_pop(phi, xx, T3, 1.0)   # current
    phi = dadi.Integration.one_pop(phi, xx, T2, nu2)    # post-bottleneck
    phi = dadi.Integration.one_pop(phi, xx, T1, nu1)    # bottleneck
    fs  = dadi.Spectrum.from_phi(phi, ns, (xx,))
    return fs


# ── SFS extraction from VCF ────────────────────────────────────────────────────

def extract_sfs_from_vcf(chrom, pop_samples, out_dir):
    """Use bcftools to count allele frequencies per SNP, build folded SFS."""
    sample_str = ",".join(pop_samples)
    n          = len(pop_samples)
    vcf        = VCF_BASE.format(chrom=chrom)

    cmd = (f"bcftools view -r {chrom} -s {sample_str} -v snps -m 2 -M 2 {vcf} | "
           f"bcftools query -f '%AC\\n' -S <(echo '{sample_str.replace(',', chr(10))}')")
    # Simpler: use bcftools stats for allele count
    cmd = (f"bcftools view -r {chrom} -s {sample_str} -v snps -m 2 -M 2 {vcf} | "
           f"bcftools query -f '[%GT]\\n' | "
           f"awk '{{c=0; for(i=1;i<=length($0);i++) if(substr($0,i,1)==\"1\") c++; print c}}'")

    print(f"Extracting SFS for {len(pop_samples)} samples on {chrom}...", flush=True)
    t0 = time.time()

    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"bcftools error: {proc.stderr[:200]}", flush=True)
        return None

    counts = np.array([int(x) for x in proc.stdout.strip().split('\n') if x], dtype=int)
    n_hap  = n * 2
    sfs    = np.zeros(n_hap + 1, dtype=float)
    for c in counts:
        if 0 < c < n_hap:
            sfs[c] += 1

    # Fold
    folded_n = n_hap // 2
    sfs_fold = np.zeros(folded_n + 1, dtype=float)
    for i in range(n_hap + 1):
        j = min(i, n_hap - i)
        sfs_fold[j] += sfs[i]
    sfs_fold[0] = 0.0   # remove monomorphic

    print(f"SFS done: {len(counts)} SNPs in {time.time()-t0:.0f}s", flush=True)

    sfs_out = f"{out_dir}/sfs_{chrom}.npy"
    np.save(sfs_out, sfs_fold)
    return sfs_fold


# ── dadi fitting ──────────────────────────────────────────────────────────────

def fit_dadi(sfs_fold, n_hap, mu, gen_time, out_dir, pop):
    """Fit a simple 2-epoch model to the observed SFS using dadi."""
    print("Fitting 2-epoch model with dadi...", flush=True)

    # Convert SFS to dadi Spectrum object (folded)
    ns    = (n_hap,)
    data  = dadi.Spectrum(sfs_fold, mask=[True] + [False]*(len(sfs_fold)-2) + [True],
                          folded=True)
    pts_l = [n_hap, n_hap + 5, n_hap + 10]

    # 2-epoch model: nu, T (relative to N_ref)
    def two_epoch(params, ns, pts):
        nu, T = params
        xx  = dadi.Numerics.default_grid(pts)
        phi = dadi.PhiManip.phi_1D(xx)
        phi = dadi.Integration.one_pop(phi, xx, T, nu)
        return dadi.Spectrum.from_phi(phi, ns, (xx,))

    func_ex = dadi.Numerics.make_extrap_log_func(two_epoch)

    # Optimize
    best_ll = -np.inf
    best_p  = None
    for _ in range(10):
        p0 = [np.random.uniform(0.01, 5.0), np.random.uniform(0.1, 10.0)]
        try:
            p_opt = dadi.Inference.optimize_log_fmin(
                p0, data, func_ex, pts_l,
                lower_bound=[0.001, 0.01], upper_bound=[100, 100],
                verbose=False, maxiter=100)
            model = func_ex(p_opt, ns, pts_l)
            ll    = dadi.Inference.ll_multinom(model, data)
            if ll > best_ll:
                best_ll = ll
                best_p  = p_opt
        except Exception as e:
            continue

    if best_p is None:
        print("dadi fitting failed", flush=True)
        return None

    nu_opt, T_opt = best_p
    theta_opt = dadi.Inference.optimal_sfs_scaling(
        func_ex(best_p, ns, pts_l), data)

    # Convert to physical units
    N_ref   = theta_opt / (4 * mu * 1e6)   # N_ref in individuals (assuming 1Mb used)
    T_years = T_opt * 2 * N_ref * gen_time

    print(f"Best 2-epoch fit: nu={nu_opt:.3f}  T={T_opt:.3f}  "
          f"N_ref={N_ref:.0f}  T={T_years/1000:.0f} Ka", flush=True)
    print(f"Log-likelihood: {best_ll:.2f}", flush=True)

    result = dict(nu=nu_opt, T=T_opt, N_ref=N_ref, T_years=T_years,
                  ll=best_ll, theta=theta_opt)
    np.save(f"{out_dir}/dadi_fit_{pop}.npy", result, allow_pickle=True)
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pop',   default='YRI', help='Population code')
    parser.add_argument('--chrom', default='chr22')
    parser.add_argument('--n',     type=int, default=20, help='Samples per pop')
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    # Load pop info
    pop_info = {}
    with open(POP_INFO_FILE) as f:
        for ln in f:
            p = ln.strip().split('\t')
            if len(p) >= 2: pop_info[p[0]] = p[1]

    # Load all DESI samples and filter by pop
    with open(SAMPLES_FILE) as f:
        all_samples = [l.strip() for l in f if l.strip()]

    pop_samples = [s for s in all_samples if pop_info.get(s, '') == args.pop]
    if not pop_samples:
        # Fallback: use first n_samples from Africa/EUR
        for target_pop in ['YRI', 'CEU', 'CHB']:
            pop_samples = [s for s in all_samples if pop_info.get(s, '') == target_pop]
            if pop_samples:
                break

    pop_samples = pop_samples[:args.n]
    print(f"Pop: {args.pop}  Samples: {len(pop_samples)}", flush=True)

    # Extract SFS
    sfs = extract_sfs_from_vcf(args.chrom, pop_samples, OUT_DIR)
    if sfs is None:
        print("SFS extraction failed", flush=True)
        sys.exit(1)

    # Fit dadi model
    n_hap = len(pop_samples) * 2
    fit_dadi(sfs, n_hap, MU, GEN_TIME, OUT_DIR, f"{args.pop}_{args.chrom}")


if __name__ == "__main__":
    main()
