"""
DESI genomic scan v3 — per-group aggregate TMRCA per 100 kb window.

Core change from v2:
  Instead of storing a per-pair TMRCA histogram, each 100 kb window
  yields ONE T estimate per population group, obtained by treating the
  total het-count across ALL pairs in the group as a single Poisson
  observation.  This reduces measurement uncertainty from ±329 Ka
  (10 kb, 1 pair) to <5 Ka (100 kb, thousands of pairs), comparable
  to FitCoal's formal error while remaining free of its panmixia
  assumption.

Outputs per chromosome  →  /tmp/desi_pchmm_<chrom>.npy
  win_t    (n_win, N_GROUPS)  – posterior-mean T per window per group (Ka)
                                  NaN when group has too few pairs
  win_n    (n_win, N_GROUPS)  – number of valid pairs per window per group
  win_pos  (n_win,)           – genomic position of each window (bp)
  win_blk  (n_win,)           – 1 Mb block id (for jackknife)
  pair_tmrca (n_pairs,)       – genome-wide posterior-mean TMRCA per pair (Ka)

Design is genome-agnostic: the statistical unit is the GENOMIC BLOCK
(1 Mb), not the chromosome.  Block jackknife provides proper SE for
any diploid species regardless of karyotype.

Called as:  python3 desi_run_chr.py <chrom>
"""

import numpy as np
import subprocess, time, sys

GEN_TIME   = 28
MU         = 1.2e-8
WINDOW     = 100_000          # 100 kb — reduces per-window noise to <5 Ka
BLOCK_SIZE = 1_000_000        # 1 Mb blocks for jackknife
MIN_PAIRS  = 10               # minimum pairs for a group estimate in a window

T_BINS_KA  = np.exp(np.linspace(np.log(10), np.log(5000), 60))
T_BINS_GEN = T_BINS_KA * 1000 / GEN_TIME
K          = len(T_BINS_KA)
BASE_RATE  = 2 * MU * WINDOW * T_BINS_GEN   # shape (K,)  [expected het/pair/window]

AFRICA_POPS = {'YRI','LWK','GWD','MSL','ESN','ASW','ACB','MKK','MAG'}
EAS_POPS    = {'CHB','JPT','CHS','CDX','KHV'}
EUR_POPS    = {'CEU','TSI','FIN','GBR','IBS'}
AMR_POPS    = {'MXL','PUR','CLM','PEL'}
SAS_POPS    = {'GIH','PJL','BEB','STU','ITU'}

GROUP_NAMES = [
    'within-AFR', 'within-EAS', 'within-EUR', 'within-AMR', 'within-SAS',
    'AFR-EAS', 'AFR-EUR', 'AFR-AMR', 'AFR-SAS',
    'EAS-EUR', 'EAS-AMR', 'EAS-SAS',
    'EUR-AMR', 'EUR-SAS',
    'AMR-SAS',
]
N_GROUPS = len(GROUP_NAMES)

_CROSS = {
    (0,1):5,(0,2):6,(0,3):7,(0,4):8,
    (1,2):9,(1,3):10,(1,4):11,
    (2,3):12,(2,4):13,
    (3,4):14,
}

VCF_BASE = ("/home/yanlin/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/"
            "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
SAMPLES_FILE  = "/home/yanlin/popghistory/docs/04_data_resource/desi_240_samples.txt"
POP_INFO_FILE = "/home/yanlin/public/1000GP/samples.info"

CHR_LENGTHS = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}

MAX_WIN_SNPS = 10_000   # cap SNPs buffered per window (memory safety)


def hap_region_idx(pop):
    if pop in AFRICA_POPS: return 0
    if pop in EAS_POPS:    return 1
    if pop in EUR_POPS:    return 2
    if pop in AMR_POPS:    return 3
    if pop in SAS_POPS:    return 4
    return -1


def build_pair_groups(hap_pop, rows, cols):
    hap_reg = np.array([hap_region_idx(p) for p in hap_pop], dtype=np.int8)
    rr, cc  = hap_reg[rows], hap_reg[cols]
    groups  = np.full(len(rows), -1, dtype=np.int8)
    for i, (r, c) in enumerate(zip(rr.tolist(), cc.tolist())):
        if r < 0 or c < 0: continue
        if r == c:
            groups[i] = r
        else:
            gi = _CROSS.get((min(r,c), max(r,c)), -1)
            if gi >= 0: groups[i] = gi
    return groups


def make_emission_log(t_bins_gen, window_bp, mu):
    """Log-emission matrix for per-pair genome-wide posterior (unchanged from v2)."""
    max_het = int(max(t_bins_gen) * 2 * mu * window_bp * 2)
    max_het = min(max_het, 3000)
    rates   = 2 * mu * window_bp * t_bins_gen
    n       = np.arange(max_het + 1)
    lf      = np.zeros(max_het + 1)
    for i in range(1, max_het + 1): lf[i] = lf[i-1] + np.log(i)
    return n[None,:]*np.log(rates[:,None]+1e-300) - rates[:,None] - lf[None,:], max_het


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 desi_run_chr.py chrN"); sys.exit(1)

    chrom    = sys.argv[1]
    out_file = f"/tmp/desi_pchmm_{chrom}.npy"
    t0       = time.time()
    print(f"=== DESI v3 (100kb aggregate): {chrom} ===", flush=True)

    with open(SAMPLES_FILE) as f:
        samples = [l.strip() for l in f if l.strip()]
    n_samp     = len(samples)
    n_hap      = n_samp * 2
    sample_str = ",".join(samples)
    print(f"Samples: {n_samp} ({n_hap} haplotypes)", flush=True)

    pop_info = {}
    try:
        with open(POP_INFO_FILE) as f:
            for ln in f:
                p = ln.strip().split('\t')
                if len(p) >= 2: pop_info[p[0]] = p[1]
    except FileNotFoundError:
        pass
    pop_labels = [pop_info.get(s, "UNK") for s in samples]
    hap_pop    = [pop_labels[i // 2] for i in range(n_hap)]

    chr_len     = CHR_LENGTHS.get(chrom, 250_000_000)
    n_win_total = chr_len // WINDOW + 1
    print(f"{chrom}: {chr_len//1_000_000} Mb → {n_win_total} windows @ {WINDOW//1000} kb each",
          flush=True)

    rows, cols = np.triu_indices(n_hap, k=1)
    valid       = ~(rows // 2 == cols // 2)
    valid_rows  = rows[valid].astype(np.int32)
    valid_cols  = cols[valid].astype(np.int32)
    n_pairs     = len(valid_rows)
    print(f"Cross-individual pairs: {n_pairs:,}", flush=True)

    print("Building emission matrix...", flush=True)
    log_E, max_het_e = make_emission_log(T_BINS_GEN, WINDOW, MU)

    print("Assigning pair groups...", flush=True)
    pair_group  = build_pair_groups(hap_pop, valid_rows, valid_cols)
    valid_gp    = pair_group >= 0
    pg_valid    = pair_group[valid_gp].astype(np.int32)

    for gi, gn in enumerate(GROUP_NAMES):
        cnt = (pair_group == gi).sum()
        if cnt > 0:
            print(f"  {gn:20s}: {cnt:,} pairs", flush=True)

    # ── streaming buffers ─────────────────────────────────────────────────────
    total_lp = np.zeros((n_pairs, K), dtype=np.float64)
    G_buf    = np.zeros((MAX_WIN_SNPS, n_hap), dtype=np.float32)

    # Accumulators for per-window results (filled in flush_window)
    win_t_list   = []   # each entry: (N_GROUPS,) float64, NaN for low-n groups
    win_n_list   = []   # each entry: (N_GROUPS,) int32
    win_pos_list = []   # each entry: int (bp start of window)
    win_blk_list = []   # each entry: int (1 Mb block id)

    win_cnt    = 0
    cur_win    = -1
    n_snps     = 0
    n_win_seen = 0
    last_log   = time.time()

    GT_BYTES = n_samp * 4
    A1_OFFS  = np.arange(0, GT_BYTES, 4, dtype=np.int32)
    A2_OFFS  = np.arange(2, GT_BYTES, 4, dtype=np.int32)

    def flush_window(G_w, win_pos):
        nonlocal n_win_seen, total_lp

        # ── pairwise het counts ──────────────────────────────────────────────
        rs  = G_w.sum(axis=0)
        GtG = G_w.T @ G_w
        hf  = rs[:, None] + rs[None, :] - 2.0 * GtG
        hp  = np.clip(hf[valid_rows, valid_cols].astype(np.int32), 0, max_het_e)

        # ── genome-wide per-pair posterior (for pair_tmrca) ──────────────────
        lp = log_E[:, hp].T                                 # (n_pairs, K)
        lZ = np.logaddexp.reduce(lp, axis=1, keepdims=True)
        total_lp += lp - lZ

        # ── per-group aggregate T estimate ───────────────────────────────────
        hp_valid  = hp[valid_gp].astype(np.float64)        # het for valid-pop pairs

        # Sum het counts and count pairs per group
        sum_het_g = np.bincount(pg_valid, weights=hp_valid, minlength=N_GROUPS)
        n_g       = np.bincount(pg_valid,                   minlength=N_GROUPS).astype(np.float64)

        # Group-level Poisson log-likelihood over T bins:
        #   log p(sum_het | T) = sum_het * log(n_g * base_rate) - n_g * base_rate
        rates_g   = n_g[:, None] * BASE_RATE[None, :]      # (N_GROUPS, K)
        log_ll    = (sum_het_g[:, None]
                     * np.log(rates_g + 1e-300)
                     - rates_g)                             # (N_GROUPS, K)

        # Normalise (uniform prior over log-spaced bins)
        log_ll   -= np.logaddexp.reduce(log_ll, axis=1, keepdims=True)
        post      = np.exp(log_ll)                          # (N_GROUPS, K)
        T_est     = post @ T_BINS_KA                        # (N_GROUPS,)

        # Mask groups below minimum pair count
        T_est[n_g < MIN_PAIRS] = np.nan

        win_t_list.append(T_est.copy())
        win_n_list.append(n_g.copy())
        win_pos_list.append(win_pos)
        win_blk_list.append(win_pos // BLOCK_SIZE)
        n_win_seen += 1

    # ── main VCF loop ─────────────────────────────────────────────────────────
    vcf = VCF_BASE.format(chrom=chrom)
    cmd = (f"bcftools view -r {chrom} -s {sample_str} -v snps -m 2 -M 2 {vcf} | "
           f"bcftools query -f '[%GT\\t]%POS\\n'")
    print("Reading VCF...", flush=True)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1 << 23)

    for raw in proc.stdout:
        raw = raw.rstrip(b'\n')
        if len(raw) < GT_BYTES + 1: continue
        try:
            pos = int(raw[GT_BYTES:])
        except ValueError:
            continue
        w = pos // WINDOW
        if w >= n_win_total: continue

        arr = np.frombuffer(raw, dtype=np.uint8, count=GT_BYTES)
        a1  = arr[A1_OFFS].astype(np.int16) - 48
        a2  = arr[A2_OFFS].astype(np.int16) - 48
        if a1.min() < 0 or a1.max() > 1 or a2.min() < 0 or a2.max() > 1: continue

        gts       = np.empty(n_hap, dtype=np.float32)
        gts[0::2] = a1;  gts[1::2] = a2

        if w != cur_win:
            if cur_win >= 0 and win_cnt > 0:
                flush_window(G_buf[:win_cnt], cur_win * WINDOW)
            win_cnt = 0
            cur_win = w

        if win_cnt < MAX_WIN_SNPS:
            G_buf[win_cnt] = gts
            win_cnt += 1
        n_snps += 1

        if time.time() - last_log > 60:
            print(f"  {n_snps//1000:,}k SNPs  {n_win_seen} windows  "
                  f"{time.time()-t0:.0f}s", flush=True)
            last_log = time.time()

    if cur_win >= 0 and win_cnt > 0:
        flush_window(G_buf[:win_cnt], cur_win * WINDOW)
    proc.wait()

    print(f"SNPs: {n_snps:,}  windows: {n_win_seen}  elapsed: {time.time()-t0:.0f}s",
          flush=True)

    # ── genome-wide pair_tmrca ─────────────────────────────────────────────────
    log_norm   = np.logaddexp.reduce(total_lp, axis=1, keepdims=True)
    post       = np.exp(total_lp - log_norm)
    pair_tmrca = post @ T_BINS_KA
    print(f"pair_tmrca: mean={np.mean(pair_tmrca):.0f} Ka  "
          f"med={np.median(pair_tmrca):.0f} Ka", flush=True)

    # ── assemble window arrays ─────────────────────────────────────────────────
    win_t   = np.array(win_t_list,   dtype=np.float64)   # (n_win, N_GROUPS)
    win_n   = np.array(win_n_list,   dtype=np.float64)   # (n_win, N_GROUPS)
    win_pos = np.array(win_pos_list, dtype=np.int64)      # (n_win,)
    win_blk = np.array(win_blk_list, dtype=np.int64)      # (n_win,)

    print(f"Window arrays: {win_t.shape}", flush=True)
    for gi, gn in enumerate(GROUP_NAMES):
        valid_mask = ~np.isnan(win_t[:, gi])
        if valid_mask.sum() > 0:
            t_vals = win_t[valid_mask, gi]
            print(f"  {gn:20s}: {valid_mask.sum()} windows  "
                  f"mean={np.mean(t_vals):.0f} Ka  "
                  f"median={np.median(t_vals):.0f} Ka  "
                  f"p10={np.percentile(t_vals,10):.0f}  "
                  f"p90={np.percentile(t_vals,90):.0f}", flush=True)

    np.save(out_file, {
        'chrom':       chrom,
        'win_t':       win_t,
        'win_n':       win_n,
        'win_pos':     win_pos,
        'win_blk':     win_blk,
        'pair_tmrca':  pair_tmrca,
        'pair_rows':   valid_rows,
        'pair_cols':   valid_cols,
        'samples':     np.array(samples),
        'hap_pop':     np.array(hap_pop),
        't_bins_ka':   T_BINS_KA,
        'group_names': GROUP_NAMES,
        'window_bp':   WINDOW,
        'block_size':  BLOCK_SIZE,
        'n_snps':      n_snps,
    }, allow_pickle=True)
    print(f"Saved → {out_file}   total: {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
