"""DESI Prototype Part 1: Parse VCF, compute pairwise het counts per 1Mb window."""
import numpy as np
import subprocess, time, sys

MU = 1.2e-8
GEN_TIME = 28
WINDOW_BP = 1_000_000
SEED = 42
np.random.seed(SEED)

VCF = "/tmp/pilot_chr22.vcf.gz"
SAMPLE_FILE = "/tmp/pilot_samples.txt"

with open(SAMPLE_FILE) as f:
    samples = [s.strip() for s in f.readlines()]
n_samples = len(samples)
print(f"Samples: {n_samples}, haplotypes: {n_samples*2}", flush=True)

cmd = f"bcftools query -f '[%GT\\t]%POS\\n' {VCF}"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, bufsize=1<<20)

win_buffers = {}
n_variants = 0
t0 = time.time()
for line in proc.stdout:
    parts = line.decode().rstrip().split('\t')
    pos = int(parts[-1])
    win_id = (pos - 1) // WINDOW_BP
    gts_raw = parts[:-1]
    haplo = []
    for gt in gts_raw:
        if '|' in gt:
            a, b = gt.split('|')
            haplo.append(int(a) if a.isdigit() else 0)
            haplo.append(int(b) if b.isdigit() else 0)
        else:
            haplo.extend([0, 0])
    if win_id not in win_buffers:
        win_buffers[win_id] = []
    win_buffers[win_id].append(haplo)
    n_variants += 1
    if n_variants % 100000 == 0:
        print(f"  {n_variants:,} variants...", flush=True)
proc.wait()
print(f"Total variants: {n_variants:,} in {time.time()-t0:.1f}s")

# Compute pairwise TMRCA per window
n_haplo = n_samples * 2
pairwise_T = []
for win_id in sorted(win_buffers.keys()):
    geno = np.array(win_buffers[win_id], dtype=np.float32)  # (n_snps, n_haplo)
    if geno.shape[0] < 10:
        continue
    col_sum = geno.sum(axis=0)
    dot_mat = geno.T @ geno
    outer_sum = col_sum[:, None] + col_sum[None, :]
    het_mat = outer_sum - 2 * dot_mat
    T_mat = het_mat / (2 * MU * WINDOW_BP)
    for i in range(n_haplo):
        for j in range(i+1, n_haplo):
            if i // 2 == j // 2:
                continue
            T_val = T_mat[i, j]
            if T_val > 0:
                pairwise_T.append(float(T_val))

pairwise_T = np.array(pairwise_T)
np.save("/tmp/desi_pairwise_T.npy", pairwise_T)
print(f"Saved {len(pairwise_T):,} pairwise T estimates")
print(f"T range: {pairwise_T.min():.0f} -- {pairwise_T.max():.0f} gen")
print(f"T range (Ka): {pairwise_T.min()*28/1000:.0f} -- {pairwise_T.max()*28/1000:.0f} Ka")
