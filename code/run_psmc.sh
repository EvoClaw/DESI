#!/bin/bash
# Run PSMC on selected YRI and CHB individuals
# Goal: directly test if Ne(t) shows a bottleneck at ~930 ka

set -euo pipefail
export PATH=$PATH:/tmp/psmc:/tmp/psmc/utils

VCF_BASE="/home/yanlin/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.CHROM.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
OUTDIR="/home/yanlin/popghistory/results/psmc"
THREADS=8
mkdir -p "$OUTDIR"

CHR_LENGTHS="chr1:248956422 chr2:242193529 chr3:198295559 chr4:190214555 chr5:181538259 chr6:170805979 chr7:159345973 chr8:145138636 chr9:138394717 chr10:133797422 chr11:135086622 chr12:133275309 chr13:114364328 chr14:107043718 chr15:101991189 chr16:90338345 chr17:83257441 chr18:80373285 chr19:58617616 chr20:64444167 chr21:46709983 chr22:50818468"

BIN_SIZE=100  # 100bp bins for PSMC

run_sample() {
    local SAMPLE=$1
    local OUTFA="$OUTDIR/${SAMPLE}.psmcfa"
    local OUTPSMC="$OUTDIR/${SAMPLE}.psmc"

    if [ -f "$OUTPSMC" ]; then
        echo "[$SAMPLE] Already done, skipping"
        return 0
    fi

    echo "[$SAMPLE] Building psmcfa from VCF..."
    python3 - "$SAMPLE" "$OUTFA" "$VCF_BASE" "$BIN_SIZE" "$CHR_LENGTHS" << 'PYEOF'
import sys, subprocess, numpy as np

sample   = sys.argv[1]
outfa    = sys.argv[2]
vcf_base = sys.argv[3]
bin_size = int(sys.argv[4])
chr_str  = sys.argv[5]

chr_lens = {}
for item in chr_str.split():
    chrom, length = item.split(':')
    chr_lens[chrom] = int(length)

chroms = list(chr_lens.keys())

with open(outfa, 'w') as fout:
    for chrom in chroms:
        vcf = vcf_base.replace('CHROM', chrom)
        chrom_len = chr_lens[chrom]
        n_bins = (chrom_len // bin_size) + 1
        het_bins = np.zeros(n_bins, dtype=bool)

        cmd = ['bcftools', 'query', '-s', sample,
               '-f', '%POS\t[%GT]\n', '-r', chrom,
               '-i', 'TYPE="SNP"', vcf]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            for line in result.stdout.split('\n'):
                if not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                pos = int(parts[0])
                gt  = parts[1]
                # heterozygous: 0|1 or 1|0
                if '|' in gt:
                    alleles = gt.split('|')
                    if alleles[0] != alleles[1]:
                        bin_idx = (pos - 1) // bin_size
                        if 0 <= bin_idx < n_bins:
                            het_bins[bin_idx] = True
        except Exception as e:
            print(f"  Warning {chrom}: {e}", file=sys.stderr)
            continue

        seq = ''.join('T' if b else 'F' for b in het_bins)
        fout.write(f'>{chrom}\n{seq}\n')
        print(f"  {chrom}: {het_bins.sum()} het bins / {n_bins} total ({chrom_len/1e6:.0f} Mb)", file=sys.stderr)

print(f"psmcfa written: {outfa}", file=sys.stderr)
PYEOF

    if [ ! -f "$OUTFA" ]; then
        echo "[$SAMPLE] ERROR: psmcfa not created" >&2
        return 1
    fi

    echo "[$SAMPLE] Running PSMC..."
    psmc -N25 -t15 -r5 -p "4+25*2+4+6" "$OUTFA" > "$OUTPSMC"
    echo "[$SAMPLE] Done! Output: $OUTPSMC"
}

export -f run_sample
export VCF_BASE OUTDIR CHR_LENGTHS BIN_SIZE

# 4 YRI + 4 CHB samples
SAMPLES="NA18497 NA18498 NA18499 NA18517 NA18526 NA18535 NA18542 NA18548"

echo "$SAMPLES" | tr ' ' '\n' | xargs -P $THREADS -I{} bash -c 'run_sample "$@"' _ {}
echo "All done!"
