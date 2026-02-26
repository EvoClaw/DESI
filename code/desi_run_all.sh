#!/bin/bash
# DESI Phase 4b: Run PCHMM on all 22 autosomes in parallel.
# Uses at most MAX_JOBS concurrent jobs.
# Usage: bash desi_run_all.sh [MAX_JOBS]

MAX_JOBS=${1:-22}   # default: all 22 in parallel (we have 128 threads)
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20
        chr21 chr22)
SCRIPT="$(dirname "$0")/desi_run_chr.py"
LOG_DIR="/tmp/desi_logs"
OUT_DIR="/tmp"

mkdir -p "$LOG_DIR"
echo "=== DESI full-genome run: $(date) ==="
echo "Chromosomes: ${CHROMS[*]}"
echo "Max parallel jobs: $MAX_JOBS"
echo ""

# Skip chromosomes that already have results
pids=()
for chrom in "${CHROMS[@]}"; do
    out="/tmp/desi_pchmm_${chrom}.npy"
    if [ -f "$out" ]; then
        echo "SKIP $chrom (result exists: $out)"
        continue
    fi

    log="$LOG_DIR/desi_${chrom}.log"
    echo "START $chrom → $log"
    python3 -u "$SCRIPT" "$chrom" > "$log" 2>&1 &
    pids+=("$!")

    # Throttle to MAX_JOBS concurrent
    while [ ${#pids[@]} -ge $MAX_JOBS ]; do
        # Wait for any job to finish
        new_pids=()
        for pid in "${pids[@]}"; do
            if kill -0 "$pid" 2>/dev/null; then
                new_pids+=("$pid")
            fi
        done
        pids=("${new_pids[@]}")
        [ ${#pids[@]} -ge $MAX_JOBS ] && sleep 10
    done
done

# Wait for all remaining jobs
echo ""
echo "All jobs launched. Waiting for completion..."
for pid in "${pids[@]}"; do
    wait "$pid"
done

echo ""
echo "=== All chromosomes done: $(date) ==="
echo ""
echo "Results:"
for chrom in "${CHROMS[@]}"; do
    out="/tmp/desi_pchmm_${chrom}.npy"
    if [ -f "$out" ]; then
        size=$(du -h "$out" | cut -f1)
        echo "  OK  $chrom  $size"
    else
        echo "  FAIL $chrom (no output)"
    fi
done
