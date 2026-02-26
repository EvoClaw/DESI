#!/bin/bash
# Monitor DESI chromosome runs and trigger analysis when milestones are reached.
# Usage: bash desi_monitor.sh  (run in background)

RESULTS_DIR="/home/yanlin/popghistory/results"
LOG_DIR="/tmp/desi_logs"
mkdir -p "$RESULTS_DIR"

SMALL_CHROMS="chr19 chr20 chr21 chr22"
MEDIUM_CHROMS="chr13 chr14 chr15 chr16 chr17 chr18"
ALL_CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"

done_small=false
done_medium=false
done_all=false

echo "[$(date)] Monitor started" | tee "$RESULTS_DIR/monitor.log"

while true; do
    sleep 120   # check every 2 minutes

    n_done=$(ls /tmp/desi_pchmm_chr*.npy 2>/dev/null | wc -l)
    n_running=$(ps aux | grep "desi_run_chr" | grep -v grep | wc -l)

    echo "[$(date)] Done: $n_done/22  Running: $n_running" >> "$RESULTS_DIR/monitor.log"

    # Milestone 1: small chromosomes done (chr19-22)
    if [ "$done_small" = false ]; then
        all_small=true
        for c in $SMALL_CHROMS; do
            [ -f "/tmp/desi_pchmm_${c}.npy" ] || { all_small=false; break; }
        done
        if [ "$all_small" = true ]; then
            done_small=true
            echo "[$(date)] Milestone 1: small chroms done — running preliminary analysis" \
                | tee -a "$RESULTS_DIR/monitor.log"
            python3 -u /home/yanlin/popghistory/code/desi_analyze.py \
                --chroms "chr19,chr20,chr21,chr22" \
                --outdir "$RESULTS_DIR/preliminary" \
                >> "$RESULTS_DIR/monitor.log" 2>&1 &
        fi
    fi

    # Milestone 2: medium chromosomes done (add chr13-18)
    if [ "$done_small" = true ] && [ "$done_medium" = false ]; then
        all_med=true
        for c in $MEDIUM_CHROMS; do
            [ -f "/tmp/desi_pchmm_${c}.npy" ] || { all_med=false; break; }
        done
        if [ "$all_med" = true ]; then
            done_medium=true
            echo "[$(date)] Milestone 2: medium chroms done — running mid analysis" \
                | tee -a "$RESULTS_DIR/monitor.log"
            python3 -u /home/yanlin/popghistory/code/desi_analyze.py \
                --chroms "chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22" \
                --outdir "$RESULTS_DIR/midway" \
                >> "$RESULTS_DIR/monitor.log" 2>&1 &
        fi
    fi

    # Milestone 3: all chromosomes done
    if [ "$n_done" -eq 22 ] && [ "$done_all" = false ]; then
        done_all=true
        echo "[$(date)] ALL 22 chromosomes done — running final analysis" \
            | tee -a "$RESULTS_DIR/monitor.log"
        python3 -u /home/yanlin/popghistory/code/desi_analyze.py \
            --outdir "$RESULTS_DIR/final" \
            >> "$RESULTS_DIR/monitor.log" 2>&1
        echo "[$(date)] Final analysis complete" | tee -a "$RESULTS_DIR/monitor.log"
        break
    fi

    # Exit if nothing running and nothing left to wait for
    if [ "$n_running" -eq 0 ] && [ "$n_done" -gt 0 ]; then
        echo "[$(date)] All jobs finished (n_done=$n_done)" \
            | tee -a "$RESULTS_DIR/monitor.log"
        if [ "$done_all" = false ]; then
            python3 -u /home/yanlin/popghistory/code/desi_analyze.py \
                --outdir "$RESULTS_DIR/final" \
                >> "$RESULTS_DIR/monitor.log" 2>&1
        fi
        break
    fi
done
