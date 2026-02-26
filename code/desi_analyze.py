"""
DESI Phase 4b Analysis — Stage 3 & 4
Loads per-chromosome PCHMM results, aggregates genome-wide,
runs K-mixture inference, and produces all key figures.

Usage:
  python3 desi_analyze.py                  # use all available /tmp/desi_pchmm_*.npy
  python3 desi_analyze.py --chroms chr22   # quick test on one chromosome
"""

import numpy as np
import argparse, os, sys, time, warnings
from collections import defaultdict
from scipy.special import gammaln
from scipy.stats import gamma as scipy_gamma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

warnings.filterwarnings('ignore')

# ── Physical constants ────────────────────────────────────────────────────────
MU       = 1.2e-8   # mutation rate per bp per generation
GEN_TIME = 28       # years per generation

# Key time windows for deep-time analysis (in Ka)
TIME_WINDOWS = {
    'recent':      (0,    200),
    'mid':         (200,  600),
    'bottleneck':  (700,  1200),   # 930 Ka window
    'deep':        (1200, 5000),
}

OUTDIR = "/home/yanlin/popghistory/results"


# ── K-mixture EM ─────────────────────────────────────────────────────────────

def em_gamma_mixture(x, K, n_iter=400, n_restarts=8, seed=42):
    """Fit K-component Gamma mixture to TMRCA values (Ka).
    Returns (pi, shape, scale, log_likelihood, BIC).
    """
    rng = np.random.default_rng(seed)
    n   = len(x)
    best_ll, best_params = -np.inf, None

    for restart in range(n_restarts):
        # Initialise from data quantiles with jitter
        qs    = np.linspace(0.1, 0.9, K)
        mu0   = np.quantile(x, rng.permutation(len(qs))[:K] / len(qs)
                            if K <= len(qs) else qs)
        mu0   = np.quantile(x, np.sort(rng.uniform(0.05, 0.95, K)))
        var0  = np.var(x) * rng.uniform(0.2, 4.0, K)
        shape = np.clip(mu0**2 / (var0 + 1e-8), 0.5, 200.0)
        scale = np.clip(var0 / (mu0 + 1e-8), 0.01, 1000.0)
        pi    = np.ones(K) / K
        ll_prev = -np.inf

        for _ in range(n_iter):
            # E-step
            log_comp = np.zeros((n, K))
            for k in range(K):
                log_comp[:, k] = (
                    (shape[k] - 1) * np.log(x + 1e-300)
                    - x / (scale[k] + 1e-10)
                    - gammaln(shape[k])
                    - shape[k] * np.log(scale[k] + 1e-10)
                    + np.log(pi[k] + 1e-300)
                )
            log_Z   = np.logaddexp.reduce(log_comp, axis=1, keepdims=True)
            resp    = np.exp(log_comp - log_Z)
            ll_curr = float(log_Z.sum())
            if ll_curr - ll_prev < 1e-4:
                break
            ll_prev = ll_curr

            # M-step
            Nk    = resp.sum(axis=0) + 1e-10
            pi    = Nk / Nk.sum()
            mu_k  = (resp * x[:, None]).sum(axis=0) / Nk
            ex2   = (resp * x[:, None]**2).sum(axis=0) / Nk
            var_k = np.clip(ex2 - mu_k**2, 1e-3, None)
            shape = np.clip(mu_k**2 / (var_k + 1e-8), 0.3, 300.0)
            scale = np.clip(var_k / (mu_k + 1e-8), 0.01, 2000.0)

        if ll_curr > best_ll:
            best_ll     = ll_curr
            best_params = (pi.copy(), shape.copy(), scale.copy())

    pi, shape, scale = best_params
    order = np.argsort(shape * scale)
    pi, shape, scale = pi[order], shape[order], scale[order]
    n_params = 3 * K - 1
    bic = -2 * best_ll + n_params * np.log(n)
    return pi, shape, scale, best_ll, bic


def select_K(x, K_max=6, n_restarts=8, seed=42, label=""):
    """BIC-based K selection. Returns {K: result_dict}, best_K."""
    results = {}
    print(f"  K-mixture BIC scan{' (' + label + ')' if label else ''}:", flush=True)
    for K in range(1, K_max + 1):
        pi, shape, scale, ll, bic = em_gamma_mixture(x, K, n_restarts=n_restarts,
                                                       seed=seed)
        means = shape * scale
        results[K] = dict(pi=pi, shape=shape, scale=scale, ll=ll, bic=bic,
                          means=means)
        stars = " ← best" if K == min(range(1, K + 1),
                                       key=lambda k: results[k]['bic']) else ""
        print(f"    K={K}: BIC={bic:,.0f}  ll={ll:,.0f}  "
              f"means={np.round(means).astype(int)}{stars}", flush=True)
    best_K = min(results, key=lambda k: results[k]['bic'])
    return results, best_K


# ── Population utilities ──────────────────────────────────────────────────────

AFRICA_POPS = {'YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ASW', 'ACB', 'MKK', 'MAG'}

def pop_region(pop):
    """Map 1kGP population code to continental region."""
    if pop in AFRICA_POPS:
        return 'AFR'
    elif pop in {'CHB', 'JPT', 'CHS', 'CDX', 'KHV'}:
        return 'EAS'
    elif pop in {'CEU', 'TSI', 'FIN', 'GBR', 'IBS'}:
        return 'EUR'
    elif pop in {'MXL', 'PUR', 'CLM', 'PEL'}:
        return 'AMR'
    elif pop in {'GIH', 'PJL', 'BEB', 'STU', 'ITU'}:
        return 'SAS'
    return 'OTH'


def pair_region_label(pop_r, pop_c):
    r1, r2 = pop_region(pop_r), pop_region(pop_c)
    if r1 == r2:
        return f'within-{r1}'
    return '-'.join(sorted([r1, r2]))


# ── Main analysis ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chroms',   default='all')
    parser.add_argument('--indir',    default='/tmp')
    parser.add_argument('--outdir',   default=OUTDIR)
    parser.add_argument('--K_max',    type=int, default=6)
    parser.add_argument('--n_restarts', type=int, default=8)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    t0 = time.time()

    # ── Load results ──────────────────────────────────────────────────────────
    all_chroms = [f'chr{i}' for i in range(1, 23)]
    if args.chroms == 'all':
        chroms = [c for c in all_chroms
                  if os.path.exists(f'{args.indir}/desi_pchmm_{c}.npy')]
    else:
        chroms = [c.strip() for c in args.chroms.split(',')]

    if not chroms:
        print("No result files found. Run desi_run_chr.py first.", flush=True)
        sys.exit(1)

    print(f"=== DESI Analysis: {len(chroms)} chromosomes ===", flush=True)
    print(f"Chromosomes: {chroms}", flush=True)

    all_tmrca   = []
    hap_pop     = None
    pair_rows_0 = None
    pair_cols_0 = None
    total_snps  = 0; total_wins = 0

    for chrom in chroms:
        fn = f'{args.indir}/desi_pchmm_{chrom}.npy'
        d  = np.load(fn, allow_pickle=True).item()
        all_tmrca.append(d['pair_tmrca'])
        if hap_pop is None:
            hap_pop     = list(d['hap_pop'])
            pair_rows_0 = d['pair_rows']
            pair_cols_0 = d['pair_cols']
            n_pairs     = len(pair_rows_0)
        total_snps += int(d['n_snps'])
        total_wins += int(d['n_windows'])
        print(f"  {chrom}: {d['n_snps']:,} SNPs  "
              f"mean TMRCA={np.mean(d['pair_tmrca']):.0f} Ka", flush=True)

    # Genome-wide mean TMRCA per pair (average across chromosomes)
    tmrca = np.vstack(all_tmrca).mean(axis=0)   # (n_pairs,)
    print(f"\nGenome-wide: {n_pairs:,} pairs  {total_snps:,} SNPs  "
          f"{total_wins:,} windows", flush=True)
    print(f"TMRCA: mean={np.mean(tmrca):.0f} Ka  "
          f"med={np.median(tmrca):.0f} Ka  "
          f"p5={np.percentile(tmrca,5):.0f} Ka  "
          f"p95={np.percentile(tmrca,95):.0f} Ka", flush=True)

    np.save(f'{args.outdir}/desi_genome_tmrca.npy', {
        'pair_tmrca': tmrca, 'pair_rows': pair_rows_0, 'pair_cols': pair_cols_0,
        'hap_pop': np.array(hap_pop), 'chroms': chroms,
        'n_snps': total_snps, 'n_windows': total_wins,
    }, allow_pickle=True)

    # ── Stage 3: Global K-mixture ─────────────────────────────────────────────
    print("\n=== Stage 3: Global K-mixture ===", flush=True)
    km_results, best_K = select_K(tmrca, K_max=args.K_max,
                                   n_restarts=args.n_restarts, label="all pairs")
    print(f"\n  → Best K = {best_K}", flush=True)
    r = km_results[best_K]
    for k in range(best_K):
        print(f"     Component {k+1}: weight={r['pi'][k]:.3f}  "
              f"mean={r['means'][k]:.0f} Ka  shape={r['shape'][k]:.2f}", flush=True)

    np.save(f'{args.outdir}/desi_kmixture_global.npy',
            {'best_K': best_K, 'results': km_results, 'chroms': chroms},
            allow_pickle=True)

    # ── Stage 4a: Population-region K-mixture ─────────────────────────────────
    print("\n=== Stage 4a: Population-region TMRCA ===", flush=True)
    region_tmrca = defaultdict(list)
    for idx in range(n_pairs):
        r_label = pair_region_label(hap_pop[pair_rows_0[idx]],
                                     hap_pop[pair_cols_0[idx]])
        region_tmrca[r_label].append(tmrca[idx])
    region_tmrca = {k: np.array(v) for k, v in region_tmrca.items()}

    region_km = {}
    for label in sorted(region_tmrca.keys()):
        vals = region_tmrca[label]
        if len(vals) < 100:
            continue
        print(f"\n  [{label}] n={len(vals):,}  "
              f"med={np.median(vals):.0f} Ka", flush=True)
        res, bK = select_K(vals, K_max=min(4, args.K_max),
                            n_restarts=6, label=label)
        region_km[label] = {'results': res, 'best_K': bK, 'n': len(vals)}

    np.save(f'{args.outdir}/desi_region_kmixture.npy', region_km, allow_pickle=True)

    # ── Stage 4b: Time-window structure test ──────────────────────────────────
    print("\n=== Stage 4b: Time-window structure test ===", flush=True)
    # For each time window, compute "structure fraction" = fraction of pairs
    # best assigned to a component with mean in the deep-time window
    tw_results = {}
    for win_name, (t_lo, t_hi) in TIME_WINDOWS.items():
        mask    = (tmrca >= t_lo) & (tmrca < t_hi)
        n_win   = mask.sum()
        frac    = n_win / len(tmrca)

        # Within-Africa fraction (expected to be highest for recent windows)
        afr_mask = np.array([pop_region(hap_pop[pair_rows_0[i]]) == 'AFR' and
                              pop_region(hap_pop[pair_cols_0[i]]) == 'AFR'
                              for i in range(n_pairs)])
        afr_frac_in_win = (mask & afr_mask).sum() / max(n_win, 1)

        tw_results[win_name] = dict(t_lo=t_lo, t_hi=t_hi, n=n_win,
                                     frac=frac, afr_frac=afr_frac_in_win)
        print(f"  {win_name:12s} [{t_lo:4d}-{t_hi:4d} Ka]: "
              f"n={n_win:,} ({100*frac:.1f}%)  "
              f"AFR-AFR fraction={100*afr_frac_in_win:.1f}%", flush=True)

    # ── Structure fraction at 930 Ka ──────────────────────────────────────────
    print("\n=== 930 Ka Structure Test ===", flush=True)
    btl_mask = (tmrca >= 700) & (tmrca <= 1200)
    n_btl    = btl_mask.sum()
    print(f"  Pairs with TMRCA in [700-1200 Ka]: {n_btl:,} "
          f"({100*n_btl/len(tmrca):.1f}%)", flush=True)

    # For pairs in this window, compute region breakdown
    print("  Region breakdown of 700-1200 Ka pairs:", flush=True)
    btl_regions = defaultdict(int)
    for idx in np.where(btl_mask)[0]:
        lbl = pair_region_label(hap_pop[pair_rows_0[idx]], hap_pop[pair_cols_0[idx]])
        btl_regions[lbl] += 1
    total_btl = sum(btl_regions.values())
    for lbl, cnt in sorted(btl_regions.items(), key=lambda x: -x[1]):
        print(f"    {lbl}: {cnt:,} ({100*cnt/total_btl:.1f}%)", flush=True)

    # ── Figures ───────────────────────────────────────────────────────────────
    print("\n=== Generating figures ===", flush=True)
    make_figures(tmrca, km_results, best_K, region_tmrca, region_km,
                 tw_results, hap_pop, pair_rows_0, pair_cols_0, args.outdir,
                 chroms)

    print(f"\n=== Analysis complete  ({time.time()-t0:.0f}s) ===", flush=True)
    print(f"Outputs: {args.outdir}/", flush=True)


# ── Figures ───────────────────────────────────────────────────────────────────

COMP_COLORS = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628']
REGION_COLORS = {
    'within-AFR': '#d73027', 'within-EUR': '#4575b4', 'within-EAS': '#91cf60',
    'within-SAS': '#fee090', 'within-AMR': '#fc8d59',
    'AFR-EUR': '#762a83', 'AFR-EAS': '#1b7837', 'AFR-SAS': '#bf812d',
    'AFR-AMR': '#80cdc1', 'EUR-EAS': '#6a3d9a',
}


def make_figures(tmrca, km_results, best_K, region_tmrca, region_km,
                 tw_results, hap_pop, pair_rows_0, pair_cols_0, outdir, chroms):

    r = km_results[best_K]

    # ── Figure 1: Global TMRCA distribution + K-mixture fit ──────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    t_min = max(tmrca.min(), 10)
    bins  = np.logspace(np.log10(t_min), np.log10(tmrca.max()), 100)
    ax.hist(tmrca, bins=bins, density=True, alpha=0.45, color='#4393c3',
            label='All pairs', zorder=2)

    t_grid  = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), 600)
    mix_pdf = np.zeros_like(t_grid)
    for k in range(best_K):
        comp = r['pi'][k] * scipy_gamma.pdf(t_grid, a=r['shape'][k], scale=r['scale'][k])
        ax.plot(t_grid, comp, color=COMP_COLORS[k], lw=1.5, ls='--',
                label=f'Comp {k+1}: {r["means"][k]:.0f} Ka (w={r["pi"][k]:.2f})')
        mix_pdf += comp
    ax.plot(t_grid, mix_pdf, 'k-', lw=2.0, label=f'K={best_K} mixture')

    # Mark 930 Ka
    ax.axvline(930, color='#d6604d', lw=1.2, ls=':', alpha=0.8, label='930 Ka')
    ax.set_xscale('log')
    ax.set_xlabel('Pairwise TMRCA (Ka)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'Genome-wide pairwise TMRCA\n({len(chroms)} chr, '
                 f'{len(tmrca):,} pairs)', fontsize=11)
    ax.legend(fontsize=8, ncol=1)
    ax.grid(True, alpha=0.25, which='both')

    # BIC panel
    ax = axes[1]
    K_vals   = sorted(km_results.keys())
    bic_vals = [km_results[k]['bic'] for k in K_vals]
    ax.plot(K_vals, bic_vals, 'o-', color='#2166ac', lw=2, ms=8)
    ax.axvline(best_K, color='red', ls='--', lw=1.5, label=f'Best K={best_K}')
    ax.set_xticks(K_vals)
    ax.set_xlabel('Number of components K', fontsize=12)
    ax.set_ylabel('BIC (lower = better)', fontsize=12)
    ax.set_title('Model selection (BIC)', fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(f'{outdir}/fig1_tmrca_global.{ext}', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved fig1_tmrca_global", flush=True)

    # ── Figure 2: Population-region TMRCA distributions ──────────────────────
    labels_sorted = sorted(region_tmrca.keys(), key=lambda l: np.median(region_tmrca[l]))
    n_panels      = len(labels_sorted)
    ncols         = min(4, n_panels)
    nrows         = (n_panels + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows),
                              sharey=False)
    axes = np.array(axes).flatten()

    for idx, label in enumerate(labels_sorted):
        ax   = axes[idx]
        vals = region_tmrca[label]
        col  = REGION_COLORS.get(label, '#999999')
        bins = np.logspace(np.log10(max(vals.min(), 10)), np.log10(vals.max()), 60)
        ax.hist(vals, bins=bins, density=True, alpha=0.65, color=col)

        # Overlay K-mixture if fitted
        if label in region_km and region_km[label]['best_K'] > 1:
            bK  = region_km[label]['best_K']
            rr  = region_km[label]['results'][bK]
            tg  = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), 400)
            mp  = np.zeros_like(tg)
            for k in range(bK):
                cp = rr['pi'][k] * scipy_gamma.pdf(tg, a=rr['shape'][k],
                                                    scale=rr['scale'][k])
                mp += cp
            ax.plot(tg, mp, 'k-', lw=1.5)

        ax.axvline(930, color='#d6604d', lw=0.8, ls=':')
        ax.set_xscale('log')
        bK_label = f"K={region_km[label]['best_K']}" if label in region_km else ""
        ax.set_title(f'{label}\nn={len(vals):,}  med={np.median(vals):.0f} Ka'
                     + (f'  {bK_label}' if bK_label else ''), fontsize=8)
        ax.set_xlabel('TMRCA (Ka)', fontsize=8)
        ax.tick_params(labelsize=7)

    for idx in range(len(labels_sorted), len(axes)):
        axes[idx].axis('off')

    plt.suptitle('TMRCA by population pair region', fontsize=11, y=1.01)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(f'{outdir}/fig2_region_tmrca.{ext}', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved fig2_region_tmrca", flush=True)

    # ── Figure 3: Time-window structure — fraction by region ─────────────────
    fig, ax = plt.subplots(figsize=(10, 5))

    t_centers = [(tw['t_lo'] + tw['t_hi']) / 2
                 for tw in [tw_results[w] for w in TIME_WINDOWS]]
    frac_vals = [tw_results[w]['frac'] * 100 for w in TIME_WINDOWS]
    afr_fracs = [tw_results[w]['afr_frac'] * 100 for w in TIME_WINDOWS]
    win_names = list(TIME_WINDOWS.keys())

    x = np.arange(len(win_names))
    w = 0.35
    ax.bar(x - w/2, frac_vals, w, label='All pairs (%)', color='#4393c3', alpha=0.8)
    ax.bar(x + w/2, afr_fracs, w, label='AFR-AFR fraction within window',
           color='#d73027', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{n}\n[{tw_results[n]["t_lo"]}-{tw_results[n]["t_hi"]} Ka]'
                        for n in win_names], fontsize=9)
    ax.set_ylabel('Fraction (%)', fontsize=11)
    ax.set_title('TMRCA time-window analysis: pair fraction and population composition',
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(f'{outdir}/fig3_time_windows.{ext}', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved fig3_time_windows", flush=True)

    # ── Figure 4: Structure fraction vs TMRCA (sliding window) ───────────────
    # Key DESI metric: across the TMRCA range, what is the "AFR vs non-AFR" signal?
    fig, ax = plt.subplots(figsize=(11, 4.5))

    t_bins  = np.logspace(np.log10(20), np.log10(4000), 60)
    t_mids  = (t_bins[:-1] + t_bins[1:]) / 2
    afr_flag  = np.array([pop_region(hap_pop[r]) == 'AFR'
                          for r in pair_rows_0], dtype=bool)
    nafr_flag = np.array([pop_region(hap_pop[r]) != 'AFR'
                          for r in pair_rows_0], dtype=bool)

    frac_afr_afr   = []
    frac_afr_nafr  = []
    frac_nafr_nafr = []

    for i in range(len(t_bins) - 1):
        mask = (tmrca >= t_bins[i]) & (tmrca < t_bins[i+1])
        if mask.sum() < 5:
            frac_afr_afr.append(np.nan)
            frac_afr_nafr.append(np.nan)
            frac_nafr_nafr.append(np.nan)
            continue
        n_total = mask.sum()
        afr_cols_flag = np.array([pop_region(hap_pop[c]) == 'AFR'
                                   for c in pair_cols_0], dtype=bool)
        aa = ((afr_flag & afr_cols_flag) & mask).sum()
        an = (((afr_flag & ~afr_cols_flag) | (~afr_flag & afr_cols_flag)) & mask).sum()
        nn = ((~afr_flag & ~afr_cols_flag) & mask).sum()
        frac_afr_afr.append(100 * aa / n_total)
        frac_afr_nafr.append(100 * an / n_total)
        frac_nafr_nafr.append(100 * nn / n_total)

    fa  = np.array(frac_afr_afr)
    fan = np.array(frac_afr_nafr)
    fnn = np.array(frac_nafr_nafr)

    ax.fill_between(t_mids, 0, fa,  alpha=0.7, color='#d73027', label='AFR–AFR')
    ax.fill_between(t_mids, fa, fa + fan, alpha=0.7, color='#762a83', label='AFR–non-AFR')
    ax.fill_between(t_mids, fa + fan, fa + fan + fnn, alpha=0.7,
                    color='#4575b4', label='non-AFR–non-AFR')

    ax.axvline(930, color='black', lw=1.5, ls='--', label='930 Ka')
    ax.set_xscale('log')
    ax.set_xlim(20, 4000)
    ax.set_ylim(0, 100)
    ax.set_xlabel('Pairwise TMRCA (Ka)', fontsize=12)
    ax.set_ylabel('Fraction of pairs (%)', fontsize=12)
    ax.set_title('Population pair composition across TMRCA range\n'
                 '(DESI reveals deep-time structure)', fontsize=11)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.25, which='both')

    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(f'{outdir}/fig4_structure_vs_tmrca.{ext}', dpi=150,
                    bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved fig4_structure_vs_tmrca", flush=True)


if __name__ == "__main__":
    main()
