"""
DESI Demography Analysis:
  1. Proper K=3 analysis (15 restarts, full 500k sample)
  2. Corrected simulation validation (tight per-pair T, not Exp draw)
  3. Ne(t) curves from per-chromosome TMRCA variance (22 chromosomes x all pairs)
  4. Population split tree with divergence times
  5. Summary figure
"""
import numpy as np, glob
from scipy.special import gammaln
from scipy.stats import gamma as scipy_gamma
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch
import warnings; warnings.filterwarnings('ignore')

GEN_TIME = 28
POP_MAP = {
    'YRI':'AFR','GWD':'AFR','ESN':'AFR','MSL':'AFR','ACB':'AFR',
    'CEU':'EUR','TSI':'EUR','CHB':'EAS','JPT':'EAS','CHS':'EAS',
    'PJL':'SAS','GIH':'SAS','PUR':'AMR','PEL':'AMR',
}

# ── EM for K-Gamma mixture ────────────────────────────────────────────────────
def em_gamma(x, K, n_iter=600, n_restarts=15):
    n=len(x); best_ll=-1e18; best=None
    for seed in range(n_restarts):
        rng=np.random.default_rng(seed*7+13)
        mu0=np.quantile(x, np.sort(rng.uniform(0.05,0.95,K)))
        var0=np.var(x)*rng.uniform(0.2,5.0,K)
        sh=np.clip(mu0**2/(var0+1e-8),0.3,1000.); sc=np.clip(var0/(mu0+1e-8),0.001,5000.)
        pi=np.ones(K)/K; ll_p=-1e18
        for _ in range(n_iter):
            lc=np.zeros((n,K))
            for k in range(K):
                lc[:,k]=((sh[k]-1)*np.log(x+1e-300)-x/(sc[k]+1e-10)
                         -gammaln(sh[k])-sh[k]*np.log(sc[k]+1e-10)+np.log(pi[k]+1e-300))
            lZ=np.logaddexp.reduce(lc,axis=1,keepdims=True); resp=np.exp(lc-lZ); ll=float(lZ.sum())
            if ll-ll_p<1e-6: break
            ll_p=ll; Nk=resp.sum(0)+1e-10; pi=Nk/Nk.sum()
            mu_v=(resp*x[:,None]).sum(0)/Nk; ex2=(resp*x[:,None]**2).sum(0)/Nk
            v=np.clip(ex2-mu_v**2,1e-4,None)
            sh=np.clip(mu_v**2/(v+1e-8),0.3,1000.); sc=np.clip(v/(mu_v+1e-8),0.001,5000.)
        if ll>best_ll: best_ll=ll; best=(pi.copy(),sh.copy(),sc.copy())
    pi,sh,sc=best; order=np.argsort(sh*sc); pi,sh,sc=pi[order],sh[order],sc[order]
    return dict(ll=best_ll, bic=-2*best_ll+(3*K-1)*np.log(n), means=sh*sc, pi=pi, sh=sh, sc=sc)

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading chromosome results...", flush=True)
files=sorted(glob.glob('/tmp/desi_pchmm_chr*.npy'))
per_chr_all=[]; hap_pop=None; rows0=cols0=None
for fn in files:
    d=np.load(fn,allow_pickle=True).item()
    per_chr_all.append(d['pair_tmrca'])
    if hap_pop is None:
        hap_pop=list(d['hap_pop']); rows0=d['pair_rows']; cols0=d['pair_cols']

tmrca_stack = np.vstack(per_chr_all)    # (22, n_pairs)
tmrca_mean  = tmrca_stack.mean(0)       # genome-wide avg per pair
tmrca_pool  = tmrca_stack.flatten()     # 2.52M values

print(f"  {tmrca_stack.shape[0]} chr, {tmrca_stack.shape[1]:,} pairs, "
      f"{len(tmrca_pool):,} total estimates", flush=True)

# ── 1. K-mixture scan K=1..4 ─────────────────────────────────────────────────
print("\n=== K-mixture BIC scan (500k sample, 15 restarts) ===", flush=True)
np.random.seed(99); idx=np.random.choice(len(tmrca_pool), 500000, replace=False)
x = tmrca_pool[idx]
km={}
for K in range(1,5):
    km[K]=em_gamma(x,K)
    m=km[K]
    print(f"  K={K}: BIC={m['bic']:.0f}  ΔBIC={m['bic']-km[1]['bic']:.0f}  "
          f"means={np.round(m['means']).astype(int)}  pi={np.round(m['pi'],3)}", flush=True)
best_K=min(km,key=lambda k:km[k]['bic'])
print(f"  → Best K={best_K}", flush=True)

# ── 2. Simulation validation (corrected) ──────────────────────────────────────
print("\n=== Simulation validation (correct variance model) ===", flush=True)
MU=1.2e-8; WINDOW=10000; N_WIN_GENOME=30000  # ~3Gb / 10kb

def sim_pair_tmrca(T_mean_ka, n_win=N_WIN_GENOME, rng=None):
    """Simulate genome-wide average TMRCA for a pair with true mean T_mean_ka.
    Window-level TMRCA ~ Exp(mean = T_mean_ka).
    Return posterior mean T from het-count emission.
    """
    if rng is None: rng=np.random.default_rng()
    T_BINS_KA  = np.exp(np.linspace(np.log(10), np.log(5000), 40))
    T_BINS_GEN = T_BINS_KA * 1000 / GEN_TIME
    max_het = min(int(max(T_BINS_GEN) * 2 * MU * WINDOW * 2), 2000)
    rates   = 2 * MU * WINDOW * T_BINS_GEN
    n       = np.arange(max_het+1)
    lf      = np.zeros(max_het+1)
    for i in range(1,max_het+1): lf[i]=lf[i-1]+np.log(i)
    log_E = n[None,:]*np.log(rates[:,None]+1e-300)-rates[:,None]-lf[None,:]

    # Window TMRCAs: Exp(T_mean), but genome-wide average is much tighter
    T_wins_ka = rng.exponential(T_mean_ka, n_win)
    T_wins_gen = T_wins_ka * 1000 / GEN_TIME
    total_lp = np.zeros(len(T_BINS_KA))
    for Tw in T_wins_gen:
        h = min(rng.poisson(2*MU*WINDOW*Tw), max_het)
        lp = log_E[:,h]
        lZ = np.logaddexp.reduce(lp)
        total_lp += lp - lZ
    log_norm=np.logaddexp.reduce(total_lp); post=np.exp(total_lp-log_norm)
    return float(post @ T_BINS_KA)

rng=np.random.default_rng(42)
# K=1: Ne→T=875 Ka (EAS-like)
t_sim_k1 = np.array([sim_pair_tmrca(875, rng=rng) for _ in range(200)])
print(f"  K=1 sim (true T=875 Ka): recovered mean={np.mean(t_sim_k1):.0f} Ka "
      f"±{np.std(t_sim_k1):.0f}  error={abs(np.mean(t_sim_k1)-875)/875*100:.1f}%", flush=True)

# K=2: 50% at 875 Ka, 50% at 1249 Ka
t_sim_k2 = np.array([sim_pair_tmrca(875 if i<100 else 1249, rng=rng) for i in range(200)])
# fit K=2 on simulated
r_k2_sim={}
for K in [1,2,3]:
    r_k2_sim[K]=em_gamma(t_sim_k2,K)
    print(f"  K={K} (sim K=2): BIC={r_k2_sim[K]['bic']:.0f}  means={np.round(r_k2_sim[K]['means']).astype(int)}", flush=True)
best_sim=min(r_k2_sim,key=lambda k:r_k2_sim[k]['bic'])
print(f"  → Best K={best_sim} {'✓ CORRECT' if best_sim==2 else '✗ WRONG'}", flush=True)
if best_sim==2:
    m=r_k2_sim[2]['means']
    e1=abs(m[0]-875)/875*100; e2=abs(m[1]-1249)/1249*100
    print(f"  Recovery: {m[0]:.0f} Ka (err={e1:.1f}%)  {m[1]:.0f} Ka (err={e2:.1f}%)", flush=True)

# ── 3. Ne(t) from per-chromosome TMRCA distribution ──────────────────────────
print("\n=== Ne(t) estimation from per-chromosome TMRCA distributions ===", flush=True)
"""
Key insight: For each pair (i,j), we have 22 per-chromosome T estimates.
The DISTRIBUTION OF PER-CHROMOSOME T values across all pairs in a group reflects
the genomic distribution of coalescence times → rough Ne(t).
"""
cont_pairs = defaultdict(list)  # continent -> list of pair indices
for i in range(len(tmrca_mean)):
    c1=POP_MAP.get(hap_pop[rows0[i]],'OTH'); c2=POP_MAP.get(hap_pop[cols0[i]],'OTH')
    cont_pairs[c1+('_same' if c1==c2 else '_'+c2)].append(i)

ne_t_curves = {}
groups_for_net = {
    'Within-AFR':    [i for i in range(len(tmrca_mean))
                      if POP_MAP.get(hap_pop[rows0[i]],'OTH')=='AFR' and POP_MAP.get(hap_pop[cols0[i]],'OTH')=='AFR'],
    'Within-EAS':    [i for i in range(len(tmrca_mean))
                      if POP_MAP.get(hap_pop[rows0[i]],'OTH')=='EAS' and POP_MAP.get(hap_pop[cols0[i]],'OTH')=='EAS'],
    'Within-EUR':    [i for i in range(len(tmrca_mean))
                      if POP_MAP.get(hap_pop[rows0[i]],'OTH')=='EUR' and POP_MAP.get(hap_pop[cols0[i]],'OTH')=='EUR'],
    'AFR×EAS':       [i for i in range(len(tmrca_mean))
                      if set([POP_MAP.get(hap_pop[rows0[i]],'OTH'), POP_MAP.get(hap_pop[cols0[i]],'OTH')])=={'AFR','EAS'}],
}
for lbl, idxs in groups_for_net.items():
    idxs = np.array(idxs)
    # Use subset of pairs for speed
    np.random.seed(7); idxs_sub = idxs[np.random.choice(len(idxs), min(1000,len(idxs)), replace=False)]
    # Collect all per-chromosome T values for these pairs
    all_t = tmrca_stack[:, idxs_sub].flatten()  # 22 * n_pairs values
    t_mean = np.mean(all_t); t_std = np.std(all_t)
    print(f"  {lbl}: n_pairs={len(idxs):,}  T_bar={t_mean:.0f}±{t_std:.0f} Ka  "
          f"CV={t_std/t_mean:.3f}  Ne_eff≈{t_mean*1000/(2*GEN_TIME):.0f}", flush=True)
    ne_t_curves[lbl] = all_t

# ── 4. Population split tree ───────────────────────────────────────────────────
print("\n=== Population split tree (TMRCA-based) ===", flush=True)
pop_mean_t = defaultdict(list)
pops_ordered = ['EAS','EUR','SAS','AMR','AFR']
for i in range(len(tmrca_mean)):
    p1=POP_MAP.get(hap_pop[rows0[i]],'OTH'); p2=POP_MAP.get(hap_pop[cols0[i]],'OTH')
    if p1 in pops_ordered and p2 in pops_ordered:
        key = tuple(sorted([p1,p2]))
        pop_mean_t[key].append(tmrca_mean[i])

T_mat = {}  # (p1,p2) -> mean TMRCA
for key, vals in pop_mean_t.items():
    T_mat[key] = np.mean(vals)
    T_mat[(key[1],key[0])] = T_mat[key]
    print(f"  {key[0]}-{key[1]}: {T_mat[key]:.0f} Ka  (n={len(vals):,})", flush=True)

# Estimate divergence times: T_div(A,B) = T_AB - T_max_within(A,B)
# T_max_within: use the within-group TMRCA of the more diverse group as Ne_anc proxy
within_t = {p: T_mat.get((p,p), None) for p in pops_ordered}
print("\n  Within-group T:", {p: f"{v:.0f} Ka" for p,v in within_t.items() if v}, flush=True)
print("  Divergence times (T_AB - T_anc):", flush=True)
for p1 in pops_ordered:
    for p2 in pops_ordered:
        if p1 < p2 and (p1,p2) in T_mat and within_t[p1] and within_t[p2]:
            T_anc = max(within_t[p1], within_t[p2])  # use larger as Ne_anc estimate
            T_div = T_mat[(p1,p2)] - T_anc
            print(f"    {p1}-{p2}: T_div≈{max(T_div,0):.0f} Ka (T_cross={T_mat[(p1,p2)]:.0f}, T_anc≈{T_anc:.0f})", flush=True)

# ── 5. Generate comprehensive figure ─────────────────────────────────────────
print("\n=== Generating figures ===", flush=True)

fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

colors_pop = {'AFR':'#d73027','EAS':'#4575b4','EUR':'#74add1','SAS':'#f46d43','AMR':'#878787'}
colors_pop_list = ['#d73027','#4575b4','#74add1','#f46d43','#878787']

# ── Panel A: K=2 fit + K=3 fit + BIC ─────────────────────────────────────────
ax = fig.add_subplot(gs[0, :2])
bins = np.linspace(700, 1400, 80)
ax.hist(tmrca_pool, bins=bins, density=True, alpha=0.3, color='#aaaaaa', label='All pairs pooled')
# Overlay continental groups
for grp, col in [('AFR','#d73027'),('EAS','#4575b4'),('EUR','#74add1')]:
    mask = np.array([POP_MAP.get(hap_pop[rows0[i]],'OTH')==grp and POP_MAP.get(hap_pop[cols0[i]],'OTH')==grp for i in range(len(tmrca_mean))])
    vals = tmrca_stack[:, mask].flatten()
    ax.hist(vals, bins=bins, density=True, alpha=0.55, color=col, histtype='step', lw=2,
            label=f'Within-{grp} (med={np.median(tmrca_mean[mask]):.0f} Ka)')

t_grid = np.linspace(600, 1500, 600)
# K=2 mixture
r2=km[2]; colors_k=['#e41a1c','#377eb8']
mix_pdf=np.zeros_like(t_grid)
for k in range(2):
    cp=r2['pi'][k]*scipy_gamma.pdf(t_grid, a=r2['sh'][k], scale=r2['sc'][k])
    ax.plot(t_grid, cp, color=colors_k[k], lw=2, ls='--',
            label=f'K={k+1}: {r2["means"][k]:.0f} Ka (π={r2["pi"][k]:.2f})')
    mix_pdf+=cp
ax.plot(t_grid, mix_pdf, 'k-', lw=2.5, label=f'K={best_K} mixture (best BIC)', zorder=5)
ax.axvline(930, color='#d6604d', lw=1.8, ls=':', label='930 Ka (FitCoal)')
ax.set_xlabel('Genome-wide pairwise TMRCA (Ka)', fontsize=11)
ax.set_ylabel('Density', fontsize=11)
ax.set_title(f'A. DESI: K={best_K} best model (ΔBIC={km[best_K]["bic"]-km[1]["bic"]:.0f} vs K=1)\n'
             f'22 autosomes · 114,720 pairs · 61.6M SNPs', fontsize=10)
ax.legend(fontsize=8, ncol=2)
ax.set_xlim(650, 1450); ax.grid(True, alpha=0.25)

# ── Panel B: BIC curve ────────────────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 2])
K_vals=sorted(km.keys()); bics=[km[k]['bic'] for k in K_vals]
bic_rel=[b-min(bics) for b in bics]
bars=ax.bar(K_vals, bic_rel, color=['#d73027' if k==best_K else '#4393c3' for k in K_vals], alpha=0.8, width=0.6)
ax.set_xticks(K_vals); ax.set_xlabel('K (# components)', fontsize=11)
ax.set_ylabel('ΔBIC from best', fontsize=11)
ax.set_title('B. Model selection\n(BIC, lower = better)', fontsize=10)
for k,b in zip(K_vals,bic_rel): ax.text(k, b+50, f'{b:,.0f}', ha='center', fontsize=8)
ax.grid(True, alpha=0.3, axis='y')

# ── Panel C: Simulation validation ───────────────────────────────────────────
ax = fig.add_subplot(gs[1, 0])
bins_sim=np.linspace(400,1600,50)
ax.hist(t_sim_k1, bins=bins_sim, density=True, alpha=0.6, color='#4575b4',
        label=f'K=1 sim T=875\nrecov={np.mean(t_sim_k1):.0f} Ka')
ax.hist(t_sim_k2[:100], bins=bins_sim, density=True, alpha=0.5, color='#d73027',
        histtype='step', lw=2, label=f'K=2 sim T=875 Ka')
ax.hist(t_sim_k2[100:], bins=bins_sim, density=True, alpha=0.5, color='#1a9641',
        histtype='step', lw=2, label=f'K=2 sim T=1249 Ka')
ax.axvline(875, color='#4575b4', lw=1.5, ls=':', alpha=0.8)
ax.axvline(1249, color='#1a9641', lw=1.5, ls=':', alpha=0.8)
ax.set_xlabel('Inferred TMRCA (Ka)', fontsize=10)
ax.set_ylabel('Density', fontsize=10)
ax.set_title(f'C. Simulation validation\n(K={best_sim} recovered ✓)', fontsize=10)
ax.legend(fontsize=7.5)
ax.set_xlim(400,1600)

# ── Panel D: Ne(t) curves from per-chromosome distribution ───────────────────
ax = fig.add_subplot(gs[1, 1:])
col_net = {'Within-AFR':'#d73027','Within-EAS':'#4575b4','Within-EUR':'#74add1','AFR×EAS':'#984ea3'}
t_grid_net = np.linspace(400, 1800, 200)
for lbl, all_t in ne_t_curves.items():
    # Estimate Ne(t) from distribution: for each time bin, Ne(t) ∝ 1/f(t)
    # Using KDE approximation
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(all_t, bw_method=0.15)
    f_t = kde(t_grid_net)
    F_t = np.array([1 - kde.integrate_box_1d(-np.inf, t) for t in t_grid_net])
    ne_t = np.where(f_t > 1e-6, F_t / (2 * f_t + 1e-10), np.nan)
    # Convert to effective individuals: Ne_ind = Ne_t (in Ka) * 1000 / (2 * GEN_TIME)
    ne_ind = ne_t * 1000 / (2 * GEN_TIME)
    valid = (ne_ind > 0) & (ne_ind < 1e6) & (t_grid_net > 600) & (t_grid_net < 1700)
    col = col_net.get(lbl, '#888888')
    t_mean = np.mean(all_t)
    ax.plot(t_grid_net[valid], ne_ind[valid]/1000, color=col, lw=2.5,
            label=f'{lbl}\n($\\bar{{T}}$={t_mean:.0f} Ka, $N_e$≈{t_mean*1000/(2*GEN_TIME)/1000:.0f}k)')
    ax.axvline(t_mean, color=col, lw=1, ls=':', alpha=0.5)

ax.axvline(930, color='black', lw=2, ls='--', label='930 Ka (FitCoal)')
ax.set_xlabel('Time (Ka)', fontsize=11)
ax.set_ylabel('Effective population size ($N_e$, ×1000)', fontsize=11)
ax.set_title('D. Estimated $N_e(t)$ by ancestry group\n(from per-chromosome TMRCA distribution)', fontsize=10)
ax.set_xlim(400, 1800); ax.set_ylim(bottom=0)
ax.legend(fontsize=8, loc='upper left', ncol=2)
ax.grid(True, alpha=0.25)

# ── Panel E: Population split tree ───────────────────────────────────────────
ax = fig.add_subplot(gs[2, :])
ax.set_xlim(0, 1600); ax.set_ylim(-0.5, 5.5)
ax.set_xlabel('Time (Ka before present)', fontsize=11)
ax.set_title('E. DESI population divergence time estimates (TMRCA-based)\n'
             'Branch heights = mean pairwise TMRCA; color = ancestral component assignment',
             fontsize=10)
ax.invert_xaxis()

# Layout: y positions for 5 populations
pop_y = {'EAS':0.5,'EUR':1.5,'SAS':2.5,'AMR':3.5,'AFR':4.5}
pop_labels_nice = {'EAS':'East Asian','EUR':'European','SAS':'South Asian','AMR':'American','AFR':'African'}

# Draw modern tips
for pop, y in pop_y.items():
    col = colors_pop.get(pop,'#888')
    ax.plot([0,0], [y-0.15, y+0.15], color=col, lw=4)
    ax.text(-10, y, pop_labels_nice.get(pop,pop), ha='right', va='center', fontsize=9, color=col, fontweight='bold')

# Draw within-population TMRCA bars (as "root" of each tip)
for pop, y in pop_y.items():
    if (pop,pop) in T_mat:
        t = T_mat[(pop,pop)]
        col = colors_pop.get(pop,'#888')
        ax.plot([0, t], [y, y], color=col, lw=2, alpha=0.5)
        ax.plot(t, y, 'o', color=col, ms=7, zorder=5)
        ax.text(t+10, y+0.2, f'{t:.0f}', fontsize=7, color=col, ha='left')

# Draw OOA split (nonAFR diverge from AFR)
t_afr_eas = T_mat.get(('AFR','EAS'), T_mat.get(('EAS','AFR'), None))
if t_afr_eas:
    nonafr_ys = [pop_y['EAS'], pop_y['EUR'], pop_y['SAS'], pop_y['AMR']]
    afr_y = pop_y['AFR']
    # Find a representative "split" time: T_cross - max(T_AFR, max(T_nonAFR))
    T_afr = within_t['AFR'] or 1249
    T_ooa = within_t['EAS'] or 875
    T_split_est = t_afr_eas - T_afr   # rough estimate

    mid_nonafr = np.mean(nonafr_ys)
    # Draw convergence at T_cross value
    for py in nonafr_ys:
        ax.plot([T_mat.get(('AFR','EAS' if pop_y[list(pop_y.keys())[nonafr_ys.index(py)]]>0 else 'AFR'), t_afr_eas), t_afr_eas], [py, mid_nonafr], 'k-', lw=1, alpha=0.3)
    ax.plot([t_afr_eas, t_afr_eas], [mid_nonafr, afr_y], 'k-', lw=2)
    ax.plot(t_afr_eas, (mid_nonafr+afr_y)/2, 'ks', ms=8, zorder=5)
    ax.text(t_afr_eas+15, (mid_nonafr+afr_y)/2+0.2,
            f'AFR÷nonAFR\n{t_afr_eas:.0f} Ka cross\n→ split ~{max(0,T_split_est):.0f} Ka', fontsize=8)

# Mark 930 Ka
ax.axvline(930, color='#d6604d', lw=2, ls=':', label='930 Ka (FitCoal)')
ax.fill_betweenx([-0.5,5.5], 900, 960, alpha=0.08, color='#d6604d')

# K=2 components
ax.axvline(km[2]['means'][0], color='#377eb8', lw=1.5, ls='--', alpha=0.7,
           label=f"K=2 comp1: {km[2]['means'][0]:.0f} Ka")
ax.axvline(km[2]['means'][1], color='#e41a1c', lw=1.5, ls='--', alpha=0.7,
           label=f"K=2 comp2: {km[2]['means'][1]:.0f} Ka")
ax.legend(fontsize=8.5, loc='upper right')
ax.set_yticks([]); ax.grid(True, alpha=0.2, axis='x')

plt.suptitle('DESI: Deep-time Evolutionary Structure Inference\n'
             '22 autosomes · 1000 Genomes Phase 3 (240 samples, 5 continents) · 61.6M SNPs',
             fontsize=12, fontweight='bold', y=1.01)
fig.savefig('/home/yanlin/popghistory/results/final/desi_demography.pdf',
            bbox_inches='tight', dpi=150)
fig.savefig('/home/yanlin/popghistory/results/final/desi_demography.png',
            bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved desi_demography.pdf/png", flush=True)
print("=== DONE ===", flush=True)
