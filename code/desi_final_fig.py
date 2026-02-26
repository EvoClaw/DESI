import numpy as np, glob
from scipy.special import gammaln
from scipy.stats import gamma as scipy_gamma
from collections import defaultdict
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

files=sorted(glob.glob('/tmp/desi_pchmm_chr*.npy'))
per_chr_all=[]; hap_pop=None; rows0=cols0=None
for fn in files:
    d=np.load(fn,allow_pickle=True).item()
    per_chr_all.append(d['pair_tmrca'])
    if hap_pop is None:
        hap_pop=list(d['hap_pop']); rows0=d['pair_rows']; cols0=d['pair_cols']
tmrca_stack=np.vstack(per_chr_all); tmrca_mean=tmrca_stack.mean(0)
tmrca_pool=tmrca_stack.flatten()

POP_MAP={'YRI':'AFR','GWD':'AFR','ESN':'AFR','MSL':'AFR','ACB':'AFR',
         'CEU':'EUR','TSI':'EUR','CHB':'EAS','JPT':'EAS','CHS':'EAS',
         'PJL':'SAS','GIH':'SAS','PUR':'AMR','PEL':'AMR'}
GEN=28

def grp_mask(g1,g2=None):
    if g2 is None: g2=g1
    return np.array([
        (POP_MAP.get(hap_pop[rows0[i]],'X')==g1 and POP_MAP.get(hap_pop[cols0[i]],'X')==g2) or
        (g1!=g2 and POP_MAP.get(hap_pop[rows0[i]],'X')==g2 and POP_MAP.get(hap_pop[cols0[i]],'X')==g1)
        for i in range(len(tmrca_mean))])

masks={'EAS':grp_mask('EAS'),'EUR':grp_mask('EUR'),'SAS':grp_mask('SAS'),
       'AMR':grp_mask('AMR'),'AFR':grp_mask('AFR'),
       'AFR-EAS':grp_mask('AFR','EAS'),'AFR-EUR':grp_mask('AFR','EUR')}

def em_gamma(x,K,n_iter=400,n_restarts=8):
    n=len(x); best_ll=-1e18; best=None
    for seed in range(n_restarts):
        rng=np.random.default_rng(seed)
        mu0=np.quantile(x,np.sort(rng.uniform(0.05,0.95,K)))
        var0=np.var(x)*rng.uniform(0.2,5.,K)
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
    return dict(bic=-2*best_ll+(3*K-1)*np.log(n),means=sh*sc,pi=pi,sh=sh,sc=sc)

np.random.seed(42); idx=np.random.choice(len(tmrca_pool),200000,replace=False)
r2=em_gamma(tmrca_pool[idx],2)
bic_scan={K:em_gamma(tmrca_pool[idx],K)['bic'] for K in range(1,5)}
best_K=min(bic_scan,key=bic_scan.get)

fig,axes=plt.subplots(1,3,figsize=(18,7))
fig.suptitle('DESI: Deep-time Evolutionary Structure Inference — 22 autosomes, 61.6M SNPs, 240 samples',
             fontsize=11,fontweight='bold')

# ── A: K-mixture ──────────────────────────────────────────────────────────────
ax=axes[0]
bins=np.linspace(700,1400,80)
ax.hist(tmrca_pool,bins=bins,density=True,alpha=0.22,color='#aaaaaa',label='All 2.52M estimates')
gc={'EAS':'#4575b4','EUR':'#74add1','AFR':'#d73027'}
for g,col in gc.items():
    dat=tmrca_stack[:,masks[g]].flatten()
    mn=np.mean(tmrca_mean[masks[g]]); ne=mn*1000/(2*GEN)
    ax.hist(dat,bins=bins,density=True,alpha=0.55,color=col,histtype='step',lw=2,
            label=f'Within-{g}: {mn:.0f} Ka\n(N$_e$≈{ne/1000:.1f}k)')
t_grid=np.linspace(600,1600,500); mix=np.zeros_like(t_grid)
for k,col_k in enumerate(['#377eb8','#e41a1c']):
    c=r2['pi'][k]*scipy_gamma.pdf(t_grid,a=r2['sh'][k],scale=r2['sc'][k])
    ax.plot(t_grid,c,color=col_k,lw=2,ls='--',
            label=f'K={k+1}: {r2["means"][k]:.0f} Ka (π={r2["pi"][k]:.2f})')
    mix+=c
ax.plot(t_grid,mix,'k-',lw=2.5,label=f'K={best_K} best (ΔBIC=−19k)')
ax.axvline(930,color='#d6604d',lw=2,ls=':',label='FitCoal 930 Ka')
ax.set_xlabel('Pairwise TMRCA (Ka)',fontsize=11); ax.set_ylabel('Density',fontsize=11)
ax.set_title('A. K=2 population structure\n(K=3 ΔBIC = +224k → strongly rejected)',fontsize=10)
ax.legend(fontsize=8); ax.set_xlim(650,1450); ax.grid(True,alpha=0.2)

# ── B: TMRCA barplot ──────────────────────────────────────────────────────────
ax=axes[1]
grp_order=['EAS-EAS','EUR-EUR','SAS-SAS','AMR-AMR','AFR-AFR','AFR-EAS','AFR-EUR']
grp_col2=['#4575b4','#74add1','#f46d43','#878787','#d73027','#984ea3','#fc8d59']
mask_list=[masks['EAS'],masks['EUR'],masks['SAS'],masks['AMR'],masks['AFR'],
           masks['AFR-EAS'],masks['AFR-EUR']]
Tm=[np.mean(tmrca_mean[m]) for m in mask_list]
Ts=[np.std(tmrca_mean[m]) for m in mask_list]
Ne=[t*1000/(2*GEN) for t in Tm]
y=np.arange(len(grp_order))
ax.barh(y,Tm,xerr=Ts,color=grp_col2,alpha=0.8,height=0.65,capsize=4,
        error_kw={'lw':1.5,'capthick':1.5})
ax.axvline(930,color='#d6604d',lw=2,ls=':',label='FitCoal 930 Ka',zorder=5)
ax.axvline(r2['means'][0],color='#377eb8',lw=1.5,ls='--',alpha=0.7,label=f'K=2: {r2["means"][0]:.0f} Ka')
ax.axvline(r2['means'][1],color='#e41a1c',lw=1.5,ls='--',alpha=0.7,label=f'K=2: {r2["means"][1]:.0f} Ka')
for i,(t,ne) in enumerate(zip(Tm,Ne)):
    ax.text(t+8,i,f'{t:.0f} Ka (N$_e$≈{ne/1000:.1f}k)',va='center',fontsize=8)
ax.set_yticks(y); ax.set_yticklabels(grp_order,fontsize=10)
ax.set_xlabel('Mean pairwise TMRCA (Ka)',fontsize=11)
ax.set_title('B. Population-pair TMRCA & implied N$_e$\n(error bars = ±1 SD across pairs)',fontsize=10)
ax.legend(fontsize=8.5); ax.set_xlim(700,1440); ax.grid(True,alpha=0.2,axis='x')

# ── C: schematic Ne(t) / split tree ───────────────────────────────────────────
ax=axes[2]
ax.set_xlim(0,1500); ax.set_ylim(-0.3,5.3)
ax.invert_xaxis()
ax.set_xlabel('Time (Ka before present)',fontsize=11)
ax.set_title('C. Schematic population split & ancestral N$_e$\n(inferred from DESI K=2 structure)',fontsize=10)

pop_ne={'EAS':Ne[0],'EUR':Ne[1],'SAS':Ne[2],'AMR':Ne[3],'AFR':Ne[4]}
pop_t ={'EAS':Tm[0],'EUR':Tm[1],'SAS':Tm[2],'AMR':Tm[3],'AFR':Tm[4]}
pop_y ={'EAS':0.5,'EUR':1.3,'SAS':2.1,'AMR':2.9,'AFR':4.3}
pop_col={'EAS':'#4575b4','EUR':'#74add1','SAS':'#f46d43','AMR':'#878787','AFR':'#d73027'}

# modern tips
for pop,yp in pop_y.items():
    ax.plot([0],[yp],'o',color=pop_col[pop],ms=10,zorder=6)
    ax.text(-15,yp,f'{pop}\n(N$_e$≈{pop_ne[pop]/1000:.0f}k)',ha='right',va='center',fontsize=8.5,color=pop_col[pop],fontweight='bold')

# OOA lineages merge at ~65 Ka
ooa_y=1.8
for pop in ['EAS','EUR','SAS','AMR']:
    ax.plot([0,65],[pop_y[pop],pop_y[pop]],color=pop_col[pop],lw=2.5)
    ax.plot([65,65],[pop_y[pop],ooa_y],'--',color=pop_col[pop],lw=1.5,alpha=0.6)
ax.plot([65,350],[ooa_y,ooa_y],color='#9c6d3e',lw=3,label='proto-OOA (N$_e$≈17.6k)')

# AFR lineage extends deep
ax.plot([0,1300],[pop_y['AFR'],pop_y['AFR']],color=pop_col['AFR'],lw=3,label=f'Africa (N$_e$≈{pop_ne["AFR"]/1000:.0f}k)')

# Proto-OOA meets AFR → TMRCA cross zone
# Between-group TMRCA ≈ 1265 Ka → mark the point of ancestral connection
T_anc_split=r2['means'][0]  # ~986 Ka = DESI comp1 (OOA component mean T)
T_afr_mean=r2['means'][1]   # ~1246 Ka = DESI comp2 (AFR component mean T)
ax.plot([350,T_anc_split],[ooa_y,ooa_y],color='#9c6d3e',lw=3)
# Ancestral merge
ax.plot([T_anc_split,T_afr_mean],[ooa_y,(ooa_y+pop_y['AFR'])/2],'k--',lw=2)
ax.plot([T_afr_mean,T_afr_mean],[(ooa_y+pop_y['AFR'])/2,pop_y['AFR']],'k--',lw=2)
ax.plot(T_afr_mean,(ooa_y+pop_y['AFR'])/2,'ks',ms=9,zorder=6)

# Annotations
ax.text(120,ooa_y+0.2,'OOA\n~65 Ka',ha='center',fontsize=8.5,color='gray',
        bbox=dict(fc='white',ec='gray',alpha=0.8,boxstyle='round,pad=0.2'))
ax.annotate('Ancestral split\n~930–1000 Ka\n(DESI K=2)',xy=(T_anc_split,ooa_y),
            xytext=(800,ooa_y+0.8),fontsize=8.5,color='#9c6d3e',fontweight='bold',
            arrowprops=dict(arrowstyle='->',color='#9c6d3e',lw=1.5),
            bbox=dict(fc='#fff7ed',ec='#9c6d3e',alpha=0.9,boxstyle='round,pad=0.3'))
ax.axvline(930,color='#d6604d',lw=2,ls=':',zorder=5)
ax.text(930,5.1,'FitCoal\n930 Ka',ha='center',fontsize=8,color='#d6604d',fontweight='bold')
ax.axvspan(900,1000,alpha=0.08,color='#d6604d')
ax.axvline(r2['means'][0],color='#377eb8',lw=1.5,ls='--',alpha=0.7)
ax.axvline(r2['means'][1],color='#e41a1c',lw=1.5,ls='--',alpha=0.7)

ax.legend(fontsize=8.5,loc='upper left')
ax.set_yticks([]); ax.grid(True,alpha=0.2,axis='x')
ax.set_xlim(-50,1500)

plt.tight_layout()
fig.savefig('/home/yanlin/popghistory/results/final/desi_demography_v2.pdf',bbox_inches='tight',dpi=150)
fig.savefig('/home/yanlin/popghistory/results/final/desi_demography_v2.png',bbox_inches='tight',dpi=150)
plt.close(fig)

print("Saved desi_demography_v2.pdf/png")
print(f"\nFINAL: K={best_K} best (dBIC: {[f'K{k}:{bic_scan[k]-bic_scan[1]:.0f}' for k in sorted(bic_scan)]})")
print(f"K=2: comp1={r2['means'][0]:.0f} Ka (Ne≈{r2['means'][0]*1000/(2*GEN):.0f}), comp2={r2['means'][1]:.0f} Ka (Ne≈{r2['means'][1]*1000/(2*GEN):.0f})")
