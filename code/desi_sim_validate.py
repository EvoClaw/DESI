"""DESI Simulation Validation (E1+E2)

E1: PCHMM calibration vs exact TMRCA from msprime tree sequences.
E2: Does OoA-only model predict the observed AFR-EAS TMRCA difference (376 Ka)?

Models:
  OoA   - Standard OoA: AFR stable Ne~15000, EAS OoA bottleneck at 60Ka, NO deep structure
  A     - FitCoal: panmictic + 930Ka bottleneck
  D     - Cobraa-like: 1.5Ma split, 300Ka admixture (80:20), OoA at 60Ka

Observed: AFR=1251.1Ka, EAS=875.4Ka, diff=375.7Ka
"""
import numpy as np
import msprime
import os
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("/home/yanlin/popghistory/results/validation", exist_ok=True)

GEN    = 28
MU     = 1.2e-8
RHO    = 1e-8
WBP    = 10_000   # 10 kb windows
SL     = 50_000_000   # 50 Mb per replicate (fast; ~sqrt(N) signal sufficient)
NAFR   = 15
NEAS   = 15
SEEDS  = [42, 123, 456]
OA, OE = 1251.1, 875.4
OD     = OA - OE

T_BINS = np.exp(np.linspace(np.log(10), np.log(5000), 40))
TBG    = T_BINS * 1000 / GEN
K      = len(T_BINS)

# Build emission matrix once
_mh = min(int(max(TBG)*2*MU*WBP*2), 2000)
_r  = 2*MU*WBP*TBG
_n  = np.arange(_mh+1)
_lf = np.zeros(_mh+1)
for _i in range(1, _mh+1):
    _lf[_i] = _lf[_i-1] + np.log(_i)
LOG_E = _n[None,:]*np.log(_r[:,None]+1e-300) - _r[:,None] - _lf[None,:]

def build_OoA():
    d = msprime.Demography()
    d.add_population(name="AFR", initial_size=15000)
    d.add_population(name="EAS", initial_size=15000)
    d.add_population(name="ANC", initial_size=15000)
    d.add_population_parameters_change(time=10000/GEN, population="EAS", initial_size=8000)
    d.add_population_parameters_change(time=50000/GEN, population="EAS", initial_size=200)
    d.add_population_split(time=60000/GEN, derived=["AFR","EAS"], ancestral="ANC")
    return d

def build_A():
    d = msprime.Demography()
    d.add_population(name="ALL", initial_size=10000)
    d.add_population_parameters_change(time=813000/GEN, population="ALL", initial_size=1280)
    d.add_population_parameters_change(time=930000/GEN, population="ALL", initial_size=10000)
    return d

def build_D():
    d = msprime.Demography()
    d.add_population(name="AFR",  initial_size=15000)
    d.add_population(name="EAS",  initial_size=15000)
    d.add_population(name="popA", initial_size=12000)
    d.add_population(name="popB", initial_size=5000)
    d.add_population(name="ANC",  initial_size=15000)
    d.add_population_parameters_change(time=10000/GEN, population="EAS", initial_size=8000)
    d.add_population_parameters_change(time=50000/GEN, population="EAS", initial_size=200)
    d.add_population_split(time=60000/GEN, derived=["AFR","EAS"], ancestral="popA")
    d.add_mass_migration(time=300000/GEN, source="popB", dest="popA", proportion=1.0)
    d.add_population_split(time=1500000/GEN, derived=["popA","popB"], ancestral="ANC")
    return d

def exact_tmrca_ka(ts):
    div = ts.divergence_matrix(span_normalise=True)
    return div / (2*MU) * GEN / 1000

def pchmm_all_pairs(ts, pops):
    n  = ts.num_samples
    G  = ts.genotype_matrix().T
    p  = ts.tables.sites.position.astype(np.int64)
    nw = int(ts.sequence_length) // WBP
    wi = p // WBP
    ok = wi < nw
    Gv = G[:, ok];  wiv = wi[ok]
    out = []
    for i in range(n):
        for j in range(i+1, n):
            hs = Gv[i] != Gv[j]
            hw = np.clip(np.bincount(wiv[hs], minlength=nw), 0, _mh)
            tlp = LOG_E[:, hw].sum(axis=1)
            tlp -= tlp.max()
            post = np.exp(tlp); post /= post.sum()
            out.append((pops[i], pops[j], float(post @ T_BINS)))
    return out

def grp_means(pairs):
    aa = [t for pi,pj,t in pairs if pi=="AFR" and pj=="AFR"]
    ee = [t for pi,pj,t in pairs if pi=="EAS" and pj=="EAS"]
    return (np.mean(aa) if aa else np.nan), (np.mean(ee) if ee else np.nan)

print("="*60)
print("DESI Simulation Validation (E1+E2)")
print(f"SL={SL//1e6:.0f}Mb  N_AFR={NAFR}  N_EAS={NEAS}  reps={len(SEEDS)}")
print(f"Observed: AFR={OA:.1f}  EAS={OE:.1f}  diff={OD:.1f} Ka")
print("="*60)

ALL = {}
COLORS = {"OoA":"#2196F3","Model_A":"#F44336","Model_D":"#4CAF50"}
SPECS = [
    ("OoA",     build_OoA, {"AFR":NAFR,"EAS":NEAS}, False),
    ("Model_A", build_A,   {"ALL":NAFR+NEAS},        True),
    ("Model_D", build_D,   {"AFR":NAFR,"EAS":NEAS},  False),
]

for mname, mfn, samp, is_pan in SPECS:
    print(f"\n--- {mname} ---")
    reps = []
    for seed in SEEDS:
        t0 = time.time()
        ts = msprime.sim_ancestry(samples=samp, demography=mfn(),
                                  sequence_length=SL, recombination_rate=RHO,
                                  random_seed=seed)
        ts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)
        # Population labels per haplotype
        if is_pan:
            pops = ["AFR"]*(2*NAFR) + ["EAS"]*(2*NEAS)
        else:
            pm = {}
            for ind in ts.individuals():
                pn = ts.population(ind.population).metadata.get("name",
                     str(ind.population))
                for nd in ind.nodes: pm[nd] = pn
            pops = [pm[i] for i in range(ts.num_samples)]
        # Exact
        tm = exact_tmrca_ka(ts)
        afr_h = [i for i,p in enumerate(pops) if p=="AFR"]
        eas_h = [i for i,p in enumerate(pops) if p=="EAS"]
        ea = np.mean([tm[i,j] for ii,i in enumerate(afr_h) for j in afr_h[ii+1:]])
        ee = np.mean([tm[i,j] for ii,i in enumerate(eas_h) for j in eas_h[ii+1:]])
        # PCHMM
        pairs = pchmm_all_pairs(ts, pops)
        pa, pe = grp_means(pairs)
        reps.append({"ea":ea,"ee":ee,"ed":ea-ee,"pa":pa,"pe":pe,"pd":pa-pe})
        print(f"  s={seed}: exact AFR={ea:.0f} EAS={ee:.0f} diff={ea-ee:.0f}  "
              f"pchmm AFR={pa:.0f} EAS={pe:.0f} diff={pa-pe:.0f}  ({time.time()-t0:.0f}s)")
    ALL[mname] = reps
    eds = [r["ed"] for r in reps]; pds = [r["pd"] for r in reps]
    print(f"  >>> exact diff={np.mean(eds):.1f}±{np.std(eds):.1f}  "
          f"pchmm diff={np.mean(pds):.1f}±{np.std(pds):.1f}  observed={OD:.1f} Ka")

np.save("/home/yanlin/popghistory/results/validation/sim_validate_results.npy",
        {"models":ALL,"obs_afr":OA,"obs_eas":OE,"obs_diff":OD})

# Figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle("DESI Simulation Validation (E1+E2)", fontsize=12)

# (a) AFR-EAS difference
ax = axes[0]
for xi, mname in enumerate(["OoA","Model_A","Model_D"]):
    reps = ALL[mname]
    eds = [r["ed"] for r in reps]; pds = [r["pd"] for r in reps]
    ax.bar(xi-0.18, np.mean(eds), 0.32, yerr=np.std(eds),
           color=COLORS[mname], alpha=0.9, capsize=4, label=f"{mname} exact")
    ax.bar(xi+0.18, np.mean(pds), 0.32, yerr=np.std(pds),
           color=COLORS[mname], alpha=0.45, hatch="//", capsize=4, label=f"{mname} PCHMM")
ax.axhline(OD, color="black", ls="--", lw=2, label=f"Observed ({OD:.0f} Ka)")
ax.axhline(0, color="grey", lw=0.5)
ax.set_xticks([0,1,2]); ax.set_xticklabels(["OoA","Model A","Model D"])
ax.set_ylabel("AFR − EAS TMRCA (Ka)")
ax.set_title("(a) Predicted vs observed AFR-EAS diff\nsolid=exact, hatched=PCHMM")
ax.legend(fontsize=7)

# (b) calibration
ax = axes[1]
for mname, reps in ALL.items():
    for r in reps:
        ax.scatter(r["ea"], r["pa"], c=COLORS[mname], s=60, marker="o", alpha=0.8)
        ax.scatter(r["ee"], r["pe"], c=COLORS[mname], s=60, marker="s", alpha=0.6)
all_v = [r[k] for rps in ALL.values() for r in rps for k in ("ea","ee","pa","pe")]
lim = [min(all_v)*0.85, max(all_v)*1.15]
ax.plot(lim, lim, "k--", lw=1.5, label="y=x")
ax.set_xlabel("Exact TMRCA (Ka)"); ax.set_ylabel("PCHMM TMRCA (Ka)")
ax.set_title("(b) PCHMM calibration (E1)\n(circle=AFR, square=EAS)")
for m in COLORS: ax.scatter([],[],c=COLORS[m],s=40,label=m)
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("/home/yanlin/popghistory/results/validation/fig_sim_validate.png",
            dpi=150, bbox_inches="tight")

print("\n"+"="*60)
print(f"{'Model':<10} {'ExAFR':>8} {'ExEAS':>8} {'ExDiff':>9} {'PhAFR':>9} {'PhEAS':>9} {'PhDiff':>9}")
for mname, reps in ALL.items():
    ea=np.mean([r["ea"]for r in reps]); ee=np.mean([r["ee"]for r in reps])
    ed=np.mean([r["ed"]for r in reps]); pa=np.mean([r["pa"]for r in reps])
    pe=np.mean([r["pe"]for r in reps]); pd=np.mean([r["pd"]for r in reps])
    print(f"  {mname:<10} {ea:>8.1f} {ee:>8.1f} {ed:>9.1f} {pa:>9.1f} {pe:>9.1f} {pd:>9.1f}")
print(f"  {'Observed':<10} {OA:>8.1f} {OE:>8.1f} {OD:>9.1f}")
print("\nSaved: results/validation/sim_validate_results.npy")
print("Saved: results/validation/fig_sim_validate.png")
