"""Fix Model D msprime sampling API and test DESI power for structured model."""
import numpy as np
import msprime
from scipy.special import gammaln
import time

GEN_TIME = 28
MU = 1.2e-8
WINDOW_BP = 1_000_000
SEQ_LEN = 50_000_000
N_SAMPLES = 30
SEED = 42
np.random.seed(SEED)

def fit_gamma_em(data, K, n_iter=200, tol=1e-7):
    n = len(data)
    if n < K*20: return None
    means = np.quantile(data, np.linspace(0.1,0.9,K))
    alphas = np.full(K, 2.0); betas = alphas/means; w = np.ones(K)/K
    log_d = np.log(data+1e-10); prev=-np.inf
    for _ in range(n_iter):
        lr = np.zeros((n,K))
        for k in range(K):
            lr[:,k] = np.log(w[k]+1e-300)+(alphas[k]-1)*log_d-betas[k]*data+alphas[k]*np.log(betas[k])-gammaln(alphas[k])
        lr -= lr.max(axis=1,keepdims=True)
        r = np.exp(lr); r /= r.sum(axis=1,keepdims=True)+1e-300
        Nk = r.sum(axis=0); w = np.maximum(Nk/n,1e-10); w /= w.sum()
        for k in range(K):
            if Nk[k]<1: continue
            mk = np.dot(r[:,k],data)/Nk[k]; ml = np.dot(r[:,k],log_d)/Nk[k]
            s = np.log(max(mk,1e-10))-ml
            if s<=0: s=1e-4
            alphas[k] = max(0.1,(3-s+np.sqrt((s-3)**2+24*s))/(12*s)); betas[k]=alphas[k]/max(mk,1e-10)
        ll = sum(np.dot(r[:,k],(alphas[k]-1)*log_d-betas[k]*data)+Nk[k]*(np.log(w[k]+1e-300)+alphas[k]*np.log(betas[k])-gammaln(alphas[k])) for k in range(K))
        if abs(ll-prev)<tol: break
        prev=ll
    bic=-2*ll+(3*K-1)*np.log(n)
    return {'w':w,'Ne':alphas/betas/2,'ll':ll,'bic':bic,'n':n}

def run_desi(ts):
    n_haplo = ts.num_samples
    windows = {}
    for var in ts.variants():
        pos = int(var.site.position)
        wid = pos // WINDOW_BP
        if wid not in windows: windows[wid]=[]
        windows[wid].append(var.genotypes.astype(np.float32))
    pT=[]
    for wid,gl in windows.items():
        if len(gl)<5: continue
        g=np.array(gl); cs=g.sum(0); d=g.T@g
        het=(cs[:,None]+cs[None,:])-2*d; T=het/(2*MU*WINDOW_BP)
        for i in range(n_haplo):
            for j in range(i+1,n_haplo):
                if i//2==j//2: continue
                if T[i,j]>0: pT.append(float(T[i,j]))
    if not pT: return None
    pT=np.array(pT)
    c=930000/GEN_TIME; lo,hi=c*0.70,c*1.30
    data=pT[(pT>=lo)&(pT<=hi)]
    if len(data)<50: return {'k_hat':None,'n':len(data)}
    f1=fit_gamma_em(data,1); f2=fit_gamma_em(data,2)
    if not f1 or not f2: return None
    kh=1 if f1['bic']<=f2['bic'] else 2
    lbf=(f2['ll']-f1['ll'])/np.log(10)
    wmin=f2['w'][np.argmin(f2['Ne'])] if kh==2 else 0.0
    return {'k_hat':kh,'w_minor':wmin,'lbf':lbf,'n':len(data)}

# Model A: panmictic + bottleneck
def model_A():
    d=msprime.Demography(); d.add_population("A",initial_size=10000)
    d.add_population_parameters_change(time=813000/GEN_TIME,population="A",initial_size=1280)
    d.add_population_parameters_change(time=930000/GEN_TIME,population="A",initial_size=10000)
    return d, {"A": N_SAMPLES}

# Model D: structured (fixed API)
def model_D():
    d=msprime.Demography()
    t_split=1500000/GEN_TIME; t_admix=300000/GEN_TIME; t_bot=930000/GEN_TIME
    d.add_population("pop1",initial_size=8000)
    d.add_population("pop2",initial_size=2000)
    d.add_population("anc", initial_size=10000)
    # At admixture time, pop1 gets 20% from pop2 (model as migration event)
    # Simplify: at t_admix, merge pop2 into pop1 (80:20)
    d.add_population_parameters_change(time=t_bot*0.87, population="pop1", initial_size=1280)
    d.add_population_parameters_change(time=t_bot,      population="pop1", initial_size=8000)
    d.add_migration_rate_change(time=t_admix-1, rate=0)
    d.add_migration_rate_change(time=t_admix, rate=0)
    d.add_population_split(time=t_admix+1, derived=["pop1","pop2"], ancestral="anc")
    d.add_population_split(time=t_split,   derived=["pop1","pop2"], ancestral="anc")
    # Sample: 24 from pop1 (80%), 6 from pop2 (20%)
    n1,n2 = int(N_SAMPLES*0.8), int(N_SAMPLES*0.2)
    return d, {"pop1": n1, "pop2": n2}

print("=== DESI Power Test: Model A vs Model D ===")
for mname, mfn in [("A_panmictic_bottleneck", model_A), ("D_structured_bottleneck", model_D)]:
    print(f"\nModel {mname}:")
    k_hats=[]
    for rep in range(5):
        seed=SEED+rep*100
        try:
            dem, samples = mfn()
            ts = msprime.sim_ancestry(samples=samples, demography=dem,
                                      sequence_length=SEQ_LEN, recombination_rate=1e-8,
                                      random_seed=seed)
            ts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)
            r = run_desi(ts)
            if r and r['k_hat']:
                print(f"  Rep{rep+1}: K={r['k_hat']} w_min={r.get('w_minor',0):.3f} lbf={r.get('lbf',0):.1f} n={r['n']}")
                k_hats.append(r['k_hat'])
        except Exception as e:
            print(f"  Rep{rep+1}: ERROR {e}")
    correct_k = 1 if "panmictic" in mname else 2
    if k_hats:
        rate=sum(1 for k in k_hats if k==correct_k)/len(k_hats)
        print(f"  K={correct_k} recovery: {rate*100:.0f}% ({sum(1 for k in k_hats if k==correct_k)}/{len(k_hats)})")
