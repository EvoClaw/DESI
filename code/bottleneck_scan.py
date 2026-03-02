"""
bottleneck_scan.py
──────────────────
参数扫描：Ne_bottleneck × duration 二维网格
在 FitCoal 式的泛交单群体框架下检验 930 ka 瓶颈参数是否与 DESI 观测兼容。

核心方法：
  - 使用 msprime 仿真 + tskit branch-mode diversity（无需 sim_mutations）
  - diversity_branch / 2 * GEN_TIME / 1000 = 每窗口群组均值 TMRCA (ka)
    （与 DESI PWL 校准一致：两者都计算 E[T_pair]）
  - 与 DESI 观测的 win_t 分布直接比较

检验统计量：
  1. P(win_t_AFR > 930 ka)   ← 主证据（observed ≈ 47%）
  2. AFR 均值 TMRCA          ← 辅助（observed ≈ 1229 ka）
  3. AFR-EAS 差              ← 辅助（observed ≈ 371 ka）
"""

import os, time
import numpy as np
import msprime
import pickle

GEN_TIME  = 28
MU        = 1.2e-8
WIN_SIZE  = 100_000
SEQ_LEN   = 15_000_000   # 15 Mb，每次仿真给约 150 窗口，够用
N_HAP_AFR = 20           # AFR 单倍体样本数
N_HAP_EAS = 20           # EAS 单倍体样本数
N_SEEDS   = 4
OUTDIR    = "/home/yanlin/popghistory/results/bottleneck_scan"
os.makedirs(OUTDIR, exist_ok=True)

# ── 参数网格 ────────────────────────────────────────────────────────
# Ne_bottleneck（log 间距），包含 FitCoal 的 1280
NE_BN_GRID   = [300, 500, 800, 1280, 2000, 4000, 8000, 20000]
# 持续时间 (ka)，包含 FitCoal 的 117 ka
DUR_KA_GRID  = [30, 50, 80, 117, 150, 200, 300]
# 祖先 Ne（瓶颈外），测试两个值
ANC_NE_LIST  = [20000, 50000]

# ── 加载 DESI 观测数据（用于对比）──────────────────────────────────
print("加载 DESI 观测数据...")
AFR_IDX, EAS_IDX = 0, 1
obs_afr, obs_eas = [], []
for chrom in [f"chr{i}" for i in range(1, 23)]:
    p = f"/tmp/desi_pchmm_{chrom}.npy"
    if os.path.exists(p):
        d = np.load(p, allow_pickle=True).item()
        wt, wn = d['win_t'], d['win_n']
        m_a = (wn[:, AFR_IDX] >= 10) & ~np.isnan(wt[:, AFR_IDX])
        m_e = (wn[:, EAS_IDX] >= 10) & ~np.isnan(wt[:, EAS_IDX])
        obs_afr.extend(wt[m_a, AFR_IDX])
        obs_eas.extend(wt[m_e, EAS_IDX])

oa = np.array(obs_afr); oe = np.array(obs_eas)
OBS_P_AFR_GT930  = float(np.mean(oa > 930))
OBS_AFR_MEAN     = float(np.mean(oa))
OBS_DIFF         = float(np.mean(oa) - np.mean(oe))
print(f"  观测 P(AFR>930)  = {OBS_P_AFR_GT930:.3f}  ({OBS_P_AFR_GT930*100:.1f}%)")
print(f"  观测 AFR 均值    = {OBS_AFR_MEAN:.0f} ka")
print(f"  观测 AFR-EAS 差  = {OBS_DIFF:.0f} ka")

# ── 核心仿真函数 ─────────────────────────────────────────────────────
def simulate_scenario(ne_bn, dur_ka, anc_ne, seed):
    """
    仿真 OoA + FitCoal 式瓶颈，返回 per-window 群组均值 TMRCA。
    使用 tskit diversity(mode="branch")，比逐树循环快 10x。
    """
    demog = msprime.Demography()
    demog.add_population(name="AFR", initial_size=15000)
    demog.add_population(name="EAS", initial_size=3000)
    demog.add_population(name="ANC", initial_size=anc_ne)
    # OoA 分化
    demog.add_population_split(
        time=65000 // GEN_TIME,
        derived=["AFR", "EAS"], ancestral="ANC"
    )
    # 瓶颈：813 ka 结束，813+dur_ka ka 开始
    t_end_gen   = int(813_000 / GEN_TIME)
    t_start_gen = int((813_000 + dur_ka * 1_000) / GEN_TIME)
    demog.add_population_parameters_change(time=t_end_gen,   population="ANC", initial_size=ne_bn)
    demog.add_population_parameters_change(time=t_start_gen, population="ANC", initial_size=anc_ne)

    ts = msprime.sim_ancestry(
        samples={"AFR": N_HAP_AFR // 2, "EAS": N_HAP_EAS // 2},
        demography=demog,
        sequence_length=SEQ_LEN,
        recombination_rate=1e-8,
        random_seed=seed,
    )

    # 样本索引
    afr_samples = list(ts.samples(population=0))   # 20 AFR haplotypes
    eas_samples = list(ts.samples(population=1))   # 20 EAS haplotypes

    # 按窗口计算 branch diversity = 2 * E[T_MRCA] (单位：代)
    n_wins   = int(SEQ_LEN // WIN_SIZE)
    windows  = [i * WIN_SIZE for i in range(n_wins + 1)]

    # ts.diversity returns shape (n_windows,) for one sample set
    # DO NOT index with [0] — that would return only the first window!
    div_afr = np.array(ts.diversity(
        sample_sets=[afr_samples], windows=windows, mode="branch"
    )).ravel()   # shape: (n_wins,)
    div_eas = np.array(ts.diversity(
        sample_sets=[eas_samples], windows=windows, mode="branch"
    )).ravel()

    # branch diversity = 2 * E[T_pair_gen]
    # T_MRCA_ka = E[T_pair_gen] * GEN_TIME / 1000
    t_afr_ka = div_afr / 2 * GEN_TIME / 1000
    t_eas_ka = div_eas / 2 * GEN_TIME / 1000

    # 去掉空窗口
    valid = (div_afr > 0) & (div_eas > 0)
    return t_afr_ka[valid], t_eas_ka[valid]


# ── 参数扫描主循环 ───────────────────────────────────────────────────
print("\n开始参数扫描...")
total = len(NE_BN_GRID) * len(DUR_KA_GRID) * len(ANC_NE_LIST) * N_SEEDS
done  = 0
results = []   # list of dicts

for anc_ne in ANC_NE_LIST:
    for ne_bn in NE_BN_GRID:
        for dur_ka in DUR_KA_GRID:
            afr_pooled, eas_pooled = [], []
            t0 = time.time()
            for seed in range(1, N_SEEDS + 1):
                a, e = simulate_scenario(ne_bn, dur_ka, anc_ne, seed)
                afr_pooled.extend(a); eas_pooled.extend(e)
                done += 1

            ta = np.array(afr_pooled); te = np.array(eas_pooled)
            row = {
                "anc_ne":   anc_ne,
                "ne_bn":    ne_bn,
                "dur_ka":   dur_ka,
                "afr_mean": float(np.mean(ta)),
                "eas_mean": float(np.mean(te)),
                "p_afr_gt930": float(np.mean(ta > 930)),
                "p_eas_gt930": float(np.mean(te > 930)),
                "diff_mean":   float(np.mean(ta) - np.mean(te)),
                "afr_cv":      float(np.std(ta) / np.mean(ta)),
                "n_windows":   len(ta),
                # KS-like: mean abs deviation from observed P(>930)
                "delta_p_gt930": abs(float(np.mean(ta > 930)) - OBS_P_AFR_GT930),
                "delta_afr_mean": abs(float(np.mean(ta)) - OBS_AFR_MEAN),
                "delta_diff": abs(float(np.mean(ta) - np.mean(te)) - OBS_DIFF),
            }
            results.append(row)
            elapsed = time.time() - t0
            print(f"  anc={anc_ne//1000}k  Ne_bn={ne_bn:>6}  dur={dur_ka:>4}ka  "
                  f"AFR={row['afr_mean']:>7.0f}ka  "
                  f"P(>930)={row['p_afr_gt930']:.2f}  "
                  f"Δ_obs={row['delta_p_gt930']:.2f}  "
                  f"({elapsed:.1f}s)")

# ── 保存结果 ─────────────────────────────────────────────────────────
out_path = os.path.join(OUTDIR, "scan_results.pkl")
with open(out_path, "wb") as f:
    pickle.dump({
        "results":        results,
        "obs_p_afr_gt930": OBS_P_AFR_GT930,
        "obs_afr_mean":    OBS_AFR_MEAN,
        "obs_diff":        OBS_DIFF,
        "ne_bn_grid":     NE_BN_GRID,
        "dur_ka_grid":    DUR_KA_GRID,
        "anc_ne_list":    ANC_NE_LIST,
    }, f)
print(f"\n结果保存至 {out_path}")

# ── 快速总结 ─────────────────────────────────────────────────────────
import pandas as pd
df = pd.DataFrame(results)

print("\n" + "=" * 75)
print("Pass/Fail 判定（兼容区：|P(AFR>930) - obs| < 0.10 且 |AFR均值 - obs| < 150ka）")
print("=" * 75)
compat = df[(df.delta_p_gt930 < 0.10) & (df.delta_afr_mean < 150)]
print(f"兼容组合数：{len(compat)} / {len(df)}")
if len(compat) > 0:
    print(compat[["anc_ne","ne_bn","dur_ka","afr_mean","p_afr_gt930","delta_p_gt930"]].to_string(index=False))

# FitCoal 具体点
fitcoal = df[(df.anc_ne == 50000) & (df.ne_bn == 1280) & (df.dur_ka == 117)]
if len(fitcoal):
    print("\nFitCoal 具体参数 (anc=50k, Ne=1280, dur=117ka):")
    print(fitcoal[["afr_mean","p_afr_gt930","diff_mean","delta_p_gt930","delta_afr_mean"]].to_string(index=False))
else:
    print("\n(FitCoal Ne=1280 + dur=117 ka 的点可能已通过最接近组合估计)")

print("\n完成。")
