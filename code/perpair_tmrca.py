"""
Per-pair per-window TMRCA analysis
===================================
不用神经网络。直接方法：

1. 从 VCF 提取选定 pair 的每 100kb 窗口 het 计数
2. T_ij(w) = n_het / (2 * mu * L)   ← 即局部 π/2μ
3. 得到每 pair 的"跨基因组 TMRCA 分布"
4. 与模拟的 FitCoal 瓶颈 vs 无瓶颈分布对比

关键测试：
  - 若 930 ka 泛交瓶颈存在：AFR/EAS pair 的分布应在 800-1050 ka 有显著峰值
  - 若是高祖先 Ne（无瓶颈）：分布应较平滑，无 930 ka 峰
"""

import subprocess, os, sys, json
import numpy as np
from collections import defaultdict

MU        = 1.2e-8
WIN       = 100_000   # 100kb 窗口
GEN_TIME  = 28

VCF_BASE = ("/home/yanlin/public/1000GP/20220422_3202_phased_SNV_INDEL_SV/"
            "1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

OUTDIR = "/home/yanlin/popghistory/results/perpair"
os.makedirs(OUTDIR, exist_ok=True)

# ── 选定样本 ─────────────────────────────────────────────────────
YRI = ['NA18497','NA18498','NA18499','NA18517','NA18518',
       'NA18521','NA18523','NA18858','NA18879','NA18906']
CHB = ['NA18526','NA18535','NA18542','NA18548','NA18563',
       'NA18565','NA18595','NA18596','NA18605','NA18611']
JPT = ['NA18939','NA18952','NA18953','NA18961','NA18963',
       'NA18974','NA18978','NA18979','NA18984','NA18997']

ALL_SAMPLES = YRI + CHB + JPT
N = len(ALL_SAMPLES)
sample_idx = {s: i for i, s in enumerate(ALL_SAMPLES)}

# ── 定义分析的 pair 组 ───────────────────────────────────────────
pairs_within_yri  = [(YRI[i], YRI[j]) for i in range(len(YRI)) for j in range(i+1, len(YRI))]
pairs_within_chb  = [(CHB[i], CHB[j]) for i in range(len(CHB)) for j in range(i+1, len(CHB))]
pairs_within_jpt  = [(JPT[i], JPT[j]) for i in range(len(JPT)) for j in range(i+1, len(JPT))]
pairs_yri_chb     = [(y, c) for y in YRI for c in CHB]

def compute_perpair_wins(chrom, sample_list, pairs):
    """
    返回：dict[pair_tuple] -> np.array of per-window TMRCA (ka)
    """
    vcf = VCF_BASE.format(chrom=chrom)
    if not os.path.exists(vcf):
        return {}
    
    # 每个样本对应两条单倍型 (hap0 = col 0|1 的左, hap1 = 右)
    # 一对样本 (s1, s2) 是异质的当且仅当它们的单倍型组合存在不同等位基因
    # 我们用最简单的方式：s1 和 s2 的基因型不同 → heterozygous site for this pair
    
    sample_str = ",".join(sample_list)
    cmd = [
        "bcftools", "query",
        "-s", sample_str,
        "-f", "%CHROM\t%POS\t[%GT\t]\n",
        "-i", 'TYPE="SNP" && N_ALT=1',
        vcf
    ]
    
    # 初始化每 pair 的窗口 het 计数
    # 先估计染色体长度
    chr_len_cmd = ["bcftools", "view", "-H", vcf]
    
    # 用字典存 win -> {pair -> het_count}
    pair_win_het = defaultdict(lambda: defaultdict(int))
    pair_win_sites = defaultdict(int)  # callable sites per window (estimated)
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        print(f"  [{chrom}] bcftools timeout!", file=sys.stderr)
        return {}
    
    n_sites = 0
    for line in result.stdout.split('\n'):
        if not line.strip(): continue
        parts = line.split('\t')
        if len(parts) < 2 + N: continue
        
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        win_id = pos // WIN
        
        gts = parts[2:2+N]  # genotypes for each sample
        if len(gts) < N: continue
        
        # 解析每个样本的等位基因
        alleles = []
        for gt in gts:
            gt = gt.strip()
            if '|' in gt:
                a = gt.split('|')
            elif '/' in gt:
                a = gt.split('/')
            else:
                alleles.append(None)
                continue
            try:
                alleles.append((int(a[0]), int(a[1])))
            except (ValueError, IndexError):
                alleles.append(None)
        
        pair_win_sites[win_id] += 1
        
        # 对每对样本：如果两个人的等位基因不完全相同 → het pair
        for (s1, s2) in pairs:
            i1 = sample_idx.get(s1); i2 = sample_idx.get(s2)
            if i1 is None or i2 is None: continue
            a1 = alleles[i1]; a2 = alleles[i2]
            if a1 is None or a2 is None: continue
            # 任意单倍型组合中有差异就算 het
            if set(a1) != set(a2) or a1 != a2:
                pair_win_het[win_id][(s1, s2)] += 1
        
        n_sites += 1
    
    if n_sites == 0:
        return {}
    
    # 转换为 TMRCA
    result_dict = {}
    all_wins = sorted(pair_win_sites.keys())
    
    for pair in pairs:
        win_tmrca = []
        for w in all_wins:
            L = pair_win_sites[w]           # 该窗口的 SNP 数（近似可调用位点）
            if L < 5: continue              # 太少 SNP 跳过
            n_het = pair_win_het[w].get(pair, 0)
            # 用 WIN (100kb) 作为分母，而非 SNP 数
            T_ka = n_het / (2 * MU * WIN) * GEN_TIME / 1000
            win_tmrca.append(T_ka)
        if win_tmrca:
            result_dict[pair] = np.array(win_tmrca)
    
    print(f"  [{chrom}] {n_sites} SNPs, {len(all_wins)} windows", file=sys.stderr)
    return result_dict


def run_all_chroms(pairs, label, chroms=None):
    """跑所有染色体，返回合并的 per-pair TMRCA 列表"""
    if chroms is None:
        chroms = [f"chr{i}" for i in range(1, 23)]
    
    all_pair_wins = defaultdict(list)
    
    for chrom in chroms:
        d = compute_perpair_wins(chrom, ALL_SAMPLES, pairs)
        for pair, wins in d.items():
            all_pair_wins[pair].extend(wins.tolist())
    
    return dict(all_pair_wins)


# ── 主分析 ────────────────────────────────────────────────────────
if __name__ == "__main__":
    # 为了快速验证，先跑 chr21 + chr22
    test_chroms = ["chr21", "chr22"]
    
    print("=" * 60)
    print("Per-pair per-window TMRCA 分析")
    print(f"测试染色体: {test_chroms}")
    print("=" * 60)
    
    results = {}
    
    for label, pairs in [
        ("within-YRI", pairs_within_yri[:15]),     # 选15对
        ("within-CHB", pairs_within_chb[:15]),
        ("within-JPT", pairs_within_jpt[:15]),
        ("YRI×CHB",    pairs_yri_chb[:20]),
    ]:
        print(f"\n计算 {label} ({len(pairs)} pairs)...")
        pd = run_all_chroms(pairs, label, chroms=test_chroms)
        
        all_wins = []
        for wins in pd.values():
            all_wins.extend(wins)
        
        arr = np.array(all_wins)
        arr = arr[(arr > 50) & (arr < 5000)]  # 过滤异常值
        
        if len(arr) == 0:
            print(f"  {label}: 无数据")
            continue
        
        results[label] = arr.tolist()
        
        # 分布统计
        bins = [0, 200, 400, 600, 700, 800, 900, 1000, 1100, 1200, 1500, 2000, 5000]
        hist, _ = np.histogram(arr, bins=bins)
        
        print(f"  均值={np.mean(arr):.0f}ka  中位={np.median(arr):.0f}ka  "
              f"标准差={np.std(arr):.0f}ka  n窗口={len(arr):,}")
        print(f"  930ka附近（800-1050ka）: {100*np.mean((arr>800)&(arr<1050)):.1f}%")
        print(f"  < 200 ka: {100*np.mean(arr<200):.1f}%  "
              f"200-600 ka: {100*np.mean((arr>=200)&(arr<600)):.1f}%  "
              f"> 1500 ka: {100*np.mean(arr>1500):.1f}%")
    
    # 保存结果
    out = {k: v for k, v in results.items()}
    with open(f"{OUTDIR}/perpair_tmrca_chr21_22.json", "w") as f:
        json.dump({k: v[:5000] for k, v in out.items()}, f)  # 限制大小
    
    print(f"\n结果保存至 {OUTDIR}/")
