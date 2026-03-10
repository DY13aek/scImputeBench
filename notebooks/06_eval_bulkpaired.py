# Set the base directory for datasets
# Change this path to point to your local data directory
BASE_DIR = "../data"

# ====== Cell Type ======
import os
import scanpy as sc
import pandas as pd

os.chdir(f"{BASE_DIR}/Chu")
chu_raw = sc.read_csv("./GSE75748_sc_cell_type_ec.csv.gz")
chu_raw = chu_raw.T.copy()

chu_raw.obs_names.name = "cell_id"
chu_raw.var_names.name = "gene_symbol"

def parse_cell_type(cell_id):
    if cell_id.startswith("H1_"):
        return "H1_ES"
    if cell_id.startswith("H9_"):
        return "H9_ES"
    if cell_id.startswith("DEC"):
        return "DEC"
    if cell_id.startswith("EC"):
        return "EC"
    if cell_id.startswith("NPC"):
        return "NPC"
    if cell_id.startswith("TB"):
        return "TB"
    if cell_id.startswith("HFF"):
        return "HFF"
    return "unknown"

chu_raw.obs["cell_type"] = [parse_cell_type(i) for i in chu_raw.obs_names]

# ====================================================
# ====================================================
# 0. 데이터 복사 (원본 보존)
adata = chu_raw.copy()
# QC 메트릭 계산 (n_genes_by_counts 등을 얻음)
sc.pp.calculate_qc_metrics(adata, inplace=True)
# 전체 유전자 수 대비 발현된 유전자 비율 계산
adata.obs['frac_genes'] = adata.obs['n_genes_by_counts'] / adata.n_vars

# ====================================================
# 1. 기본 필터링 (품질 낮은 세포/유전자 제거)
# 참고: Huang et al. 2025 기준 (min_genes=200, min_cells=10)
# 1. 저장할 경로 설정
save_dir = f'{BASE_DIR}/Chu/Raw'
# 2. Scanpy의 그림 저장 경로 지정
sc.settings.figdir = save_dir
# 세포 필터링: 유전자가 200개 미만으로 발현된 죽은 세포 제거
sc.pp.filter_cells(adata, min_genes=200)
# 유전자 필터링: 3개 미만의 세포에서만 발견된 희귀 유전자 제거
sc.pp.filter_genes(adata, min_cells=3)

print(f"기본 필터링 후 크기: {adata.shape}") #(1018, 19097) -> (1018, 17559)

# ====================================================
# 2. 품질 관리 (QC) 메트릭 계산
# 미토콘드리아(MT) 유전자 비율 등을 계산하여 스트레스 받은 세포 식별
# ====================================================
# 미토콘드리아 유전자 식별 (Human 데이터는 보통 'MT-'로 시작)
adata.var['mt'] = adata.var_names.str.startswith('MT-') 

# QC 메트릭 계산
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 3. 그림 그리기 (show=False, save='파일식별자.png')
# save 매개변수에 문자열을 넣으면 'figures/violin' + '파일식별자.png' 형태로 저장됩니다.
# 주의: Scanpy는 기본적으로 figdir 아래 'figures' 폴더를 한 번 더 만듭니다.
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
             jitter=0.4, multi_panel=True, 
             show=False, save='_qc_metrics.png')
adata = adata[adata.obs['pct_counts_mt'] < 25, :] #(1006 × 17559)

adata.write(os.path.join(save_dir, 'chu.h5ad'))

# ====== Time Course ======
import os
import scanpy as sc
import pandas as pd

os.chdir(f"{BASE_DIR}/Chu_time")
chu_raw = sc.read_csv("./GSE75748_sc_time_course_ec.csv.gz")
chu_raw = chu_raw.T.copy()

chu_raw.obs_names.name = "cell_id"
chu_raw.var_names.name = "gene_symbol"

# 1. Index에서 시간 정보 추출 (예: 'H9.00hb4s_001' -> '00h')
# '.'으로 나누고 첫 3글자(00h, 12h 등)만 가져오거나, 정규식을 씁니다.
time_points = [idx.split('.')[1][:3] for idx in chu_raw.obs_names]

# 2. obs['cell_type'] 컬럼에 할당
chu_raw.obs['cell_type'] = time_points

# 3. 제대로 들어갔는지 확인
print(chu_raw.obs['cell_type'].value_counts())

# ====================================================
# 0. 데이터 복사 (원본 보존)
adata = chu_raw.copy()
# QC 메트릭 계산 (n_genes_by_counts 등을 얻음)
sc.pp.calculate_qc_metrics(adata, inplace=True)
# 전체 유전자 수 대비 발현된 유전자 비율 계산
adata.obs['frac_genes'] = adata.obs['n_genes_by_counts'] / adata.n_vars

# ====================================================
# 1. 기본 필터링 (품질 낮은 세포/유전자 제거)
# 참고: Huang et al. 2025 기준 (min_genes=200, min_cells=10)
# 1. 저장할 경로 설정
save_dir = f'{BASE_DIR}/Chu_time'
# 2. Scanpy의 그림 저장 경로 지정
sc.settings.figdir = save_dir
# 세포 필터링: 유전자가 200개 미만으로 발현된 죽은 세포 제거
sc.pp.filter_cells(adata, min_genes=200)
# 유전자 필터링: 3개 미만의 세포에서만 발견된 희귀 유전자 제거
sc.pp.filter_genes(adata, min_cells=3)

print(f"기본 필터링 후 크기: {adata.shape}") #(1018, 19097) -> (1018, 17559)

# ====================================================
# 2. 품질 관리 (QC) 메트릭 계산
# 미토콘드리아(MT) 유전자 비율 등을 계산하여 스트레스 받은 세포 식별
# ====================================================
# 미토콘드리아 유전자 식별 (Human 데이터는 보통 'MT-'로 시작)
adata.var['mt'] = adata.var_names.str.startswith('MT-') 

# QC 메트릭 계산
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# 3. 그림 그리기 (show=False, save='파일식별자.png')
# save 매개변수에 문자열을 넣으면 'figures/violin' + '파일식별자.png' 형태로 저장됩니다.
# 주의: Scanpy는 기본적으로 figdir 아래 'figures' 폴더를 한 번 더 만듭니다.
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
             jitter=0.4, multi_panel=True, 
             show=False, save='_qc_metrics.png')
adata = adata[adata.obs['pct_counts_mt'] < 25, :] #(1006 × 17559)

adata.write(os.path.join(save_dir, 'chu.h5ad'))

# ==========scRNA-seq==========
# 필요한 라이브러리 임포트
import os
import scanpy as sc
import numpy as np
import pandas as pd

os.chdir(f"{BASE_DIR}/Chu/Raw")
raw = sc.read_h5ad("chu.h5ad")
cell_idx, gene_idx = raw.obs.index.tolist(), raw.var.index.tolist()
#raw.raw = raw.copy()  # 원본 데이터 보존을 위해 raw에 복사본 저장
#sc.pp.normalize_total(raw, target_sum=1e4)
#sc.pp.log1p(raw)
#raw.write_h5ad("chu.h5ad")  # 정규화 및 로그 변환된 데이터를 저장

# imputed datasets
os.chdir(f"{BASE_DIR}/Chu/imputed")
# === 1. DCA : Done ===
#dca_df = pd.read_csv("./DCA/mean.tsv", sep="\t", header=0, index_col=0).T
#dca = sc.AnnData(dca_df)
dca = sc.read_h5ad("./DCA/dca.h5ad")
#dca.obs, dca.var = raw.obs.copy(), raw.var.loc[dca.var_names].copy()
#dca.raw = dca.copy()  # 원본 데이터 보존을 위해 raw에 복사본 저장
#sc.pp.normalize_total(dca, target_sum=1e4)
#sc.pp.log1p(dca)
#dca.write_h5ad("./DCA/dca.h5ad")# 최종 저장 및 불러오기

# === 2. DrImpute: Done ===
#drimpute_df = pd.read_csv("./DrImpute/chu_DrImpute.csv", sep=",", index_col=0).T
#drimpute = sc.AnnData(drimpute_df)
#drimpute.obs = raw.obs.copy()
#drimpute.write_h5ad("./DrImpute/drimpute.h5ad")
drimpute = sc.read_h5ad("./DrImpute/drimpute.h5ad")

# === 3. MAGIC: Done ===
magic = sc.read_h5ad("./MAGIC/chu_MAGIC.h5ad")

# === 4. SAVER: Done ===
#saver_df = pd.read_csv("./SAVER/chu_saver_imputed.csv", sep=",", index_col=0).T
#saver = sc.AnnData(saver_df)
#saver.obs, saver.var = raw.obs.copy(), raw.var.copy()
saver = sc.read_h5ad("./SAVER/saver.h5ad")
#saver.raw = saver.copy()  # 원본 데이터 보존을 위해 raw에 복사본 저장
#sc.pp.normalize_total(saver, target_sum=1e4)
#sc.pp.log1p(saver)
#saver.write_h5ad("./SAVER/saver.h5ad")

# === 5. scIDPMs  ===
#scidpms_df = pd.read_csv("./scIDPMs/imputed.csv", sep=",", index_col=None, header=None)
#scidpms = sc.AnnData(scidpms_df)
#scidpms.obs, scidpms.var = raw.obs.copy(), raw.var.copy()
scidpms = sc.read_h5ad("./scIDPMs/scidpms.h5ad")
#scidpms.raw = scidpms.copy()  # 원본 데이터 보존을 위해 raw에 복사본 저장
#sc.pp.normalize_total(scidpms, target_sum=1e4)
#sc.pp.log1p(scidpms)
#scidpms.write_h5ad("./scIDPMs/scidpms.h5ad")

# === 6. scIGANs ===
#scigans_df = pd.read_csv("./scIGANs/scIGANs_chu.txt", sep="\t", index_col=0).T
#scigans = sc.AnnData(scigans_df)
#scigans.obs, scigans.var = raw.obs.copy(), raw.var.copy()
scigans = sc.read_h5ad("./scIGANs/scigans.h5ad")
#scigans.raw = scigans.copy()  # 원본 데이터 보존을 위해 raw에 복사본 저장
#sc.pp.normalize_total(scigans, target_sum=1e4)
#sc.pp.log1p(scigans)
#scigans.write_h5ad("./scIGANs/scigans.h5ad")

# === 7. scMultiGAN: Done ===
#scmultigan_df = pd.read_csv("./scMultiGAN/Restored/chu_scMultiGAN.tsv", sep="\t", header=None, index_col=None).T
#scmultigan = sc.AnnData(scmultigan_df)
#scmultigan.obs, scmultigan.var = raw.obs.copy(), raw.var.copy()
scmultigan = sc.read_h5ad("./scMultiGAN/scmultigan.h5ad")
#scmultigan.raw = scmultigan.copy()  # 원본 데이터 보존을 위해 raw에 복사본 저장
#sc.pp.normalize_total(scmultigan, target_sum=1e4)
#sc.pp.log1p(scmultigan)
#scmultigan.write_h5ad("./scMultiGAN/scmultigan.h5ad")

# === 8. scSTD: Done ===
#scstd_df = pd.read_csv("./scSTD/Result/chu_imputed.txt", sep=",", index_col=None, header=None)
#scstd = sc.AnnData(scstd_df)
#scstd.obs, scstd.var = raw.obs.copy(), raw.var.copy()
scstd = sc.read_h5ad("./scSTD/scstd.h5ad")
#scstd.raw = scstd.copy()  # 원본 데이터 보존을 위해 raw에 복사본 저장
#sc.pp.normalize_total(scstd, target_sum=1e4)
#sc.pp.log1p(scstd)
#scstd.write_h5ad("./scSTD/scstd.h5ad")

# === 9. scVI: Done ===
scvi = sc.read_h5ad("./scVI/scvi.h5ad")
scvi.obs_names, scvi.var_names = cell_idx, gene_idx
print("scvi shape:", scvi.shape)

# === 10. ALRA: Done ===
#alra_df = pd.read_csv("./ALRA/chu_ALRA_imputed.csv", index_col=0)
#alra = sc.AnnData(alra_df)
#alra.obs, alra.var = raw.obs.copy(), raw.var.copy()
#alra.obs_names, alra.var_names = cell_idx, gene_idx
#alra.write_h5ad("./ALRA/alra.h5ad")
alra = sc.read_h5ad("./ALRA/alra.h5ad")


# ==========bulkRNA-seq==========
os.chdir(f"{BASE_DIR}/Chu/")
bulk = sc.read_csv("./GSE75748_bulk_cell_type_ec.csv.gz").T

chu_data = {
    "RAW"       : raw,
    "ALRA"      : alra,
    "DCA"       : dca,
    "DrImpute"  : drimpute,
    "MAGIC"     : magic,
    "SAVER"     : saver,
    "scIDPMs"   : scidpms,
    "scIGANs"   : scigans,
    "scMultiGAN": scmultigan,
    "scSTD"     : scstd,
    "scVI"      : scvi,
}

bulk = 

# === Bulk 데이터에서 DESeq2로 DEG 추출 (환경: commonR)===
# ===================================================================
# Bulk DEG Analysis - DESeq2
# Contrast: DEC vs H1 (positive LFC = higher in DEC)
# ===================================================================

library(DESeq2)
library(tidyverse)

setwd(f"{BASE_DIR}/Chu")

# ---------- 데이터 로드 ----------
data_path <- f"{BASE_DIR}/Chu/GSE75748_bulk_cell_type_ec.csv.gz"
raw_counts <- read.csv(data_path, row.names = 1, check.names = FALSE)

# ---------- H1, DEC 샘플 추출 ----------
target_samples <- c("H1_rep1", "H1_rep2", "H1_rep3", "H1_rep4", "DEC_rep1", "DEC_rep2")
counts_subset <- raw_counts[, target_samples]

cat("Count matrix dimensions:", dim(counts_subset), "\n")
cat("Library sizes:\n")
print(colSums(counts_subset))

# ---------- 메타데이터 생성 ----------
sample_info <- data.frame(
  condition = factor(c(rep("H1", 4), rep("DEC", 2)), levels = c("H1", "DEC")),
  row.names = target_samples
)

# ---------- Low count 필터 ----------
# 최소 2개 샘플에서 10 이상 발현된 유전자만 유지
keep <- rowSums(counts_subset >= 10) >= 2
counts_filtered <- counts_subset[keep, ]
cat("Genes after filtering:", nrow(counts_filtered), "/", nrow(counts_subset), "\n")

# ---------- DESeq2 ----------
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered),
  colData   = sample_info,
  design    = ~ condition
)

dds <- DESeq(dds)

cat("\nSize factors:\n")
print(sizeFactors(dds))

# ---------- 결과 추출 (DEC vs H1) ----------
# log2FC > 0  →  DEC에서 높음
# log2FC < 0  →  H1에서 높음
res <- results(dds,
               contrast = c("condition", "DEC", "H1"),
               alpha    = 0.05)

cat("\nDESeq2 Results Summary:\n")
summary(res)

# ---------- 결과 정리 ----------
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  mutate(
    significant = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
    diff_status = case_when(
      significant & log2FoldChange > 0  ~ "UP",    # DEC > H1
      significant & log2FoldChange < 0  ~ "DOWN",  # H1 > DEC
      TRUE                               ~ "NS"
    )
  ) %>%
  arrange(desc(abs(log2FoldChange)))

# ---------- 요약 ----------
cat("\nDEG counts:\n")
print(table(res_df$diff_status))

cat("\nTop 10 DEGs by |log2FC|:\n")
print(res_df %>%
        filter(significant) %>%
        head(10) %>%
        select(gene, log2FoldChange, padj, diff_status))

# ---------- 저장 ----------
write.csv(res_df,
          "./DEG/BULK_DEG_H1_vs_DEC_results.csv",
          row.names = FALSE)

cat("\n✅ Analysis complete! Results saved to BULK_DEG_H1_vs_DEC_results.csv\n")

# Single-cell DEG Analysis - Wilcoxon Rank-Sum Test
# Contrast: DEC vs H1 (positive LFC = higher in DEC)
#           → bulk DESeq2와 방향 통일
# Input:    log1p normalized data (sc.pp.normalize_total + sc.pp.log1p)
# Filter:   min_pct=0.1 (dropout 노이즈 제거)
# ===================================================================

import numpy as np
import pandas as pd
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import os

# ===================================================================
# DEG 함수
# ===================================================================

def wilcoxon_deg_analysis(adata, group_col, group1, group2,
                           pval_threshold=0.05, min_pct=0.1):
    """
    Wilcoxon rank-sum test for DEG analysis on log1p-normalized data.
    
    Parameters
    ----------
    adata         : AnnData  (log1p normalized)
    group_col     : str      obs column with group labels
    group1        : str      reference group   (e.g. 'H1_ES')  → denominator of LFC
    group2        : str      comparison group  (e.g. 'DEC')    → numerator   of LFC
    pval_threshold: float    FDR threshold for significance flag
    min_pct       : float    min fraction of cells expressing gene
                             in at least one group (dropout filter)
    
    Returns
    -------
    pd.DataFrame sorted by |log2FC| descending
    log2FC > 0  →  higher in group2 (DEC)
    log2FC < 0  →  higher in group1 (H1_ES)
    """
    
    mask1 = adata.obs[group_col] == group1   # H1_ES
    mask2 = adata.obs[group_col] == group2   # DEC
    
    print(f"  {group1}: {mask1.sum()} cells")
    print(f"  {group2}: {mask2.sum()} cells")
    
    expr1 = adata[mask1].X
    expr2 = adata[mask2].X
    
    if hasattr(expr1, 'toarray'):
        expr1 = expr1.toarray()
    if hasattr(expr2, 'toarray'):
        expr2 = expr2.toarray()
    
    # ------------------------------------------------------------------
    # pct 필터: 두 그룹 중 하나라도 min_pct 이상 발현된 유전자만 테스트
    # 효과:
    #   - Raw  : dropout 많음 → 필터에서 많이 걸림 (테스트 유전자 ↓)
    #   - Imputed : dropout 채워짐 → 더 많은 유전자 통과 (테스트 유전자 ↑)
    #   → 이 차이 자체가 imputation 효과를 보여주는 결과가 됨
    # ------------------------------------------------------------------
    pct1 = (expr1 > 0).mean(axis=0)   # shape: (n_genes,)
    pct2 = (expr2 > 0).mean(axis=0)
    expressed_mask = (pct1 > min_pct) | (pct2 > min_pct)
    expressed_idx  = np.where(expressed_mask)[0]
    
    print(f"  Genes after pct>{min_pct} filter: {expressed_mask.sum()} / {adata.n_vars}")
    
    results = []    
    
    for i in expressed_idx:
        gene       = adata.var_names[i]
        g1         = expr1[:, i]
        g2         = expr2[:, i]
    
        # Wilcoxon rank-sum test
        try:
            stat, pval = ranksums(g1, g2)
        except Exception:
            pval, stat = 1.0, 0.0
    
        # Log2FC: group2 / group1  →  DEC / H1_ES (bulk와 방향 동일)
        # expm1: log1p → count-like 복원 후 평균 → pseudo-count +1 로 0 방지
        mean2 = np.mean(np.expm1(g2))
        mean1 = np.mean(np.expm1(g1))
        log2fc = np.log2((mean2 + 1) / (mean1 + 1))
    
        results.append({
            'gene'          : gene,
            'pvalue'        : pval,
            'log2FoldChange': log2fc,
            'pct1'          : pct1[i],   # H1_ES 발현 비율
            'pct2'          : pct2[i],   # DEC   발현 비율
            'mean_H1'       : np.mean(g1),
            'mean_DEC'      : np.mean(g2),
        })
    
    results_df = pd.DataFrame(results)
    
    # FDR correction (Benjamini-Hochberg)
    _, results_df['padj'], _, _ = multipletests(results_df['pvalue'], method='fdr_bh')
    
    results_df['significant'] = results_df['padj'] < pval_threshold
    
    # |log2FC| 내림차순 정렬 (bulk DESeq2와 동일 정렬 기준)
    results_df = results_df.sort_values('log2FoldChange',
                                        key=lambda x: np.abs(x),
                                        ascending=False).reset_index(drop=True)
    
    return results_df


# ===================================================================
# 데이터 준비
# ===================================================================

chu_data = {
    "RAW"       : raw,
    "ALRA"      : alra,
    "DCA"       : dca,
    "DrImpute"  : drimpute,
    "MAGIC"     : magic,
    "SAVER"     : saver,
    "scIDPMs"   : scidpms,
    "scIGANs"   : scigans,
    "scMultiGAN": scmultigan,
    "scSTD"     : scstd,
    "scVI"      : scvi,
}

# ===================================================================
# 실행
# ===================================================================

print("=" * 60)
print("Single-cell DEG Analysis (DEC vs H1_ES) - Wilcoxon")
print("log2FC > 0  →  higher in DEC  (bulk와 방향 동일)")
print("=" * 60)

all_deg_results = []

for method_name, adata in chu_data.items():
    print(f"\n{'=' * 60}")
    print(f"Method: {method_name}")
    print(f"{'=' * 60}")
    
    deg_result = wilcoxon_deg_analysis(
        adata,
        group_col = 'cell_type',
        group1    = 'H1_ES',   # denominator
        group2    = 'DEC',     # numerator  → LFC 방향: DEC/H1 = bulk와 동일
        min_pct   = 0.1
    )
    
    deg_result['Method'] = method_name
    
    sig   = deg_result['significant']
    up    = (sig & (deg_result['log2FoldChange'] > 0)).sum()
    down  = (sig & (deg_result['log2FoldChange'] < 0)).sum()
    
    print(f"\n  Results:")
    print(f"    Total significant DEGs : {sig.sum()}")
    print(f"    UP   (DEC > H1_ES)     : {up}")
    print(f"    DOWN (H1_ES > DEC)     : {down}")
    
    all_deg_results.append(deg_result)

# ===================================================================
# 저장
# ===================================================================

os.chdir(f"{BASE_DIR}/Chu/DEG")

combined_deg = pd.concat(all_deg_results, ignore_index=True)

combined_deg = combined_deg[[
    'Method', 'gene', 'log2FoldChange', 'pvalue', 'padj',
    'significant', 'pct1', 'pct2', 'mean_H1', 'mean_DEC'
]]

combined_deg.to_csv("SC_DEG_All_Methods_DEC_vs_H1_Wilcoxon.csv", index=False)

# ===================================================================
# 요약 통계
# ===================================================================

print("\n" + "=" * 80)
print("SUMMARY STATISTICS")
print("=" * 80)

summary = combined_deg.groupby('Method').agg(
    Total_genes   = ('gene', 'count'),
    DEG_count     = ('significant', 'sum'),
    Genes_tested  = ('gene', 'count'),   # pct 필터 통과 유전자 수
).copy()

up_counts   = combined_deg[combined_deg['significant'] & (combined_deg['log2FoldChange'] > 0)].groupby('Method').size()
down_counts = combined_deg[combined_deg['significant'] & (combined_deg['log2FoldChange'] < 0)].groupby('Method').size()

summary['UP_count']   = up_counts
summary['DOWN_count'] = down_counts
summary['DEG_ratio']  = (summary['DEG_count'] / summary['Total_genes'] * 100).round(2)

print(summary[['Total_genes', 'DEG_count', 'DEG_ratio', 'UP_count', 'DOWN_count']])

summary.to_csv("SC_DEG_Summary_Statistics.csv")

print("\n✅ All SC DEG analyses complete!")
print("✅ Combined results : SC_DEG_All_Methods_DEC_vs_H1_Wilcoxon.csv")
print("✅ Summary          : SC_DEG_Summary_Statistics.csv")

# DEG Comparison Metrics - FINAL
# bulk : DESeq2 Wald stat 랭킹
# SC   : Wilcoxon + pct>0.1 필터 + LFC 랭킹
#
# 최종 메트릭 5종:
#   1. F1                  (DEG 목록 일치도 + false positive)
#   2. SCC                 (전체 LFC 패턴 상관)
#   3. AUPRC               (랭킹 품질, threshold-free)
#   4. Top-K Overlap       (K = 사용자 지정)
#   5. Top-K Direction     (K = 사용자 지정)
# ===================================================================

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score
import os

os.chdir(f"{BASE_DIR}/Chu/DEG")

# ===================================================================
# ★ 여기서 Top-K 값을 자유롭게 지정하세요 ★
# ===================================================================

TOP_K_VALUES = [100, 200, 300, 400, 500]   # 원하는 K 값 목록

# ===================================================================
# 데이터 로드
# ===================================================================

bulk_deg = pd.read_csv("BULK_DEG_H1_vs_DEC_results.csv")
sc_deg   = pd.read_csv("SC_DEG_All_Methods_DEC_vs_H1_Wilcoxon.csv")

# rank score
bulk_deg['rank_score'] = bulk_deg['stat']          # Wald stat
sc_deg['rank_score']   = sc_deg['log2FoldChange']  # LFC

# bulk significant 정보
bulk_sig_genes = set(bulk_deg[bulk_deg['significant']]['gene'])
bulk_dir_dict  = (bulk_deg[bulk_deg['significant']]
                  .set_index('gene')['log2FoldChange']
                  .apply(lambda x: 'UP' if x > 0 else 'DOWN').to_dict())

print(f"Bulk significant DEGs: {len(bulk_sig_genes)} (padj<0.05, |log2FC|>1)")
print(f"Top-K values         : {TOP_K_VALUES}")

# ===================================================================
# Metric 1 & 2: F1 + SCC
# ===================================================================

def calculate_f1_scc(bulk_deg, sc_method, bulk_sig_genes):
    sc_sig_genes  = set(sc_method[sc_method['significant']]['gene'])
    overlap_count = len(bulk_sig_genes & sc_sig_genes)
    recall    = overlap_count / len(bulk_sig_genes) if bulk_sig_genes else 0
    precision = overlap_count / len(sc_sig_genes)   if sc_sig_genes  else 0
    f1        = (2 * precision * recall / (precision + recall)
                 if (precision + recall) > 0 else 0)
    merged = (bulk_deg[['gene', 'log2FoldChange']]
              .merge(sc_method[['gene', 'log2FoldChange']],
                     on='gene', suffixes=('_bulk', '_sc'))
              .dropna())
    scc, _ = (spearmanr(merged['log2FoldChange_bulk'], merged['log2FoldChange_sc'])
               if len(merged) > 1 else (np.nan, np.nan))
    return {
        'SC_DEGs'  : len(sc_sig_genes),
        'Overlap'  : overlap_count,
        'Recall'   : recall,
        'Precision': precision,
        'F1'       : f1,
        'SCC'      : scc,
    }

# ===================================================================
# Metric 3: AUPRC
# ===================================================================

def calculate_auprc(bulk_deg, sc_method):
    merged = (bulk_deg[['gene', 'significant']]
              .merge(sc_method[['gene', 'rank_score']], on='gene')
              .dropna())
    if len(merged) > 1 and merged['significant'].sum() > 0:
        auprc = average_precision_score(
            merged['significant'].astype(int),
            merged['rank_score'].abs()
        )
    else:
        auprc = np.nan
    return {'AUPRC': auprc}

# ===================================================================
# Metric 4 & 5: Top-K Overlap + Top-K Direction
# ===================================================================

def calculate_topk(bulk_deg, sc_method, bulk_dir_dict, k_values):
    common_genes = set(bulk_deg['gene']) & set(sc_method['gene'])
    bulk_c = bulk_deg[bulk_deg['gene'].isin(common_genes)].copy()
    sc_c   = sc_method[sc_method['gene'].isin(common_genes)].copy()
    sc_dir_dict = (sc_c.set_index('gene')['log2FoldChange']
                       .apply(lambda x: 'UP' if x > 0 else 'DOWN').to_dict())
    out = {'Common_genes': len(common_genes)}
    for k in k_values:
        eff_k = min(k, len(common_genes))
        bulk_topk_genes = set(bulk_c.loc[
            bulk_c['rank_score'].abs().nlargest(eff_k).index, 'gene'])
        sc_topk_genes   = set(sc_c.loc[
            sc_c['rank_score'].abs().nlargest(eff_k).index, 'gene'])
        overlap_genes = bulk_topk_genes & sc_topk_genes
        overlap_count = len(overlap_genes)
        dir_match = sum(
            1 for g in overlap_genes
            if bulk_dir_dict.get(g) == sc_dir_dict.get(g)
        )
        dir_concordance = dir_match / overlap_count if overlap_count > 0 else np.nan
        out[f'Top{k}_Overlap']   = overlap_count
        out[f'Top{k}_Ratio']     = round(overlap_count / eff_k, 4)
        out[f'Top{k}_Direction'] = round(dir_concordance, 4)
    return out

# ===================================================================
# 실행
# ===================================================================

print("\n" + "=" * 80)
print("Final DEG Comparison  |  bulk: Wald stat  /  SC: LFC")
print("=" * 80)

all_results = []

for method_name in sc_deg['Method'].unique():
    sc_method = sc_deg[sc_deg['Method'] == method_name].copy()
    row = {'Method': method_name}
    row.update(calculate_f1_scc(bulk_deg, sc_method, bulk_sig_genes))
    row.update(calculate_auprc(bulk_deg, sc_method))
    row.update(calculate_topk(bulk_deg, sc_method, bulk_dir_dict, TOP_K_VALUES))
    all_results.append(row)
    topk_str = " | ".join(
        f"Top{k}={row[f'Top{k}_Overlap']}" for k in TOP_K_VALUES
    )
    print(f"  {method_name:<12} F1={row['F1']:.3f} | SCC={row['SCC']:.3f} | "
          f"AUPRC={row['AUPRC']:.3f} | {topk_str}")

# ===================================================================
# 결과 정리
# ===================================================================

results_df = pd.DataFrame(all_results)

topk_cols = []
for k in TOP_K_VALUES:
    topk_cols += [f'Top{k}_Overlap', f'Top{k}_Ratio', f'Top{k}_Direction']

col_order = (
    ['Method', 'SC_DEGs', 'Overlap', 'Recall', 'Precision', 'F1',
     'SCC', 'AUPRC', 'Common_genes']
    + topk_cols
)
results_df = results_df[col_order]

# ===================================================================
# Tool 랭킹
# ===================================================================

rank_metrics = (
    ['F1', 'SCC', 'AUPRC']
    + [f'Top{k}_Overlap'   for k in TOP_K_VALUES]
    + [f'Top{k}_Direction' for k in TOP_K_VALUES]
)

rank_df = results_df[['Method']].copy()
for m in rank_metrics:
    rank_df[f'rank_{m}'] = (results_df[m]
                             .rank(ascending=False, method='min')
                             .astype(int))

rank_df['Mean_Rank']  = rank_df[[f'rank_{m}' for m in rank_metrics]].mean(axis=1)
rank_df['Final_Rank'] = rank_df['Mean_Rank'].rank(ascending=True, method='min').astype(int)

final_df = results_df.merge(rank_df, on='Method').sort_values('Final_Rank')

# ===================================================================
# 저장
# ===================================================================

final_df.to_csv("DEG_Comparison_Final.csv", index=False)
rank_df.sort_values('Final_Rank').to_csv("DEG_Tool_Ranking.csv", index=False)

# ===================================================================
# 출력
# ===================================================================

print("\n" + "=" * 80)
print("ESSENTIAL METRICS")
print("=" * 80)
print(results_df[['Method', 'SC_DEGs', 'Overlap', 'F1',
                   'SCC', 'AUPRC']].to_string(index=False))

print("\n" + "=" * 80)
print("TOP-K OVERLAP & DIRECTION")
print("=" * 80)
topk_disp_cols = ['Method'] + [
    col for k in TOP_K_VALUES
    for col in [f'Top{k}_Overlap', f'Top{k}_Direction']
]
print(results_df[topk_disp_cols].to_string(index=False))

print("\n" + "=" * 80)
print("TOOL RANKING")
print("=" * 80)
rank_display = (rank_df.sort_values('Final_Rank')
                [['Method', 'Final_Rank', 'Mean_Rank'] +
                 [f'rank_{m}' for m in rank_metrics]]
                .rename(columns={f'rank_{m}': m for m in rank_metrics}))
print(rank_display.to_string(index=False))

print("\nTop 3:")
for _, row in rank_display.head(3).iterrows():
    print(f"  #{int(row['Final_Rank'])} {row['Method']}")

print("Bottom 3:")
for _, row in rank_display.tail(3).iterrows():
    print(f"  #{int(row['Final_Rank'])} {row['Method']}")

print("\n✅ Saved: DEG_Comparison_Final.csv")
print("✅ Saved: DEG_Tool_Ranking.csv")

# ===================================================================
# DEG Comparison Visualization (v3)
#
# Figure 1: Heatmap       — 전체 메트릭 overview
# Figure 2: Top-K Line    — Top-K Overlap 라인 플롯
# Figure 3: LFC Scatter   — 모든 method LFC 상관 (한 장)
# Figure 4a: Bar chart    — 종합 랭킹 (higher is better)
# Figure 4b: Radar chart  — 종합 랭킹
# ===================================================================

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr
import os

# ===================================================================
# ★ 여기서 Top-K 값을 자유롭게 지정하세요 ★
# (DEG_Comparison_Final.py 와 동일하게 맞춰주세요)
# ===================================================================

TOP_K_VALUES = [100, 200, 300, 400, 500]

# ===================================================================
# 경로 설정
# ===================================================================

DEG_DIR = f"{BASE_DIR}/Chu/DEG"
FIG_DIR = f"{BASE_DIR}/Chu/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# ===================================================================
# 데이터 로드
# ===================================================================

os.chdir(DEG_DIR)
final_df  = pd.read_csv("DEG_Comparison_Final.csv")
rank_df   = pd.read_csv("DEG_Tool_Ranking.csv")
bulk_deg  = pd.read_csv("BULK_DEG_H1_vs_DEC_results.csv")
sc_deg    = pd.read_csv("SC_DEG_All_Methods_DEC_vs_H1_Wilcoxon.csv")

# 정렬: Final_Rank 기준
final_df = final_df.sort_values('Final_Rank').reset_index(drop=True)
rank_df  = rank_df.sort_values('Final_Rank').reset_index(drop=True)

# TOP_K_VALUES CSV에 실제로 존재하는 K만 사용
available_k = [k for k in TOP_K_VALUES
               if f'Top{k}_Overlap' in final_df.columns]
if available_k != TOP_K_VALUES:
    print(f"⚠ TOP_K_VALUES 중 CSV에 없는 K 제외: "
          f"{set(TOP_K_VALUES) - set(available_k)}")
TOP_K_VALUES = available_k

print(f"Using Top-K values: {TOP_K_VALUES}")

# ===================================================================
# 공통 설정
# ===================================================================

METHOD_ORDER = final_df['Method'].tolist()
RAW_IDX      = METHOD_ORDER.index('RAW')

PALETTE = plt.cm.tab10.colors
method_colors = {}
imputed_idx = 0
for m in METHOD_ORDER:
    if m == 'RAW':
        method_colors[m] = 'black'
    else:
        method_colors[m] = PALETTE[imputed_idx % len(PALETTE)]
        imputed_idx += 1

# ===================================================================
# Figure 1: Heatmap
# ===================================================================

def plot_heatmap(final_df, method_order, raw_idx, k_values):
    df = final_df.set_index('Method').copy()
    
    # Top-K 전체 평균 컬럼
    df['TopK_Overlap_Mean']   = df[[f'Top{k}_Overlap'   for k in k_values]].mean(axis=1)
    df['TopK_Direction_Mean'] = df[[f'Top{k}_Direction' for k in k_values]].mean(axis=1)
    
    heatmap_metrics = ['F1', 'SCC', 'AUPRC',
                       'TopK_Overlap_Mean', 'TopK_Direction_Mean']
    k_label = f"K={','.join(str(k) for k in k_values)}"
    metric_labels = ['F1', 'SCC', 'AUPRC', 'Top-K\nOverlap', 'Top-K\nDirection']
    
    df_plot = df.loc[method_order, heatmap_metrics]
    df_norm = (df_plot - df_plot.min()) / (df_plot.max() - df_plot.min())
    
    fig, ax = plt.subplots(figsize=(9, 6))
    im = ax.imshow(df_norm.values, cmap='Blues', aspect='auto', vmin=0, vmax=1)
    
    for i in range(len(method_order)):
        for j in range(len(heatmap_metrics)):
            val      = df_plot.values[i, j]
            norm_val = df_norm.values[i, j]
            text_color = 'white' if norm_val > 0.6 else 'black'
            ax.text(j, i, f'{val:.3f}', ha='center', va='center',
                    fontsize=9, color=text_color)
            
    for j in range(len(heatmap_metrics)):
        ax.add_patch(plt.Rectangle(
            (j - 0.5, raw_idx - 0.5), 1, 1,
            fill=False, edgecolor='red', linewidth=2, zorder=3
        ))
        
    ax.set_xticks(range(len(heatmap_metrics)))
    ax.set_xticklabels(metric_labels, fontsize=10)
    ax.set_yticks(range(len(method_order)))
    ax.set_yticklabels(method_order, fontsize=10)
    ax.get_yticklabels()[raw_idx].set_fontweight('bold')
    ax.get_yticklabels()[raw_idx].set_color('red')
    
    cbar = plt.colorbar(im, ax=ax, label='Normalized Score', shrink=0.8)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    
    ax.set_title(
        f'DEG Benchmark: All Metrics Overview\n'
        f'(Top-K Overlap & Direction = mean across {k_label} | red = RAW)',
        fontsize=11, fontweight='bold', pad=15
    )
    plt.tight_layout()
    path = os.path.join(FIG_DIR, "Fig1_DEG_Heatmap.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {path}")


# ===================================================================
# Figure 2: Top-K Line Plot
# ===================================================================

def plot_topk_line(final_df, method_order, method_colors, k_values):
    k_cols = [f'Top{k}_Overlap' for k in k_values]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    for method in method_order:
        row    = final_df[final_df['Method'] == method].iloc[0]
        values = [row[c] for c in k_cols]
        ls     = '--' if method == 'RAW' else '-'
        zorder = 5 if method == 'RAW' else 2
        ax.plot(k_values, values,
                marker='o', label=method,
                color=method_colors[method],
                linewidth=1.5,
                linestyle=ls,
                zorder=zorder,
                markersize=5)
        
    ax.set_xlabel('K', fontsize=12)
    ax.set_ylabel('Overlap with Bulk Top-K', fontsize=12)
    ax.set_title(
        'Top-K Overlap: SC vs Bulk (Wald stat ranking)\n(dashed = RAW)',
        fontsize=12, fontweight='bold'
    )
    ax.set_xticks(k_values)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    path = os.path.join(FIG_DIR, "Fig2_TopK_Line.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {path}")


# ===================================================================
# Figure 3: LFC Scatter — 모든 tool, 한 장
# ===================================================================

def plot_lfc_scatter(bulk_deg, sc_deg, method_order):
    n_methods = len(method_order)
    n_cols    = 4
    n_rows    = int(np.ceil(n_methods / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(n_cols * 4, n_rows * 4))
    axes_flat = axes.flatten() if n_methods > 1 else [axes]
    
    for ax, method in zip(axes_flat, method_order):
        sc_m   = sc_deg[sc_deg['Method'] == method][['gene', 'log2FoldChange']].copy()
        merged = (bulk_deg[['gene', 'log2FoldChange', 'significant']]
                  .merge(sc_m, on='gene', suffixes=('_bulk', '_sc'))
                  .dropna())
        
        non_sig = merged[~merged['significant']]
        sig_up  = merged[merged['significant'] & (merged['log2FoldChange_bulk'] > 0)]
        sig_dn  = merged[merged['significant'] & (merged['log2FoldChange_bulk'] <= 0)]
        
        ax.scatter(non_sig['log2FoldChange_bulk'], non_sig['log2FoldChange_sc'],
                   alpha=0.15, s=3, color='gray', rasterized=True)
        ax.scatter(sig_up['log2FoldChange_bulk'], sig_up['log2FoldChange_sc'],
                   alpha=0.6, s=8, color='red',  label='UP DEG', rasterized=True)
        ax.scatter(sig_dn['log2FoldChange_bulk'], sig_dn['log2FoldChange_sc'],
                   alpha=0.6, s=8, color='blue', label='DOWN DEG', rasterized=True)
        
        sc_lim   = max(abs(merged['log2FoldChange_sc']).quantile(0.99) * 1.1, 0.5)
        bulk_lim = abs(merged['log2FoldChange_bulk']).max() * 1.05
        ax.set_ylim(-sc_lim, sc_lim)
        ax.plot([-bulk_lim, bulk_lim], [-bulk_lim, bulk_lim],
                'k--', linewidth=0.8, alpha=0.4)
        ax.axhline(0, color='gray', linewidth=0.5, alpha=0.4)
        ax.axvline(0, color='gray', linewidth=0.5, alpha=0.4)
        
        scc, _ = spearmanr(merged['log2FoldChange_bulk'], merged['log2FoldChange_sc'])
        ax.text(0.05, 0.95, f'SCC = {scc:.3f}',
                transform=ax.transAxes, fontsize=9,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        ax.set_title(method, fontsize=10, fontweight='bold')
        ax.set_xlabel('Bulk log2FC', fontsize=8)
        ax.set_ylabel('SC log2FC', fontsize=8)
        ax.tick_params(labelsize=7)
        
    for ax in axes_flat[n_methods:]:
        ax.set_visible(False)
        
    handles = [
        mpatches.Patch(color='red',  label='UP DEG (bulk sig)'),
        mpatches.Patch(color='blue', label='DOWN DEG (bulk sig)'),
        mpatches.Patch(color='gray', label='Non-DEG'),
    ]
    fig.legend(handles=handles, loc='lower right',
               bbox_to_anchor=(1.0, 0.0), fontsize=10, framealpha=0.9)
    fig.suptitle('LFC Correlation: SC vs Bulk RNA-seq (all methods)',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    
    path = os.path.join(FIG_DIR, "Fig3_LFC_Scatter.png")
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {path}")


# ===================================================================
# Figure 4a: Bar Chart — 종합 랭킹 (higher is better)
# ===================================================================

def plot_ranking_bar(rank_df, method_colors):
    df        = rank_df.sort_values('Final_Rank').copy()
    rank_cols = [c for c in rank_df.columns if c.startswith('rank_')]
    N         = len(df)
    
    df['HiB_Score'] = N + 1 - df['Mean_Rank']
    df['rank_std']  = df[rank_cols].std(axis=1)
    df = df.sort_values('HiB_Score', ascending=False).reset_index(drop=True)
    
    fig, ax = plt.subplots(figsize=(10, 5))
    bars = ax.bar(
        df['Method'], df['HiB_Score'],
        color=[method_colors[m] for m in df['Method']],
        edgecolor='black', linewidth=0.8,
        yerr=df['rank_std'], capsize=4,
        error_kw={'linewidth': 1.2}
    )
    
    raw_bar_idx = df['Method'].tolist().index('RAW')
    bars[raw_bar_idx].set_edgecolor('red')
    bars[raw_bar_idx].set_linewidth(2.5)
    
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(i, row['HiB_Score'] + row['rank_std'] + 0.05,
                f"#{int(row['Final_Rank'])}", ha='center', va='bottom', fontsize=9)
        
    ax.set_ylabel('Score (higher = better)', fontsize=11)
    ax.set_title(
        'DEG Benchmark: Tool Ranking (higher is better)\n'
        '(error bar = rank std across metrics | red border = RAW)',
        fontsize=12, fontweight='bold'
    )
    ax.set_ylim(0, df['HiB_Score'].max() + df['rank_std'].max() + 1.5)
    
    raw_patch     = mpatches.Patch(facecolor='black', edgecolor='red',
                                    linewidth=2, label='RAW (no imputation)')
    imputed_patch = mpatches.Patch(facecolor='steelblue', label='Imputed methods')
    ax.legend(handles=[raw_patch, imputed_patch], fontsize=9)
    plt.tight_layout()
    
    path = os.path.join(FIG_DIR, "Fig4a_Ranking_Bar.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {path}")


# ===================================================================
# Figure 4b: Radar Chart
# ===================================================================

def plot_ranking_radar(final_df, method_order, method_colors, k_values):
    df = final_df.set_index('Method').copy()
    
    df['TopK_Overlap']   = df[[f'Top{k}_Overlap'   for k in k_values]].mean(axis=1)
    df['TopK_Direction'] = df[[f'Top{k}_Direction' for k in k_values]].mean(axis=1)
    
    radar_metrics = ['F1', 'SCC', 'AUPRC', 'TopK_Overlap', 'TopK_Direction']
    radar_labels  = ['F1', 'SCC', 'AUPRC', 'Top-K\nOverlap', 'Top-K\nDirection']
    
    df_plot = df.loc[method_order, radar_metrics].copy()
    df_norm = (df_plot - df_plot.min()) / (df_plot.max() - df_plot.min())
    
    N      = len(radar_metrics)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]
    
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    
    for method in method_order:
        values = df_norm.loc[method, radar_metrics].tolist()
        values += values[:1]
        ls     = '--' if method == 'RAW' else '-'
        zorder = 5 if method == 'RAW' else 2
        alpha  = 0.25 if method == 'RAW' else 0.05
        
        ax.plot(angles, values,
                color=method_colors[method],
                linewidth=1.5,
                linestyle=ls,
                label=method,
                zorder=zorder)
        ax.fill(angles, values,
                color=method_colors[method], alpha=alpha)
        
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(radar_labels, fontsize=11)
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0.25', '0.50', '0.75', '1.00'], fontsize=8)
    ax.grid(True, alpha=0.3)
    k_label = f"K={','.join(str(k) for k in k_values)}"
    ax.set_title(
        f'DEG Benchmark: Method Profile\n'
        f'(Top-K mean across {k_label} | dashed = RAW)',
        fontsize=12, fontweight='bold', pad=20
    )
    ax.legend(bbox_to_anchor=(1.3, 1.1), loc='upper left', fontsize=9)
    plt.tight_layout()
    path = os.path.join(FIG_DIR, "Fig4b_Ranking_Radar.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {path}")


# ===================================================================
# 실행
# ===================================================================

print("=" * 60)
print("Generating DEG Benchmark Figures (v3)...")
print("=" * 60)

plot_heatmap(final_df, METHOD_ORDER, RAW_IDX, TOP_K_VALUES)
plot_topk_line(final_df, METHOD_ORDER, method_colors, TOP_K_VALUES)
plot_lfc_scatter(bulk_deg, sc_deg, METHOD_ORDER)
plot_ranking_bar(rank_df, method_colors)
plot_ranking_radar(final_df, METHOD_ORDER, method_colors, TOP_K_VALUES)

print("\n" + "=" * 60)
print(f"✅ All figures saved to: {FIG_DIR}")
print("  Fig1_DEG_Heatmap.png")
print("  Fig2_TopK_Line.png")
print("  Fig3_LFC_Scatter.png")
print("  Fig4a_Ranking_Bar.png")
print("  Fig4b_Ranking_Radar.png")
print("=" * 60)

# ===================================================================
# DEG Comparison — Individual Metric Figures
#
# FigA: F1 / SCC / AUPRC — 3-col barplot (1장)
# FigB: Top-K Overlap / Top-K Direction — 2-col line plot (1장)
# ===================================================================

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# ===================================================================
# ★ DEG_Comparison_Final.py 와 동일하게 맞추세요 ★
# ===================================================================

TOP_K_VALUES = [100, 200, 300, 400, 500]

# ===================================================================
# 경로 설정
# ===================================================================

DEG_DIR = f"{BASE_DIR}/Chu/DEG"
FIG_DIR = f"{BASE_DIR}/Chu/figures"
os.makedirs(FIG_DIR, exist_ok=True)

os.chdir(DEG_DIR)
final_df = pd.read_csv("DEG_Comparison_Final.csv")
final_df = final_df.sort_values('Final_Rank').reset_index(drop=True)

# TOP_K_VALUES 중 CSV에 실제로 있는 K만 사용
TOP_K_VALUES = [k for k in TOP_K_VALUES
                if f'Top{k}_Overlap' in final_df.columns]
print(f"Using Top-K values: {TOP_K_VALUES}")

# ===================================================================
# 색상 설정
# ===================================================================

METHOD_ORDER = final_df['Method'].tolist()

# RAW를 맨 앞으로
if 'RAW' in METHOD_ORDER:
    METHOD_ORDER = ['RAW'] + [m for m in METHOD_ORDER if m != 'RAW']

PALETTE = plt.cm.tab10.colors
method_colors = {}
imp_idx = 0
for m in METHOD_ORDER:
    if m == 'RAW':
        method_colors[m] = 'black'
    else:
        method_colors[m] = PALETTE[imp_idx % len(PALETTE)]
        imp_idx += 1

colors_ordered = [method_colors[m] for m in METHOD_ORDER]

# ===================================================================
# FigA: F1 / SCC / AUPRC — 3-col barplot
# ===================================================================

def plot_scalar_metrics(final_df, method_order, colors_ordered):
    metrics = [
        ('F1',    'F1 Score',  'F1'),
        ('SCC',   'Spearman Correlation (SCC)', 'SCC (LFC)'),
        ('AUPRC', 'AUPRC',     'AUPRC'),
    ]
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
    
    for ax, (col, ylabel, title) in zip(axes, metrics):
        vals  = final_df.set_index('Method').loc[method_order, col].values
        bars  = ax.bar(method_order, vals,
                       color=colors_ordered,
                       edgecolor='black', linewidth=0.7)
        
        # RAW 기준선 (가로 점선)
        raw_val = vals[method_order.index('RAW')]
        ax.axhline(raw_val, color='black', linewidth=1.2,
                   linestyle='--', alpha=0.6, zorder=3, label='RAW baseline')
        
        # 값 텍스트
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    v + max(vals) * 0.01,
                    f'{v:.3f}', ha='center', va='bottom', fontsize=7.5)
            
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_xticks(range(len(method_order)))
        ax.set_xticklabels(method_order, rotation=35, ha='right', fontsize=9)
        ax.set_ylim(0, max(vals) * 1.15)
        ax.grid(axis='y', alpha=0.3)
    
    fig.suptitle('DEG Benchmark: F1 / SCC / AUPRC per Method\n(red = RAW)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    path = os.path.join(FIG_DIR, "FigA_DEG_F1_SCC_AUPRC.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {path}")


# ===================================================================
# FigB: Top-K Overlap / Top-K Direction — 2-col line plot
# ===================================================================

def plot_topk_lines(final_df, method_order, method_colors, k_values):
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    
    panels = [
        (axes[0], 'Overlap',   'Overlap with Bulk Top-K',        'Top-K Overlap'),
        (axes[1], 'Direction', 'Direction Concordance (0–1)',     'Top-K Direction'),
    ]
    
    for ax, suffix, ylabel, title in panels:
        k_cols = [f'Top{k}_{suffix}' for k in k_values]
        
        # RAW 기준값 (각 K별) → 가로 점선용
        raw_row = final_df[final_df['Method'] == 'RAW'].iloc[0]
        raw_vals = [raw_row[c] for c in k_cols]
        
        for method in method_order:
            row    = final_df[final_df['Method'] == method].iloc[0]
            values = [row[c] for c in k_cols]
            ls     = '--' if method == 'RAW' else '-'
            zorder = 5 if method == 'RAW' else 2
            
            ax.plot(k_values, values,
                    marker='o',
                    label=method,
                    color=method_colors[method],
                    linewidth=1.5,
                    linestyle=ls,
                    zorder=zorder,
                    markersize=5)
            
        ax.set_xlabel('K', fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xticks(k_values)
        ax.grid(True, alpha=0.3)
        
        if suffix == 'Direction':
            # y축: 데이터 범위에 맞게 최적화
            all_vals = final_df[k_cols].values.flatten()
            all_vals = all_vals[~np.isnan(all_vals)]
            margin = (all_vals.max() - all_vals.min()) * 0.5
            ymin = max(0, all_vals.min() - margin)
            ymax = 1.0005   # 1.0 클리핑 제거
            ax.set_ylim(ymin, ymax)
            
    # 공통 범례 (오른쪽 바깥)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               bbox_to_anchor=(1.01, 0.5), loc='center left', fontsize=9)
    
    k_label = f"K = {', '.join(str(k) for k in k_values)}"
    fig.suptitle(f'DEG Benchmark: Top-K Overlap & Direction\n'
                 f'({k_label} | dashed = RAW)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    path = os.path.join(FIG_DIR, "FigB_DEG_TopK_Overlap_Direction.png")
    plt.savefig(path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Saved: {path}")


# ===================================================================
# 실행
# ===================================================================

print("=" * 60)
print("Generating individual metric figures...")
print("=" * 60)

plot_scalar_metrics(final_df, METHOD_ORDER, colors_ordered)
plot_topk_lines(final_df, METHOD_ORDER, method_colors, TOP_K_VALUES)

print("\n✅ Done.")
print("  FigA_DEG_F1_SCC_AUPRC.png")
print("  FigB_DEG_TopK_Overlap_Direction.png")

import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
import os
import pandas as pd
import math

# =====================
# 1. 설정
# =====================
save_dir = f"{BASE_DIR}/Chu/figures"
os.makedirs(save_dir, exist_ok=True)

filename = "LFC_Scatter_Bulk_vs_SC_MAST_Top1000_AllMethods.png"
save_path = os.path.join(save_dir, filename)

methods = list(sc_deg.keys())
n_methods = len(methods)

ncols = 4
nrows = math.ceil(n_methods / ncols)

# =====================
# 2. Figure
# =====================
fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(4.5 * ncols, 4 * nrows),
    squeeze=False
)
colors = {'UP': '#E41A1C', 'DOWN': '#377EB8'}

# =====================
# 3. Loop per method
# =====================
for ax, method in zip(axes.flatten(), methods):
    sc_df = sc_deg[method]
    # --- merge ---
    df = bulk_deg.merge(
        sc_df,
        left_index=True,
        right_index=True
    )
    # --- bulk 기준 정렬 ---
    sort_col = "padj"
    # --- Top 500 UP / DOWN (bulk 기준) ---
    up = df[df['diff_status_x'] == 'UP'].sort_values(sort_col).head(500)
    down = df[df['diff_status_x'] == 'DOWN'].sort_values(sort_col).head(500)
    top_df = pd.concat([up, down])
    # --- 컬럼 명시 ---
    x_col = 'log2FoldChange'   # bulk
    y_col = 'avg_log2FC'       # single-cell
    
    top_df = top_df.dropna(subset=[x_col, y_col])
    if len(top_df) == 0:
        ax.set_title(f"{method} (no data)")
        ax.axis("off")
        continue
    # --- axis limit ---
    max_x = np.max(np.abs(top_df[x_col]))
    max_y = np.max(np.abs(top_df[y_col]))
    limit_x = max_x * 1.1
    limit_y = max_y * 1.1
    # --- correlation ---
    p_r, _ = stats.pearsonr(top_df[x_col], top_df[y_col])
    s_r, _ = stats.spearmanr(top_df[x_col], top_df[y_col])
    # --- scatter ---
    sns.scatterplot(
        data=top_df,
        x=x_col,
        y=y_col,
        hue='diff_status_x',   # ← 여기
        palette=colors,
        alpha=0.6,
        s=18,
        edgecolor=None,
        legend=False,
        ax=ax
    )
    sns.regplot(
        data=top_df,
        x=x_col,
        y=y_col,
        scatter=False,
        color='black',
        line_kws={"linestyle": "--", "linewidth": 1.2},
        ax=ax
    )
    ax.axhline(0, color='grey', linewidth=0.7)
    ax.axvline(0, color='grey', linewidth=0.7)
    ax.set_xlim(-limit_x, limit_x)
    ax.set_ylim(-limit_y, limit_y)
    ax.set_title(
    f"{method_display[method]}\nPearson r={p_r:.2f}, Spearman r={s_r:.2f}",
    fontsize=11
    )
    ax.set_xlabel("Bulk LFC")
    ax.set_ylabel("SC LFC")

# =====================
# 4. 빈 subplot 제거
# =====================
for ax in axes.flatten()[len(methods):]:
    ax.axis("off")

plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches="tight")
plt.close()
print(f"Saved to: {save_path}")


# === GSEA 수행 및 결과 저장 ===
import numpy as np
import pandas as pd
import gseapy as gp
import os

# 저장 경로 설정
GSEA_BASE_DIR = os.path.join(BASE_DIR, "GSEA_results")
if not os.path.exists(GSEA_BASE_DIR):
    os.makedirs(GSEA_BASE_DIR)

# 사용할 Gene Set (GO Biological Process 등)
# 라이브러리에 따라 'GO_Biological_Process_2023' 등을 선택하세요.
gene_set = 'GO_Biological_Process_2023'

# --- [STEP 1] Bulk GSEA 결과 (Ground Truth) 준비 ---
print("Running GSEA for Bulk data...")
# Bulk Sign(P)로 랭킹 생성
bulk_rnk = bulk_deg[['sign_p']].sort_values(by='sign_p', ascending=False).dropna()
bulk_gsea = gp.prerank(rnk=bulk_rnk, gene_sets=gene_set, permutation_num=1000, outdir=None).res2d
bulk_gsea.to_csv(os.path.join(GSEA_BASE_DIR, "GSEA_Bulk_Reference.csv"))

# --- [STEP 2] Method별 Single-cell GSEA 실행 ---
sc_gsea_results = {}

for method, df in sc_deg.items():
    display_name = method_display.get(method, method)
    print(f"Running GSEA for {display_name}...")
    
    # 1. SC Sign(P) 계산 및 랭킹 생성
    temp_df = df.copy()
    temp_df['sc_sign_p'] = -np.log10(temp_df['p_val'] + 1e-300) * np.sign(temp_df['avg_log2FC'])
    
    # 유전자 이름(index)과 Sign(P)만 추출하여 정렬
    sc_rnk = temp_df[['sc_sign_p']].sort_values(by='sc_sign_p', ascending=False).dropna()
    
    # 2. GSEA Prerank 실행
    try:
        pre_res = gp.prerank(rnk=sc_rnk, 
                             gene_sets=gene_set, 
                             permutation_num=1000, 
                             outdir=None, 
                             seed=42).res2d
        
        # 결과 저장
        save_path = os.path.join(GSEA_BASE_DIR, f"GSEA_{method}_results.csv")
        pre_res.to_csv(save_path)
        sc_gsea_results[method] = pre_res
        
    except Exception as e:
        print(f"Error running GSEA for {display_name}: {e}")

print(">>> All GSEA runs completed.")

# ==========Python: scRNA-seq==========
# 필요한 라이브러리 임포트
import os
import scanpy as sc
import numpy as np
import pandas as pd

os.chdir(f"{BASE_DIR}/Chu_time/Raw")
raw = sc.read_h5ad("chu.h5ad")
cell_idx, gene_idx = raw.obs.index.tolist(), raw.var.index.tolist()
sc.pp.normalize_total(raw, target_sum=1e4)
sc.pp.log1p(raw)
sc.pp.highly_variable_genes(raw, n_top_genes=2000)
raw.write_h5ad("./chu_norm.h5ad")# 최종 저장 및 불러오기

# imputed datasets
os.chdir(f"{BASE_DIR}/Chu_time/imputed")

# === 1. DCA : Done ===
#dca_df = pd.read_csv("./DCA/mean.tsv", sep="\t", header=0, index_col=0).T
#dca = sc.AnnData(dca_df)
#dca.obs, dca.var = raw.obs.copy(), raw.var.loc[dca.var_names].copy()
dca = sc.read_h5ad("./DCA/dca.h5ad")
sc.pp.normalize_total(dca, target_sum=1e4)
sc.pp.log1p(dca)
sc.pp.highly_variable_genes(dca, n_top_genes=2000)
dca.write_h5ad("./DCA/dca.h5ad")# 최종 저장 및 불러오기

# === 2. DrImpute: Done ===
#drimpute_df = pd.read_csv("./DrImpute/chu_DrImpute.csv", sep=",", index_col=0).T
#drimpute = sc.AnnData(drimpute_df)
#drimpute.obs, drimpute.var = raw.obs.copy(), raw.var.copy()
#drimpute.write_h5ad("./DrImpute/drimpute.h5ad")
drimpute = sc.read_h5ad("./DrImpute/drimpute.h5ad")
sc.pp.highly_variable_genes(drimpute, n_top_genes=2000)

# === 3. MAGIC: Done ===
magic = sc.read_h5ad("./MAGIC/chu_MAGIC.h5ad")
sc.pp.highly_variable_genes(magic, n_top_genes=2000)
magic.write_h5ad("./MAGIC/chu_MAGIC.h5ad")

# === 4. SAVER: Done ===
#saver_df = pd.read_csv("./SAVER/chu_saver_imputed.csv", sep=",", index_col=0).T
#saver = sc.AnnData(saver_df)
#saver.obs, saver.var = raw.obs.copy(), raw.var.copy()
saver = sc.read_h5ad("./SAVER/saver.h5ad")
sc.pp.normalize_total(saver, target_sum=1e4)
sc.pp.log1p(saver)
sc.pp.highly_variable_genes(saver, n_top_genes=2000)
#saver.write_h5ad("./SAVER/saver.h5ad")

# === 5. scIDPMs  ===
#scidpms_df = pd.read_csv("./scIDPMs/imputed.csv", sep=",", index_col=None, header=None)
#scidpms = sc.AnnData(scidpms_df)
#scidpms.obs, scidpms.var = raw.obs.copy(), raw.var.copy()
scidpms = sc.read_h5ad("./scIDPMs/scidpms.h5ad")
sc.pp.normalize_total(scidpms, target_sum=1e4)
sc.pp.log1p(scidpms)
sc.pp.highly_variable_genes(scidpms, n_top_genes=2000)
#scidpms.write_h5ad("./scIDPMs/scidpms.h5ad")

# === 6. scIGANs ===
#scigans_df = pd.read_csv("./scIGANs/scIGANs_chu.txt", sep="\t", index_col=0).T
#scigans = sc.AnnData(scigans_df)
#scigans.obs, scigans.var = raw.obs.copy(), raw.var.copy()
scigans = sc.read_h5ad("./scIGANs/scigans.h5ad")
sc.pp.normalize_total(scigans, target_sum=1e4)
sc.pp.log1p(scigans)
sc.pp.highly_variable_genes(scigans, n_top_genes=2000)
#scigans.write_h5ad("./scIGANs/scigans.h5ad")

# === 8. scMultiGAN: Done ===
#scmultigan_df = pd.read_csv("./scMultiGAN/Restored/chu_time_scMultiGAN.tsv", sep="\t", header=None, index_col=None).T
#scmultigan = sc.AnnData(scmultigan_df)
#scmultigan.obs, scmultigan.var = raw.obs.copy(), raw.var.copy()
scmultigan = sc.read_h5ad("./scMultiGAN/scmultigan.h5ad")
sc.pp.normalize_total(scmultigan, target_sum=1e4)
sc.pp.log1p(scmultigan)
sc.pp.highly_variable_genes(scmultigan, n_top_genes=2000)
#scmultigan.write_h5ad("./scMultiGAN/scmultigan.h5ad")

# === 9. scSTD: Done ===
#scstd_df = pd.read_csv("./scSTD/Result/chu_imputed.txt", sep=",", index_col=None, header=None)
#scstd = sc.AnnData(scstd_df)
#scstd.obs, scstd.var = raw.obs.copy(), raw.var.copy()
scstd = sc.read_h5ad("./scSTD/scstd.h5ad")
sc.pp.normalize_total(scstd, target_sum=1e4)
sc.pp.log1p(scstd)
sc.pp.highly_variable_genes(scstd, n_top_genes=2000)
#scstd.write_h5ad("./scSTD/scstd.h5ad")

# === 10. scVI: Done ===
scvi = sc.read_h5ad("./scVI/scvi.h5ad")
scvi.obs_names, scvi.var_names = cell_idx, gene_idx
print("scvi shape:", scvi.shape)
sc.pp.highly_variable_genes(scvi, n_top_genes=2000)
scvi.write_h5ad("./scVI/scvi.h5ad")

# === 11. ALRA: Done ===
alra_df = pd.read_csv("./ALRA/chu_ALRA_imputed.csv", index_col=0)
alra = sc.AnnData(alra_df)
alra.obs, alra.var = raw.obs.copy(), raw.var.copy()
alra.obs_names, alra.var_names = cell_idx, gene_idx
alra.write_h5ad("./ALRA/alra.h5ad")
alra = sc.read_h5ad("./ALRA/alra.h5ad")
sc.pp.highly_variable_genes(alra, n_top_genes=2000)

# === 12. afMF: Done ===
#afmf_df = pd.read_csv("./afMF/chu_afMF_imputed.csv", index_col=0)
#afmf = sc.AnnData(afmf_df)
#afmf.obs, afmf.var = raw.obs.copy(), raw.var.copy()
#afmf.obs_names, afmf.var_names = cell_idx, gene_idx
#afmf.write_h5ad("./afMF/afmf.h5ad")
afmf = sc.read_h5ad("./afMF/afmf.h5ad")
sc.pp.highly_variable_genes(afmf, n_top_genes=2000)

library(monocle3)
library(anndata)
library(ggplot2)
library(patchwork) # ncol 설정을 위해 필수

# 1. 경로 설정
input_path <- f"{BASE_DIR}/Chu_time/imputed/"
figure_path <- f"{BASE_DIR}/Chu_time/figures/"
if(!dir.exists(figure_path)) dir.create(figure_path)

# 2. 분석할 파일 리스트 (필요에 따라 파일명 수정)
file_list <- c("raw" = "../Raw/chu_norm.h5ad", 
               "DCA" = "DCA/dca.h5ad", 
               "DrImpute" = "DrImpute/drimpute.h5ad",
               "MAGIC" = "MAGIC/chu_MAGIC.h5ad",
               "SAVER" = "SAVER/saver.h5ad",
               "scIDPMs" = "scIDPMs/scidpms.h5ad",
               "scIGANs" = "scIGANs/scigans.h5ad",
               "scMultiGAN" = "scMultiGAN/scmultigan.h5ad",
               "scSTD" = "scSTD/scstd.h5ad",
               "scVI" = "scVI/scvi.h5ad",
               "ALRA" = "ALRA/alra.h5ad",
               "afMF" = "afMF/afmf.h5ad"
               )

plot_list <- list()

# 3. 루프 시작
for (name in names(file_list)) {
  message(paste0("Processing: ", name))
  
  # 데이터 로드
  adata <- read_h5ad(paste0(input_path, file_list[[name]]))
  
  # CDS 생성 및 이름 정제
  exp_mat <- t(adata$X)
  rownames(exp_mat) <- as.character(rownames(adata$var))
  colnames(exp_mat) <- as.character(rownames(adata$obs))
  
  cds <- new_cell_data_set(as(exp_mat, "dgCMatrix"),
                           cell_metadata = as.data.frame(adata$obs),
                           gene_metadata = data.frame(gene_short_name = rownames(exp_mat), row.names = rownames(exp_mat)))
  
  # 전처리 (기존 로직 보존)
 cds <- preprocess_cds(
  cds,
  num_dim = 50,
  norm_method = "none",
  method = "PCA",
  use_genes = rownames(adata$var)[adata$var$highly_variable]
)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  
  # Trajectory 학습 (기존 단순화 설정 적용)
  cds <- cluster_cells(cds, resolution=1e-3) 
  cds <- learn_graph(cds, use_partition = FALSE, 
                     learn_graph_control = list(minimal_branch_len = 10))
  cds <- order_cells(cds, root_cells = colnames(cds)[cds@colData$cell_type == "00h"])
  # 시각화 (요청하신 라벨 제거 및 사이즈 적용)
  p <- plot_cells(cds, color_cells_by = "cell_type", 
                  label_cell_groups = FALSE,
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  graph_label_size = 0,      
                  cell_size = 2.5,           # ncol=4이므로 조금 줄였습니다
                  trajectory_graph_segment_size = 1.2) + 
       theme_classic() +
       ggtitle(name)
       #theme(legend.position = "none") # 개별 그림의 범례는 일단 제거 (깔끔하게 보이기 위함)

  plot_list[[name]] <- p
}

# 4. 결합 및 저장 (ncol = 4)
final_plot <- wrap_plots(plot_list, ncol = 4) + plot_layout(guides = 'collect')

ggsave(filename = paste0(figure_path, "All_Models_Trajectory.png"), 
       plot = final_plot, width = 24, height = 15) # ncol=4에 맞게 가로세로 비율 조정

print(paste0("All plots saved at: ", figure_path))

library(anndata)
library(ggplot2)
library(reshape2)
library(patchwork)

# =========================
# 1. 경로 및 설정
# =========================
input_path  <- f"{BASE_DIR}/Chu_time/imputed/"
figure_path <- f"{BASE_DIR}/Chu_time/figures/"
if (!dir.exists(figure_path)) dir.create(figure_path)

# 분석할 파일 리스트
# 2. 분석할 파일 리스트 (필요에 따라 파일명 수정)
file_list <- c("raw" = "../Raw/chu_norm.h5ad", 
               "DCA" = "DCA/dca.h5ad", 
               "DrImpute" = "DrImpute/drimpute.h5ad",
               "MAGIC" = "MAGIC/chu_MAGIC.h5ad",
               "SAVER" = "SAVER/saver.h5ad",
               "scIDPMs" = "scIDPMs/scidpms.h5ad",
               "scIGANs" = "scIGANs/scigans.h5ad",
               "scMultiGAN" = "scMultiGAN/scmultigan.h5ad",
               "scSTD" = "scSTD/scstd.h5ad",
               "scVI" = "scVI/scvi.h5ad",
               "ALRA" = "ALRA/alra.h5ad",
               "afMF" = "afMF/afmf.h5ad"
               )

all_pt_data <- data.frame()

# =========================
# 2. 분석 루프
# =========================
for (name in names(file_list)) {
  message(paste0("Calculating Pseudotime for: ", name))
  
  # 데이터 로드
  adata <- read_h5ad(paste0(input_path, file_list[[name]]))
  exp_mat <- t(adata$X)
  
  cds <- new_cell_data_set(
    as(exp_mat, "dgCMatrix"),
    cell_metadata = as.data.frame(adata$obs),
    gene_metadata = data.frame(gene_short_name = rownames(exp_mat), row.names = rownames(exp_mat))
  )
  
  # 의사시간 계산을 위한 최소 파이프라인 (기존의 '나은' 로직 적용)
  cds <- preprocess_cds(cds, num_dim = 50, norm_method = "none")
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds, resolution = 1e-3)
  cds <- learn_graph(cds, use_partition = FALSE, 
                     learn_graph_control = list(minimal_branch_len = 10))
  
  # 뿌리 세포 설정 (00h)
  cds <- order_cells(cds, root_cells = colnames(cds)[cds@colData$cell_type == "00h"])
  
  # 의사시간 및 메타데이터 추출
  pt_df <- data.frame(
    pseudotime = pseudotime(cds),
    cell_type = colData(cds)$cell_type,
    model = name
  )
  
  # 무한대 값(연결 안된 세포) 제외
  pt_df <- pt_df %>% filter(is.finite(pseudotime))
  
  all_pt_data <- rbind(all_pt_data, pt_df)
}

# 타임포인트 순서 정렬
all_pt_data$cell_type <- factor(all_pt_data$cell_type, 
                                 levels = c("00h", "12h", "24h", "36h", "72h", "96h"))

# =========================
# 3. Boxplot 생성 및 PNG 저장
# =========================
message("Saving Pseudotime Boxplot...")

p <- ggplot(all_pt_data, aes(x = cell_type, y = pseudotime, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~model, scales = "free_y", ncol = 3) + # 모델별로 박스플롯 배치
  theme_bw() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Pseudotime Distribution by Actual Time Points",
       subtitle = "Check the separation between 72h and 96h",
       x = "Actual Time (Cell Type)",
       y = "Monocle3 Pseudotime") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  )

# PNG 저장
ggsave(
  filename = paste0(figure_path, "Pseudotime_Distribution_Boxplot.png"),
  plot = p,
  width = 18,
  height = 10,
  dpi = 300
)

print(paste0("Plot saved at: ", figure_path, "Pseudotime_Distribution_Boxplot.png"))
