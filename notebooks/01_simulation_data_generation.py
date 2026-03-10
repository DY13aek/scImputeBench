# Set the base directory for datasets
# Change this path to point to your local data directory
BASE_DIR = "../data"

# ─────────────────────────────────────────────────────────────
# Headless 설정 (X11 불필요)
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(splatter); library(Matrix)
  library(SingleCellExperiment); library(scater)
  library(scran); library(ggplot2)
  library(zellkonverter)  # H5AD 저장에 필요
})

options(bitmapType = "cairo")

# 저장 디렉토리: Raw (요청대로 변경)
out_dir <- "../data/splatter/Raw"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────
# H5AD 저장 (raw counts→X, obs에 cell_type) 유틸
# ─────────────────────────────────────────────────────────────
save_h5ad_raw <- function(sce, label, dir = out_dir) {
  fn <- file.path(dir, paste0(label, ".h5ad"))
  if (file.exists(fn)) unlink(fn)

  # obs에 cell_type 추가 (Group 복사)
  SummarizedExperiment::colData(sce)$cell_type <-
    as.character(SummarizedExperiment::colData(sce)$Group)

  # 직렬화 안전화
  sce2 <- sce
  S4Vectors::metadata(sce2) <- list()
  cd <- S4Vectors::DataFrame(
    Group     = as.character(SummarizedExperiment::colData(sce)$Group),
    cell_type = as.character(SummarizedExperiment::colData(sce)$cell_type),
    row.names = colnames(sce)
  )
  SummarizedExperiment::colData(sce2) <- cd
  rd <- S4Vectors::DataFrame(row.names = rownames(sce))
  SummarizedExperiment::rowData(sce2) <- rd

  # altExp 제거 (안전)
  if (methods::is(sce2, "SingleCellExperiment")) {
    ane <- SingleCellExperiment::altExpNames(sce2)
    if (length(ane) > 0) {
      for (nm in ane) {
        SingleCellExperiment::altExp(sce2, nm) <- NULL
      }
    }
  }

  # raw counts를 X로 저장
  zellkonverter::writeH5AD(sce2, file = fn, X_name = "counts")
  message("💾 saved: ", fn)
}

# ─────────────────────────────────────────────────────────────
# 공통 설정
# ─────────────────────────────────────────────────────────────
set.seed(101)
n_genes <- 10000; n_cells <- 5000

params <- newSplatParams(
  nGenes       = n_genes,  # default 10000
  batchCells   = n_cells,   # default 100
  group.prob   = c(0.30, 0.25, 0.15, 0.10, 0.08, 0.06, 0.04, 0.02),

  #mean.shape   = 0.6, # default 0.6
  #mean.rate    = 0.3, # default 0.3

  #lib.loc      = 11, # default 11
  #lib.scale    = 0.2, # default 0.2

  #bcv.common   = 0.1, # default 0.1
  #de.facLoc    = 0.1, # default 0.1
  #de.prob      = 0.2, # default 0.2
  #out.prob     = 0.05, # default 0.05

  dropout.type = "experiment"
)

dropout_cases <- c(sim_low = 1, sim_mod = 2, sim_high = 3)
zero_rate <- function(mat) 1 - Matrix::nnzero(mat) / (nrow(mat) * ncol(mat))

# ─────────────────────────────────────────────────────────────
# GT 시뮬레이션 (dropout 없음)
# ─────────────────────────────────────────────────────────────
message("=== [GT] (dropout = none) ===")
sim_gt <- splatSimulate(params, method = "groups", dropout.type = "none", verbose = FALSE)
X_gt   <- counts(sim_gt)
zr_gt  <- zero_rate(X_gt)
message(sprintf("GT Zero-Rate (biological only): %.4f (%.2f%%)", zr_gt, 100*zr_gt))

# ─────────────────────────────────────────────────────────────
# 카운트 행렬 통계 유틸
# ─────────────────────────────────────────────────────────────
stats_from_counts <- function(X) {
  nz_by_cell  <- Matrix::colSums(X > 0)
  umi_by_cell <- Matrix::colSums(X)
  total_genes_detected <- sum(Matrix::rowSums(X > 0) > 0)

  list(
    med_genes_per_cell     = as.numeric(stats::median(nz_by_cell)),
    med_umi_per_cell       = as.numeric(stats::median(umi_by_cell)),
    total_genes_detected   = as.integer(total_genes_detected)
  )
}

# ─────────────────────────────────────────────────────────────
# 각 케이스 실행 (UMAP 없음, 시각화 없음)
# ─────────────────────────────────────────────────────────────
run_one <- function(mid, label, zr_gt, params) {
  message(sprintf("\n=== [%s] dropout.mid = %g ===", label, mid))
  params_case <- setParams(params, update = list(dropout.mid = mid))
  sim <- splatSimulate(params_case, method = "groups", verbose = FALSE)

  X        <- counts(sim)
  zr_total <- zero_rate(X)
  zr_tech  <- max(0, zr_total - zr_gt)

  s <- stats_from_counts(X)

  list(
    label                 = label,
    dropout.mid           = mid,
    zero_total_pct        = 100 * zr_total,
    zero_tech_approx_pct  = 100 * zr_tech,
    n_genes               = nrow(X),
    n_cells               = ncol(X),
    med_genes_per_cell    = s$med_genes_per_cell,
    med_umi_per_cell      = s$med_umi_per_cell,
    total_genes_detected  = s$total_genes_detected,
    sce                   = sim
  )
}

results <- lapply(names(dropout_cases), function(lbl) {
  run_one(mid = dropout_cases[[lbl]], label = lbl, zr_gt = zr_gt, params = params)
})
names(results) <- names(dropout_cases)

# ─────────────────────────────────────────────────────────────
# 요약 CSV 저장 (GT 포함, 새 지표 포함)
# ─────────────────────────────────────────────────────────────
gt_stats <- stats_from_counts(X_gt)
gt_entry <- list(
  label                 = "GT",
  dropout.mid           = NA_real_,
  zero_total_pct        = 100 * zr_gt,
  zero_tech_approx_pct  = 0,
  n_genes               = nrow(X_gt),
  n_cells               = ncol(X_gt),
  med_genes_per_cell    = gt_stats$med_genes_per_cell,
  med_umi_per_cell      = gt_stats$med_umi_per_cell,
  total_genes_detected  = gt_stats$total_genes_detected,
  sce                   = sim_gt
)

dropout_list <- lapply(names(dropout_cases), function(lbl) {
  run_one(mid = dropout_cases[[lbl]], label = lbl, zr_gt = zr_gt, params = params)
})
names(dropout_list) <- names(dropout_cases)

# GT를 맨 위로
results <- c(list(GT = gt_entry), dropout_list)


summary_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    label                   = x$label,
    dropout.mid             = x$dropout.mid,
    zero_total_pct          = x$zero_total_pct,
    zero_tech_approx_pct    = x$zero_tech_approx_pct,
    n_genes                 = x$n_genes,
    n_cells                 = x$n_cells,
    median_genes_per_cell   = x$med_genes_per_cell,
    median_UMI_per_cell     = x$med_umi_per_cell,
    total_genes_detected    = x$total_genes_detected,
    stringsAsFactors        = FALSE
  )
}))

csv_fn <- file.path(out_dir, "summary_dropout_cases.csv")
write.csv(summary_df, csv_fn, row.names = FALSE)
message("\n--- 완료 ---")
print(summary_df)
message("📄 saved: ", csv_fn)

# ─────────────────────────────────────────────────────────────
# H5AD 저장: GT / sim_low / sim_mod / sim_high (raw counts를 X로)
# ─────────────────────────────────────────────────────────────
save_h5ad_raw(sim_gt,                        "GT",  out_dir)
save_h5ad_raw(results[["sim_low"]]$sce,      "sim1", out_dir)
save_h5ad_raw(results[["sim_mod"]]$sce,      "sim2", out_dir)
save_h5ad_raw(results[["sim_high"]]$sce,     "sim3", out_dir)




# ============================
# UMAP (TkAgg, no saving)
# ============================
import os
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("TkAgg")  # 이 서버는 TkAgg OK
import matplotlib.pyplot as plt

# ----- 경로 -----
in_dir = f"{BASE_DIR}/splatter/Raw"

# ----- 라벨 통일: 'gt' (우선순위: cell_type > Group) -----
def ensure_gt_label(ad):
    if "gt" in ad.obs:
        ad.obs["gt"] = ad.obs["gt"].astype("category")
        return ad
    if "cell_type" in ad.obs:
        ad.obs["gt"] = ad.obs["cell_type"].astype("category")
    elif "Group" in ad.obs:
        ad.obs["gt"] = ad.obs["Group"].astype("category")
    else:
        ad.obs["gt"] = "unknown"
    ad.obs["gt"] = ad.obs["gt"].astype("category")
    return ad

# ----- 표준 Scanpy 파이프라인 -----
def scanpy_standard_pipeline(ad, n_top_genes=3000, n_pcs=50, n_neighbors=15,
                             min_dist=0.5, seed=42):
    sc.settings.seed = seed
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=n_top_genes)
    ad = ad[:, ad.var["highly_variable"]].copy()
    sc.pp.scale(ad, max_value=10)
    n_comps = int(min(n_pcs, ad.n_vars, max(2, ad.n_obs - 1)))
    sc.tl.pca(ad, n_comps=n_comps)
    sc.pp.neighbors(ad, n_neighbors=n_neighbors, n_pcs=n_comps, metric="euclidean")
    sc.tl.umap(ad, min_dist=min_dist)
    return ad

# ----- 로드 -----
files = {
    "GT":   "GT.h5ad",
    "sim1": "sim1.h5ad",
    "sim2": "sim2.h5ad",
    "sim3": "sim3.h5ad",
}

ads = {}
for key, fname in files.items():
    ad = sc.read_h5ad(os.path.join(in_dir, fname))  # X = counts
    ad = ensure_gt_label(ad)
    ads[key] = scanpy_standard_pipeline(ad)

# ----- 2×2 서브플롯으로 바로 띄우기 -----
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
axes = axes.flatten()
titles = ["GT", "sim1", "sim2", "sim3"]

for ax, name in zip(axes, titles):
    sc.pl.umap(
        ads[name],
        color="gt",
        title=f"{name}",
        legend_loc="on data",
        legend_fontsize=8,
        ax=ax,
        show=False,   # 전체 그린 뒤 한 번에 show()
    )

plt.tight_layout()
plt.show()



