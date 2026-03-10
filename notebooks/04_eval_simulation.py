# Set the base directory for datasets
# Change this path to point to your local data directory
BASE_DIR = "../data"

# 필요한 라이브러리 임포트
import os
import scanpy as sc
import numpy as np
import pandas as pd

os.chdir(f"{BASE_DIR}/splatter/Raw")
gt = sc.read_h5ad("GT.h5ad")
cell_idx, gene_idx = gt.obs.index.tolist(), gt.var.index.tolist()
sc.pp.normalize_total(gt, target_sum=1e4)
sc.pp.log1p(gt)
sc.pp.highly_variable_genes(gt, n_top_genes=2000)

sim1 = sc.read_h5ad("sim1.h5ad")
sc.pp.normalize_total(sim1, target_sum=1e4)
sc.pp.log1p(sim1)
sc.pp.highly_variable_genes(sim1, n_top_genes=2000)

sim2 = sc.read_h5ad("sim2.h5ad")
sc.pp.normalize_total(sim2, target_sum=1e4)
sc.pp.log1p(sim2)
sc.pp.highly_variable_genes(sim2, n_top_genes=2000)

sim3 = sc.read_h5ad("sim3.h5ad")
sc.pp.normalize_total(sim3, target_sum=1e4)
sc.pp.log1p(sim3)
sc.pp.highly_variable_genes(sim3, n_top_genes=2000)

 # imputed datasets
os.chdir(f"{BASE_DIR}/splatter/imputed")

# === 1. DCA (latent.tsv) 말고 일단 mean으로 가져옴 ===
#dca1_df = pd.read_csv("./DCA/sim1/mean.tsv", sep="\t", header=0, index_col=0).T
#dca1 = sc.AnnData(dca1_df)
#dca1.obs, dca1.var = sim1.obs.copy(), sim1.var.loc[dca1.var_names].copy()
#dca1.write_h5ad("./DCA/dca1.h5ad")# 최종 저장 및 불러오기
dca1 = sc.read_h5ad("./DCA/dca1.h5ad")
sc.pp.normalize_total(dca1, target_sum=1e4)
sc.pp.log1p(dca1)
sc.pp.highly_variable_genes(dca1, n_top_genes=2000)

#dca2_df = pd.read_csv("./DCA/sim2/mean.tsv", sep="\t", header=0, index_col=0).T
#dca2 = sc.AnnData(dca2_df)
#dca2.obs, dca2.var = sim2.obs.copy(), sim2.var.loc[dca2.var_names].copy()
#dca2.write_h5ad("./DCA/dca2.h5ad")# 최종 저장 및 불러오기
dca2 = sc.read_h5ad("./DCA/dca2.h5ad")
sc.pp.normalize_total(dca2, target_sum=1e4)
sc.pp.log1p(dca2)
sc.pp.highly_variable_genes(dca2, n_top_genes=2000)

#dca3_df = pd.read_csv("./DCA/sim3/mean.tsv", sep="\t", header=0, index_col=0).T
#dca3 = sc.AnnData(dca3_df)
#dca3.obs, dca3.var = sim3.obs.copy(), sim3.var.loc[dca3.var_names].copy()
#dca3.write_h5ad("./DCA/dca3.h5ad")# 최종 저장 및 불러오기
dca3 = sc.read_h5ad("./DCA/dca3.h5ad")
sc.pp.normalize_total(dca3, target_sum=1e4)
sc.pp.log1p(dca3)
sc.pp.highly_variable_genes(dca3, n_top_genes=2000)

# === 2. DrImpute ===
#drimpute1_df = pd.read_csv("./DrImpute/sim1_DrImpute.csv", sep=",", index_col=0).T
#drimpute1 = sc.AnnData(drimpute1_df)
#drimpute1.obs, drimpute1.var = sim1.obs.copy(), sim1.var.copy()
#drimpute1.write_h5ad("./DrImpute/drimpute1.h5ad")
drimpute1 = sc.read_h5ad("./DrImpute/drimpute1.h5ad")
sc.pp.highly_variable_genes(drimpute1, n_top_genes=2000)

#drimpute2_df = pd.read_csv("./DrImpute/sim2_DrImpute.csv", sep=",", index_col=0).T
#drimpute2 = sc.AnnData(drimpute2_df)
#drimpute2.obs, drimpute2.var = sim2.obs.copy(), sim2.var.copy()
#drimpute2.write_h5ad("./DrImpute/drimpute2.h5ad")
drimpute2 = sc.read_h5ad("./DrImpute/drimpute2.h5ad")
sc.pp.highly_variable_genes(drimpute2, n_top_genes=2000)

#drimpute3_df = pd.read_csv("./DrImpute/sim3_DrImpute.csv", sep=",", index_col=0).T
#drimpute3 = sc.AnnData(drimpute3_df)
#drimpute3.obs, drimpute3.var = sim3.obs.copy(), sim3.var.copy()
#drimpute3.write_h5ad("./DrImpute/drimpute3.h5ad")
drimpute3 = sc.read_h5ad("./DrImpute/drimpute3.h5ad")
sc.pp.highly_variable_genes(drimpute3, n_top_genes=2000)

# === 3. MAGIC ===
magic1 = sc.read_h5ad("./MAGIC/sim1_MAGIC.h5ad")
sc.pp.highly_variable_genes(magic1, n_top_genes=2000)

magic2 = sc.read_h5ad("./MAGIC/sim2_MAGIC.h5ad")
sc.pp.highly_variable_genes(magic2, n_top_genes=2000)

magic3 = sc.read_h5ad("./MAGIC/sim3_MAGIC.h5ad")
sc.pp.highly_variable_genes(magic3, n_top_genes=2000)

# === 4. SAVER ===
#saver1_df = pd.read_csv("./SAVER/sim1_saver_imputed.csv", sep=",", index_col=0).T
#saver1 = sc.AnnData(saver1_df)
#saver1.obs, saver1.var = sim1.obs.copy(), sim1.var.copy()
#saver1.write_h5ad("./SAVER/saver1.h5ad")
saver1 = sc.read_h5ad("./SAVER/saver1.h5ad")
sc.pp.normalize_total(saver1, target_sum=1e4)
sc.pp.log1p(saver1)
sc.pp.highly_variable_genes(saver1, n_top_genes=2000)

#saver2_df = pd.read_csv("./SAVER/sim2_saver_imputed.csv", sep=",", index_col=0).T
#saver2 = sc.AnnData(saver2_df)
#saver2.obs, saver2.var = sim2.obs.copy(), sim2.var.copy()
#saver2.write_h5ad("./SAVER/saver2.h5ad")
saver2 = sc.read_h5ad("./SAVER/saver2.h5ad")
sc.pp.normalize_total(saver2, target_sum=1e4)
sc.pp.log1p(saver2)
sc.pp.highly_variable_genes(saver2, n_top_genes=2000)

#saver3_df = pd.read_csv("./SAVER/sim3_saver_imputed.csv", sep=",", index_col=0).T
#saver3 = sc.AnnData(saver3_df)
#saver3.obs, saver3.var = sim3.obs.copy(), sim3.var.copy()
#saver3.write_h5ad("./SAVER/saver3.h5ad")
saver3 = sc.read_h5ad("./SAVER/saver3.h5ad")
sc.pp.normalize_total(saver3, target_sum=1e4)
sc.pp.log1p(saver3)
sc.pp.highly_variable_genes(saver3, n_top_genes=2000)

# === 5. scIDPMs ===
#scidpms1_df = pd.read_csv("./scIDPMs/sim1/imputed.csv", sep=",", index_col=None, header=None)
#scidpms1 = sc.AnnData(scidpms1_df)
#scidpms1.obs, scidpms1.var = sim1.obs.copy(), sim1.var.copy()
#scidpms1.write_h5ad("./scIDPMs/scidpms1.h5ad")
scidpms1 = sc.read_h5ad("./scIDPMs/scidpms1.h5ad")
sc.pp.normalize_total(scidpms1, target_sum=1e4)
sc.pp.log1p(scidpms1)
sc.pp.highly_variable_genes(scidpms1, n_top_genes=2000)

#scidpms2_df = pd.read_csv("./scIDPMs/sim2/imputed.csv", sep=",", index_col=None, header=None)
#scidpms2 = sc.AnnData(scidpms2_df)
#scidpms2.obs, scidpms2.var = sim2.obs.copy(), sim2.var.copy()
#scidpms2.write_h5ad("./scIDPMs/scidpms2.h5ad")
scidpms2 = sc.read_h5ad("./scIDPMs/scidpms2.h5ad")
sc.pp.normalize_total(scidpms2, target_sum=1e4)
sc.pp.log1p(scidpms2)
sc.pp.highly_variable_genes(scidpms2, n_top_genes=2000)

#scidpms3_df = pd.read_csv("./scIDPMs/sim3/imputed.csv", sep=",", index_col=None, header=None)
#scidpms3 = sc.AnnData(scidpms3_df)
#scidpms3.obs, scidpms3.var = sim3.obs.copy(), sim3.var.copy()
#scidpms3.write_h5ad("./scIDPMs/scidpms3.h5ad")
scidpms3 = sc.read_h5ad("./scIDPMs/scidpms3.h5ad")
sc.pp.normalize_total(scidpms3, target_sum=1e4)
sc.pp.log1p(scidpms3)
sc.pp.highly_variable_genes(scidpms3, n_top_genes=2000)

# === 6. scIGANs ===
#scigans1_df = pd.read_csv("./scIGANs/scIGANs_sim1.txt", sep="\t", index_col=0).T
#scigans1 = sc.AnnData(scigans1_df)
#scigans1.obs, scigans1.var = sim1.obs.copy(), sim1.var.copy()
#scigans1.write_h5ad("./scIGANs/scigans1.h5ad")
scigans1 = sc.read_h5ad("./scIGANs/scigans1.h5ad")
sc.pp.normalize_total(scigans1, target_sum=1e4)
sc.pp.log1p(scigans1)
sc.pp.highly_variable_genes(scigans1, n_top_genes=2000)

#scigans2_df = pd.read_csv("./scIGANs/scIGANs_sim2.txt", sep="\t", index_col=0).T
#scigans2 = sc.AnnData(scigans2_df)
#scigans2.obs, scigans2.var = sim2.obs.copy(), sim2.var.copy()
#scigans2.write_h5ad("./scIGANs/scigans2.h5ad")
scigans2 = sc.read_h5ad("./scIGANs/scigans2.h5ad")
sc.pp.normalize_total(scigans2, target_sum=1e4)
sc.pp.log1p(scigans2)
sc.pp.highly_variable_genes(scigans2, n_top_genes=2000)

#scigans3_df = pd.read_csv("./scIGANs/scIGANs_sim3.txt", sep="\t", index_col=0).T
#scigans3 = sc.AnnData(scigans3_df)
#scigans3.obs, scigans3.var = sim3.obs.copy(), sim3.var.copy()
#scigans3.write_h5ad("./scIGANs/scigans3.h5ad")
scigans3 = sc.read_h5ad("./scIGANs/scigans3.h5ad")
sc.pp.normalize_total(scigans3, target_sum=1e4)
sc.pp.log1p(scigans3)
sc.pp.highly_variable_genes(scigans3, n_top_genes=2000)

# === 8. scMultiGAN ===
#scmultigan1_df = pd.read_csv("./scMultiGAN/sim1/Restored/sim1_scMultiGAN.tsv", sep="\t", header=None, index_col=None).T
#scmultigan1 = sc.AnnData(scmultigan1_df)
#scmultigan1.obs, scmultigan1.var = sim1.obs.copy(), sim1.var.copy()
#scmultigan1.write_h5ad("./scMultiGAN/scmultigan1.h5ad")
scmultigan1 = sc.read_h5ad("./scMultiGAN/scmultigan1.h5ad")
sc.pp.normalize_total(scmultigan1, target_sum=1e4)
sc.pp.log1p(scmultigan1)
sc.pp.highly_variable_genes(scmultigan1, n_top_genes=2000)

#scmultigan2_df = pd.read_csv("./scMultiGAN/sim2/Restored/sim2_scMultiGAN.tsv", sep="\t", header=None, index_col=None).T
#scmultigan2 = sc.AnnData(scmultigan2_df)
#scmultigan2.obs, scmultigan2.var = sim2.obs.copy(), sim2.var.copy()
#scmultigan2.write_h5ad("./scMultiGAN/scmultigan2.h5ad")
scmultigan2 = sc.read_h5ad("./scMultiGAN/scmultigan2.h5ad")
sc.pp.normalize_total(scmultigan2, target_sum=1e4)
sc.pp.log1p(scmultigan2)
sc.pp.highly_variable_genes(scmultigan2, n_top_genes=2000)

#scmultigan3_df = pd.read_csv("./scMultiGAN/sim3/Restored/sim3_scMultiGAN.tsv", sep="\t", header=None, index_col=None).T
#scmultigan3 = sc.AnnData(scmultigan3_df)
#scmultigan3.obs, scmultigan3.var = sim3.obs.copy(), sim3.var.copy()
#scmultigan3.write_h5ad("./scMultiGAN/scmultigan3.h5ad")
scmultigan3 = sc.read_h5ad("./scMultiGAN/scmultigan3.h5ad")
sc.pp.normalize_total(scmultigan3, target_sum=1e4)
sc.pp.log1p(scmultigan3)
sc.pp.highly_variable_genes(scmultigan3, n_top_genes=2000)


# === 9. scSTD ===
#scstd1_df = pd.read_csv("./scSTD/Result/sim1_imputed.txt", sep=",", index_col=None, header = None)
#scstd1 = sc.AnnData(scstd1_df)
#scstd1.obs, scstd1.var = sim1.obs.copy(), sim1.var.copy()
#scstd1.write_h5ad("./scSTD/scstd1.h5ad")
scstd1 = sc.read_h5ad("./scSTD/scstd1.h5ad")
sc.pp.normalize_total(scstd1, target_sum=1e4)
sc.pp.log1p(scstd1)
sc.pp.highly_variable_genes(scstd1, n_top_genes=2000)

#scstd2_df = pd.read_csv("./scSTD/Result/sim2_imputed.txt", sep=",", index_col=None, header = None)
#scstd2 = sc.AnnData(scstd2_df)
#scstd2.obs, scstd2.var = sim2.obs.copy(), sim2.var.copy()
#scstd2.write_h5ad("./scSTD/scstd2.h5ad")
scstd2 = sc.read_h5ad("./scSTD/scstd2.h5ad")
sc.pp.normalize_total(scstd2, target_sum=1e4)
sc.pp.log1p(scstd2)
sc.pp.highly_variable_genes(scstd2, n_top_genes=2000)

#scstd3_df = pd.read_csv("./scSTD/Result/sim3_imputed.txt", sep=",", index_col=None, header = None)
#scstd3 = sc.AnnData(scstd3_df)
#scstd3.obs, scstd3.var = sim3.obs.copy(), sim3.var.copy()
#scstd3.write_h5ad("./scSTD/scstd3.h5ad")
scstd3 = sc.read_h5ad("./scSTD/scstd3.h5ad")
sc.pp.normalize_total(scstd3, target_sum=1e4)
sc.pp.log1p(scstd3)
sc.pp.highly_variable_genes(scstd3, n_top_genes=2000)

# === 10. scVI ===
scvi1 = sc.read_h5ad("./scVI/sim1_scVI.h5ad")
scvi1.obs_names, scvi1.var_names = cell_idx, gene_idx
print("scvi1 shape:", scvi1.shape)
sc.pp.highly_variable_genes(scvi1, n_top_genes=2000)

scvi2 = sc.read_h5ad("./scVI/sim2_scVI.h5ad")
scvi2.obs_names, scvi2.var_names = cell_idx, gene_idx
print("scvi2 shape:", scvi2.shape)
sc.pp.highly_variable_genes(scvi2, n_top_genes=2000)

scvi3 = sc.read_h5ad("./scVI/sim3_scVI.h5ad")
scvi3.obs_names, scvi3.var_names = cell_idx, gene_idx
print("scvi shape:", scvi3.shape)
sc.pp.highly_variable_genes(scvi3, n_top_genes=2000)

# === 11. ALRA ===
#alra1_df =  pd.read_csv("./ALRA/sim1_ALRA_imputed.csv", index_col= 0)
#alra1 = sc.AnnData(alra1_df)
#alra1.obs, alra1.var = sim1.obs.copy(), sim1.var.copy()
#alra1.obs_names, alra1.var_names = cell_idx, gene_idx
#alra1.write_h5ad("./ALRA/alra1.h5ad")
alra1 = sc.read_h5ad("./ALRA/alra1.h5ad")
sc.pp.highly_variable_genes(alra1, n_top_genes=2000)

#alra2_df =  pd.read_csv("./ALRA/sim2_ALRA_imputed.csv", index_col= 0)
#alra2 = sc.AnnData(alra2_df)
#alra2.obs, alra2.var = sim2.obs.copy(), sim2.var.copy()
#alra2.obs_names, alra2.var_names = cell_idx, gene_idx
#alra2.write_h5ad("./ALRA/alra2.h5ad")
alra2 = sc.read_h5ad("./ALRA/alra2.h5ad")
sc.pp.highly_variable_genes(alra2, n_top_genes=2000)

#alra3_df =  pd.read_csv("./ALRA/sim3_ALRA_imputed.csv", index_col= 0)
#alra3 = sc.AnnData(alra3_df)
#alra3.obs, alra3.var = sim3.obs.copy(), sim3.var.copy()
#alra3.obs_names, alra3.var_names = cell_idx, gene_idx
#alra3.write_h5ad("./ALRA/alra3.h5ad")
alra3 = sc.read_h5ad("./ALRA/alra3.h5ad")
sc.pp.highly_variable_genes(alra3, n_top_genes=2000)

# Clustering & Evaluation utilities ((260212 수정본))
# =========================================================

import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

# =========================================================
# 1. Purity
# =========================================================
def purity_score(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    
    labels_true = np.unique(y_true)
    labels_pred = np.unique(y_pred)
    
    C = np.zeros((len(labels_true), len(labels_pred)), dtype=int)
    for i, lt in enumerate(labels_true):
        for j, lp in enumerate(labels_pred):
            C[i, j] = np.sum((y_true == lt) & (y_pred == lp))
    
    return np.sum(np.max(C, axis=0)) / np.sum(C)


# =========================================================
# 2. Leiden resolution 탐색 (cluster 수 고정)
# =========================================================
def find_resolution_for_k(ad, target_k=8, low=0.1, high=2.0, max_iter=20):
    for i in range(max_iter):
        res = (low + high) / 2
        # 경고 메시지의 권장 사항(flavor, n_iterations, directed)을 반영합니다.
        sc.tl.leiden(
            ad, 
            resolution=res, 
            key_added="_leiden_tmp",
            flavor="igraph", 
            n_iterations=2, 
            directed=False
        )
        
        k = ad.obs["_leiden_tmp"].nunique()
        print(f"Iteration {i+1}: res={res:.4f}, k={k}") # 진행 상황 확인용
    
        if k > target_k:
            high = res
        elif k < target_k:
            low = res
        else:
            return res
    return res

# =========================================================
# 3. Clustering + metric 계산
# =========================================================
def run_clustering_metrics(
    adata,
    gt_key="gt",
    target_k=8,
    n_neighbors=15,
    n_pcs=50,
    metric="euclidean",
    use_umap=True,
):
    ad = adata.copy()
    
    # 전처리
    sc.pp.scale(ad, max_value=10)
    n_pcs_eff = min(n_pcs, ad.n_vars - 1) if ad.n_vars > 1 else 1
    sc.tl.pca(ad, n_comps=n_pcs_eff)
    sc.pp.neighbors(ad, n_neighbors=n_neighbors, n_pcs=n_pcs_eff, metric=metric)
    
    # ---------------------------------------------------------
    # [수정 포인트] resolution만 가져오는 게 아니라 결과 자체를 고정
    # ---------------------------------------------------------
    res = find_resolution_for_k(ad, target_k=target_k)
    
    # find_resolution_for_k 내부에서 마지막으로 k=8을 만든 결과가 
    # 이미 ad.obs["_leiden_tmp"]에 들어있습니다. 이걸 그대로 복사합니다.
    ad.obs["leiden"] = ad.obs["_leiden_tmp"].copy()
    
    # 만약 끝내 k=8을 못 찾고 max_iter가 끝났을 경우를 대비해 경고만 남깁니다.
    current_k = ad.obs["leiden"].nunique()
    if current_k != target_k:
        print(f"      [Warning] Could not find exact k={target_k}. Found k={current_k} instead.")
    # ---------------------------------------------------------
    
    if use_umap:
        sc.tl.umap(ad)
    
    y_true = ad.obs[gt_key].astype(str).values
    y_pred = ad.obs["leiden"].astype(str).values
    
    scores = {
        "ARI": adjusted_rand_score(y_true, y_pred),
        "NMI": normalized_mutual_info_score(y_true, y_pred),
        "Purity": purity_score(y_true, y_pred),
    }
    
    return ad, scores

# =========================================================
# 4. score dict → tidy dataframe
# =========================================================
def tidy_scores(scores_dict):
    return (
        pd.DataFrame(scores_dict)
        .T[["ARI", "NMI", "Purity"]]
    )


# 2. Clustering 실행 및 점수 계산 (260212 수정본)
# =========================================================

print("--- 1. 클러스터링용 데이터 구성 시작 ---")

# -----------------------------
# Ground Truth label 설정
# -----------------------------
if "Group" not in gt.obs:
    raise KeyError("gt.obs['Group']를 찾을 수 없습니다.")

gt_labels = gt.obs["Group"].astype(str)
k_true = gt_labels.nunique()
print(f"Ground Truth 레이블 로드 완료. (k={k_true})")

if k_true != 8:
    print(f"경고: Ground Truth k={k_true}, target_k=8")

# -----------------------------
# simulation별 AnnData 묶기
# -----------------------------
sims = {
    "sim1": {
        "RAW": sim1,
        "ALRA": alra1,
        "DrImpute": drimpute1,
        "MAGIC": magic1,
        "SAVER": saver1,
        "scVI": scvi1,
        "DCA": dca1,
        "scIGANs": scigans1,
        "scMultiGAN": scmultigan1,
        "scSTD": scstd1,
        "scIDPMs": scidpms1,
    },
    "sim2": {
        "RAW": sim2,
        "ALRA": alra2,
        "DrImpute": drimpute2,
        "MAGIC": magic2,
        "SAVER": saver2,
        "scVI": scvi2,
        "DCA": dca2,
        "scIGANs": scigans2,
        "scMultiGAN": scmultigan2,
        "scSTD": scstd2,
        "scIDPMs": scidpms2,
    },
    "sim3": {
        "RAW": sim3,
        "ALRA": alra3,
        "DrImpute": drimpute3,
        "MAGIC": magic3,
        "SAVER": saver3,
        "scVI": scvi3,
        "DCA": dca3,
        "scIGANs": scigans3,
        "scMultiGAN": scmultigan3,
        "scSTD": scstd3,
        "scIDPMs": scidpms3,
    },
}

print("Simulation 객체 묶기 완료.")

# =========================================================
# 2. Clustering & Evaluation
# =========================================================
print("\n--- 2. 클러스터링 및 평가 시작 ---")

all_scores = {}
all_ads = {}

for sim_name, tools in sims.items():
    print(f"\n[{sim_name}] 처리 중...")
    all_scores[sim_name] = {}
    all_ads[sim_name] = {}
    
    for tool_name, adata in tools.items():
        print(f"  -> {tool_name} 실행 중...")
    
        # -----------------------------
        # GT label 주입 (정합은 이미 완료된 상태)
        # -----------------------------
        ad = adata.copy()
        ad.obs["gt"] = gt_labels.loc[ad.obs_names].values
    
        # -----------------------------
        # clustering + metric
        # -----------------------------
        ad_res, scores = run_clustering_metrics(
            adata=ad,
            gt_key="gt",
            target_k=8,
            n_neighbors=15,
            n_pcs=50,
            metric="euclidean",
            use_umap=True,
        )
    
        all_scores[sim_name][tool_name] = scores
        all_ads[sim_name][tool_name] = ad_res
    
        print(f"     {tool_name} scores: {scores}")

print("\n--- 모든 클러스터링 및 평가 완료 ---")

# =========================================================
# 3. Score table 정리
# =========================================================
score_tables = {
    sim_name: tidy_scores(score_dict)
    for sim_name, score_dict in all_scores.items()
}

for sim_name, df in score_tables.items():
    print(f"\n--- {sim_name} Scores ---")
    print(df)


# 3. UMAP 시각화 + 저장 (260212 수정본)
# =========================================================
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns

SAVE_DIR = f"{BASE_DIR}/splatter/imputed/figures"
os.makedirs(SAVE_DIR, exist_ok=True)

print("\n--- 3. UMAP 시각화 시작 ---")

sim_names = ["sim1", "sim2", "sim3"]
tools_to_plot = list(all_ads["sim1"].keys())
print("UMAP에 포함될 모델:", tools_to_plot)


def plot_umap_grid(sim_name, tools, all_ads, n_cols=5):
    n_items = len(tools)
    n_rows = int(np.ceil(n_items / n_cols))
    
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows)
    )
    axes = axes.flatten()
    fig.suptitle(f"UMAP - {sim_name} (colored by GT)", fontsize=18)
    
    for i, tool in enumerate(tools):
        ad = all_ads[sim_name][tool]
        sc.pl.umap(
            ad,
            color="gt",
            title=tool,
            ax=axes[i],
            show=False,
            legend_loc="on data",
            legend_fontsize=6,
        )
    
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_path = os.path.join(SAVE_DIR, f"UMAP_{sim_name}.png")
    plt.savefig(save_path, dpi=300)
    print(f"[저장됨] {save_path}")
    plt.show()


for sim in sim_names:
    print(f"  -> {sim} UMAP 생성 중...")
    plot_umap_grid(sim, tools_to_plot, all_ads)

print("UMAP 시각화 완료.")


# =========================================================
# 4. Clustering metric Barplot + 저장 (RAW 기준선)
# =========================================================
print("\n--- 4. 지표 Barplot 시각화 시작 ---")

for sim_name, df_scores in score_tables.items():
    
    df_plot = df_scores.reset_index().rename(columns={"index": "tool"})
    metrics = ["ARI", "NMI", "Purity"]
    
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5))
    fig.suptitle(f"Clustering Metrics - {sim_name.upper()}", fontsize=16)
    
    for i, metric in enumerate(metrics):
        ax = axes[i]
        
        sns.barplot(
            data=df_plot,
            x="tool",
            y=metric,
            ax=ax,
            palette='deep'
        )
    
        if "RAW" in df_plot["tool"].values:
            raw_val = df_plot.loc[df_plot["tool"] == "RAW", metric].values[0]
            ax.axhline(
                raw_val,
                color="black",
                linestyle="--",
                linewidth=1,
                alpha=0.8,
            )
    
        ax.set_title(metric)
        ax.set_xlabel("")
        ax.set_ylabel("Score")
        ax.set_ylim(0, 1.05)
        ax.tick_params(axis="x", rotation=45)
    
    plt.tight_layout()
    save_path = os.path.join(SAVE_DIR, f"Barplot_{sim_name}.png")
    plt.savefig(save_path, dpi=300)
    print(f"[저장됨] {save_path}")
    plt.show()

print("지표 Barplot 시각화 완료.")


# FC-SCC 분석용 유틸리티
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
import scanpy as sc

SAVE_DIR = f"{BASE_DIR}/splatter/imputed/figures"
os.makedirs(SAVE_DIR, exist_ok=True)
# 필요한 라이브러리 임포트
import os
import scanpy as sc
import numpy as np
import pandas as pd

# GT raw cnt로 다시 불러오기
os.chdir(f"{BASE_DIR}/splatter/Raw")
gt = sc.read_h5ad("GT.h5ad")
cell_idx, gene_idx = gt.obs.index.tolist(), gt.var.index.tolist()
# -------------------------------------------------------
# GT raw count → log2(CPM+1) 변환
# -------------------------------------------------------
def raw_to_log2cpm(adata):
    """raw count AnnData → log2(CPM+1) numpy array"""
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    lib_size = X.sum(axis=1, keepdims=True)
    lib_size = np.where(lib_size == 0, 1, lib_size)
    cpm = X / lib_size * 1e6
    return np.log2(cpm + 1)

# -------------------------------------------------------
# GT log2CPM 기반 signal gene 선별
# -------------------------------------------------------
def get_signal_genes(gt_log2, groups, lfc_cutoff=1.0):
    """
    GT에서 클러스터 간 |log2FC| >= lfc_cutoff인 gene index 반환
    (모든 pair 중 하나라도 통과하면 signal gene)
    """
    groups = np.asarray(groups)
    uniq = np.unique(groups) 
    
    pb = {}
    for g in uniq:
        idx = np.where(groups == g)[0]
        pb[g] = np.median(gt_log2[idx, :], axis=0)
    
    signal_mask = np.zeros(gt_log2.shape[1], dtype=bool)
    for g1, g2 in itertools.combinations(uniq, 2):
        lfc = np.abs(pb[g1] - pb[g2])
        signal_mask |= (lfc >= lfc_cutoff)
    
    print(f"Signal genes selected: {signal_mask.sum()} / {len(signal_mask)}")
    return signal_mask

# -------------------------------------------------------
# signal gene subset으로 Spearman FC-SCC 계산
# -------------------------------------------------------
def compute_fc_scc(gt_mat, imp_mat, groups):
    """이미 signal gene subset으로 잘린 matrix를 받음"""
    groups = np.asarray(groups)
    uniq   = np.unique(groups)
    
    pb_gt  = {g: np.median(gt_mat[groups == g, :],  axis=0) for g in uniq}
    pb_imp = {g: np.median(imp_mat[groups == g, :], axis=0) for g in uniq}
    
    fc_scc = {}
    for g1, g2 in itertools.combinations(uniq, 2):
        lfc_true = pb_gt[g1]  - pb_gt[g2]
        lfc_imp  = pb_imp[g1] - pb_imp[g2]
        r, _ = spearmanr(lfc_true, lfc_imp)
        fc_scc[(g1, g2)] = 0.0 if np.isnan(r) else r
    
    return fc_scc

# -------------------------------------------------------
# Main 분석 함수
# -------------------------------------------------------
def run_fc_scc_heatmap(sim_name, sim_data, gt_raw_adata, lfc_cutoff=1.0):
    print(f"\n--- FC-SCC analysis: {sim_name} (log2FC cutoff={lfc_cutoff}) ---")
    
    gt_log2  = raw_to_log2cpm(gt_raw_adata)
    groups   = np.asarray(gt_raw_adata.obs["Group"].astype(str))
    gt_genes = gt_raw_adata.var.index.tolist()
    
    # 모든 tool의 공통 gene 교집합 (GT 순서 기준)
    common_genes = set(gt_genes)
    for adata in sim_data.values():
        common_genes &= set(adata.var.index.tolist())
    common_genes = [g for g in gt_genes if g in common_genes]  # GT 순서 유지
    gt_col_idx   = [gt_genes.index(g) for g in common_genes]
    print(f"Common genes: {len(common_genes)}")
    
    # GT signal gene mask (공통 gene 공간에서 한 번만 계산)
    gt_sub      = gt_log2[:, gt_col_idx]
    signal_mask = get_signal_genes(gt_sub, groups, lfc_cutoff=lfc_cutoff)
    gt_signal   = gt_sub[:, signal_mask]
    
    # 각 tool별 FC-SCC 계산
    results = {}
    for tool_name, adata in sim_data.items():
        imp_genes     = adata.var.index.tolist()
        imp_gene_dict = {g: i for i, g in enumerate(imp_genes)}
        imp_col_idx   = [imp_gene_dict[g] for g in common_genes]
    
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        imp_signal = np.asarray(X, dtype=np.float64)[:, imp_col_idx][:, signal_mask]
    
        results[tool_name] = compute_fc_scc(gt_signal, imp_signal, groups)
    
    # DataFrame 정리
    pair_labels = sorted(next(iter(results.values())).keys())
    pair_str    = [f"{a}_vs_{b}" for a, b in pair_labels]
    
    fc_df = pd.DataFrame(index=list(results.keys()), columns=pair_str, dtype=float)
    for tool, fc_dict in results.items():
        for (a, b), v in fc_dict.items():
            fc_df.loc[tool, f"{a}_vs_{b}"] = v
    
    # Heatmap 저장
    fig, ax = plt.subplots(figsize=(max(12, 0.75 * len(pair_str)),
                                    max(6,  0.5  * len(results) + 3)))
    sns.heatmap(fc_df, annot=True, fmt=".2f",
                cmap="YlOrRd", vmin=0.3, vmax=1.0, ax=ax)
    ax.set_title(f"FC-SCC ({sim_name}) | GT log2FC≥{lfc_cutoff} signal genes")
    ax.set_ylabel("Tool")
    ax.set_xlabel("Celltype pair")
    plt.tight_layout()
    
    out_fn = os.path.join(SAVE_DIR, f"{sim_name}_FC_SCC_heatmap.png")
    fig.savefig(out_fn, dpi=300)
    plt.close(fig)
    print(f"[saved] {out_fn}")
    
    return fc_df

# -------------------------------------------------------
# 실행
# -------------------------------------------------------
all_results = {}
for sim_name, sim_data in sims.items():
    fc_df = run_fc_scc_heatmap(
        sim_name    = sim_name,
        sim_data    = sim_data,
        gt_raw_adata= gt,
        lfc_cutoff  = 1.0
    )
    all_results[sim_name] = fc_df
    print(f"\n[{sim_name}] FC-SCC DataFrame:")
    print(fc_df)

# 결과 집계 및 툴별 랭킹
# -------------------------------------------------------

# 1. sim별 fc_df의 평균(pair 평균) → 툴별 mean SCC
summary_rows = []
for sim_name, fc_df in all_results.items():
    mean_per_tool = fc_df.mean(axis=1)  # 각 tool의 pair 평균
    for tool, val in mean_per_tool.items():
        summary_rows.append({"sim": sim_name, "tool": tool, "mean_SCC": val})

summary_df = pd.DataFrame(summary_rows)

# 2. sim 전체 평균
overall_mean = (
    summary_df.groupby("tool")["mean_SCC"]
    .mean()
    .reset_index()
    .rename(columns={"mean_SCC": "overall_mean_SCC"})
    .sort_values("overall_mean_SCC", ascending=False)
    .reset_index(drop=True)
)
overall_mean["rank"] = overall_mean.index + 1

print("\n=== Tool Ranking (FC-SCC, averaged over all sims & pairs) ===")
print(overall_mean.to_string(index=False))

# 3. 시각화: 툴별 mean SCC (sim별 breakdown + overall)
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 왼쪽: sim별 barplot
pivot_df = summary_df.pivot(index="tool", columns="sim", values="mean_SCC")
pivot_df = pivot_df.loc[overall_mean["tool"]]  # overall 랭킹 순서로 정렬
pivot_df.plot(kind="bar", ax=axes[0], colormap="Set2", edgecolor="black", linewidth=0.5)
axes[0].set_title("Mean FC-SCC per Tool (by Sim)")
axes[0].set_ylabel("Mean Spearman r")
axes[0].set_xlabel("")
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=45, ha="right")
axes[0].legend(title="Sim", bbox_to_anchor=(1.01, 1), loc="upper left")
axes[0].axhline(0, color="gray", linewidth=0.8, linestyle="--")

# 오른쪽: overall mean + 랭킹 표시
bars = axes[1].barh(
    overall_mean["tool"][::-1],
    overall_mean["overall_mean_SCC"][::-1],
    color="steelblue", edgecolor="black", linewidth=0.5
)
for bar, val in zip(bars, overall_mean["overall_mean_SCC"][::-1]):
    axes[1].text(bar.get_width() + 0.005, bar.get_y() + bar.get_height() / 2,
                 f"{val:.3f}", va="center", fontsize=9)
axes[1].set_title("Overall Mean FC-SCC Ranking")
axes[1].set_xlabel("Mean Spearman r")
axes[1].axvline(0, color="gray", linewidth=0.8, linestyle="--")

plt.tight_layout()
out_fn = os.path.join(SAVE_DIR, "tool_ranking_FC_SCC.png")
fig.savefig(out_fn, dpi=300)
plt.close(fig)
print(f"\n[saved] {out_fn}")

# Cell-wise PCC 분석 함수 정의
os.chdir(f"{BASE_DIR}/splatter/Raw")
gt = sc.read_h5ad("GT.h5ad")
cell_idx, gene_idx = gt.obs.index.tolist(), gt.var.index.tolist()
sc.pp.normalize_total(gt, target_sum=1e4)
sc.pp.log1p(gt)

import gc
def run_cell_pcc_boxplot(sim_name, sim_data, gt_log1p_adata):
    """
    GT: log1p norm (imputed와 스케일 통일)
    Imputed: log1p norm
    → 공통 gene 교집합 기준 cell-wise Pearson r
    """
    print(f"\n--- Cell-wise PCC analysis: {sim_name} ---")
    
    X_gt = gt_log1p_adata.X
    if hasattr(X_gt, "toarray"):
        X_gt = X_gt.toarray()
    gt_log1p = np.asarray(X_gt, dtype=np.float64)
    gt_genes = gt_log1p_adata.var.index.tolist()
    
    # 모든 tool 공통 gene 교집합
    common_genes = set(gt_genes)
    for adata in sim_data.values():
        common_genes &= set(adata.var.index.tolist())
    common_genes  = [g for g in gt_genes if g in common_genes]
    gt_col_idx    = [gt_genes.index(g) for g in common_genes]
    gt_sub        = gt_log1p[:, gt_col_idx]
    print(f"Common genes: {len(common_genes)}")
    
    pcc_records = []
    
    for tool_name, adata in sim_data.items():
        imp_genes     = adata.var.index.tolist()
        imp_gene_dict = {g: i for i, g in enumerate(imp_genes)}
        imp_col_idx   = [imp_gene_dict[g] for g in common_genes]
    
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        imp_sub = np.asarray(X, dtype=np.float64)[:, imp_col_idx]
    
        n_cells = gt_sub.shape[0]
        corrs = np.array([
            np.corrcoef(gt_sub[i], imp_sub[i])[0, 1]
            for i in range(n_cells)
        ])
        corrs = np.where(np.isnan(corrs), 0.0, corrs)
    
        pcc_records.append(pd.DataFrame({"PCC": corrs, "Tool": tool_name}))
        print(f"  {tool_name}: median PCC = {np.median(corrs):.4f}")
    
    final_pcc_df = pd.concat(pcc_records, ignore_index=True)
    tool_order   = list(sim_data.keys())
    
    raw_med = None
    if "RAW" in final_pcc_df["Tool"].values:
        raw_med = final_pcc_df.loc[final_pcc_df["Tool"] == "RAW", "PCC"].median()
    
    fig, ax = plt.subplots(figsize=(max(6, len(tool_order) * 0.9), 5.8))
    sns.boxplot(
        data=final_pcc_df,
        x="Tool", y="PCC",
        order=tool_order,
        hue="Tool", hue_order=tool_order,
        dodge=False, palette="tab20",
        showfliers=False, ax=ax
    )
    if raw_med is not None:
        ax.axhline(raw_med, color="black", linestyle="--",
                   linewidth=1, alpha=0.8, label=f"RAW median ({raw_med:.3f})")
        ax.legend(fontsize=8)
    
    ax.set_title(f"Cell-wise PCC with GT - {sim_name}")
    ax.set_ylabel("Pearson r (cell-wise)")
    ax.set_xlabel("")
    ax.set_xticks(range(len(tool_order)))
    ax.set_xticklabels(tool_order, rotation=45, ha="right")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    
    out_fn = os.path.join(SAVE_DIR, f"{sim_name}_Cell_PCC_Boxplot.png")
    fig.savefig(out_fn, dpi=300)
    plt.close(fig)
    print(f"[saved] {out_fn}")
    
    del gt_sub, pcc_records
    gc.collect()
    
    return final_pcc_df

# 실행
pcc_results = {}
for sim_name, sim_data in sims.items():
    pcc_df = run_cell_pcc_boxplot(
        sim_name        = sim_name,
        sim_data        = sim_data,
        gt_log1p_adata  = gt
    )
    pcc_results[sim_name] = pcc_df

# -------------------------------------------------------
# Cell-wise PCC 랭킹
# -------------------------------------------------------
pcc_summary_rows = []
for sim_name, pcc_df in pcc_results.items():
    mean_per_tool = pcc_df.groupby("Tool")["PCC"].median()
    for tool, val in mean_per_tool.items():
        pcc_summary_rows.append({"sim": sim_name, "tool": tool, "median_PCC": val})

pcc_summary_df = pd.DataFrame(pcc_summary_rows)

pcc_ranking = (
    pcc_summary_df.groupby("tool")["median_PCC"]
    .mean()
    .reset_index()
    .rename(columns={"median_PCC": "overall_median_PCC"})
    .sort_values("overall_median_PCC", ascending=False)
    .reset_index(drop=True)
)
pcc_ranking["rank"] = pcc_ranking.index + 1

print("\n=== Tool Ranking (Cell-wise PCC) ===")
print(pcc_ranking.to_string(index=False))

# Pseudobulk SCC Scatter 분석 함수
def run_pseudobulk_scc_scatter(sim_name, sim_data, gt_log1p_adata):
    """
    Pseudobulk SCC Scatter:
    - 전체 gene을 클러스터 전체 median으로 aggregation (클러스터 구분 없이)
    - x축 = GT pseudobulk, y축 = Imputed pseudobulk
    - 한 점 = gene 하나, density 컬러링
    - 파일 하나에 tool별 subplot
    """
    print(f"\n--- Pseudobulk SCC Scatter: {sim_name} ---")
   
    X_gt = gt_log1p_adata.X
    if hasattr(X_gt, "toarray"):
        X_gt = X_gt.toarray()
    gt_log1p = np.asarray(X_gt, dtype=np.float64)
    gt_genes  = gt_log1p_adata.var.index.tolist()
   
    # 공통 gene 교집합
    common_genes = set(gt_genes)
    for adata in sim_data.values():
        common_genes &= set(adata.var.index.tolist())
    common_genes = [g for g in gt_genes if g in common_genes]
    gt_col_idx   = [gt_genes.index(g) for g in common_genes]
    gt_sub       = gt_log1p[:, gt_col_idx]
   
    # GT pseudobulk: 전체 cell median
    pb_gt = np.median(gt_sub, axis=0)  # shape: (n_genes,)
   
    # tool별 SCC 계산
    n_tools = len(sim_data)
    ncols   = 4
    nrows   = int(np.ceil(n_tools / ncols))
   
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 3.5, nrows * 3.5))
    axes = np.array(axes).flatten()
   
    summary_rows = []
   
    for i, (tool_name, adata) in enumerate(sim_data.items()):
        ax = axes[i]
   
        imp_genes     = adata.var.index.tolist()
        imp_gene_dict = {g: i for i, g in enumerate(imp_genes)}
        imp_col_idx   = [imp_gene_dict[g] for g in common_genes]
   
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        imp_sub = np.asarray(X, dtype=np.float64)[:, imp_col_idx]
        pb_imp  = np.median(imp_sub, axis=0)  # shape: (n_genes,)
   
        # density
        xy = np.vstack([pb_gt, pb_imp])
        try:
            kde     = gaussian_kde(xy)
            density = kde(xy)
        except Exception:
            density = np.ones(len(pb_gt))
   
        sc = ax.scatter(pb_gt, pb_imp, c=density, cmap="YlOrRd",
                        s=2, alpha=0.8, rasterized=True)
        plt.colorbar(sc, ax=ax, label="density")
   
        r, _ = spearmanr(pb_gt, pb_imp)
        r = 0.0 if np.isnan(r) else r
   
        ax.set_title(f"{tool_name}\nSCC = {r:.3f}", fontsize=9)
        ax.set_xlabel("GT pseudobulk (log1p)", fontsize=7)
        ax.set_ylabel(f"{tool_name} pseudobulk (log1p)", fontsize=7)
   
        lims = [min(pb_gt.min(), pb_imp.min()), max(pb_gt.max(), pb_imp.max())]
        ax.plot(lims, lims, "k--", linewidth=0.7, alpha=0.5)
   
        summary_rows.append({"sim": sim_name, "tool": tool_name, "pseudobulk_SCC": r})
   
    # 남은 subplot 비우기
    for j in range(n_tools, len(axes)):
        axes[j].set_visible(False)
   
    fig.suptitle(f"Pseudobulk SCC Scatter - {sim_name}", fontsize=12)
    plt.tight_layout()
   
    out_fn = os.path.join(SAVE_DIR, f"{sim_name}_Pseudobulk_SCC_Scatter.png")
    fig.savefig(out_fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out_fn}")
   
    return pd.DataFrame(summary_rows)

# -------------------------------------------------------
# 실행
# -------------------------------------------------------
pb_scc_scatter_results = {}
all_summary_rows = []

for sim_name, sim_data in sims.items():
    summary_df = run_pseudobulk_scc_scatter(
        sim_name       = sim_name,
        sim_data       = sim_data,
        gt_log1p_adata = gt
    )
    pb_scc_scatter_results[sim_name] = summary_df
    all_summary_rows.append(summary_df)

# 전체 요약 + 랭킹
all_summary = pd.concat(all_summary_rows, ignore_index=True)

ranking = (
    all_summary.groupby("tool")["pseudobulk_SCC"]
    .mean()
    .reset_index()
    .rename(columns={"pseudobulk_SCC": "mean_SCC"})
    .sort_values("mean_SCC", ascending=False)
    .reset_index(drop=True)
)
ranking["rank"] = ranking.index + 1

print("\n=== Pseudobulk SCC per Sim ===")
print(all_summary.pivot(index="tool", columns="sim", values="pseudobulk_SCC")
      .round(3).to_string())

print("\n=== Tool Ranking (Pseudobulk SCC) ===")
print(ranking.to_string(index=False))

# 6) Cell-wise AUPRC 분석 함수 정의
from sklearn.metrics import average_precision_score

def run_cell_auprc_boxplot(sim_name, sim_data, gt_log1p_adata, n_sample_cells=500):
    """
    Cell-wise AUPRC:
    - GT log1p > 0 → 1 (발현), = 0 → 0 (dropout) 으로 이진화
    - Imputed log1p 값을 score로 사용
    - 각 cell마다 AUPRC 계산
    """
    print(f"\n--- Cell-wise AUPRC: {sim_name} ---")
    
    X_gt = gt_log1p_adata.X
    if hasattr(X_gt, "toarray"):
        X_gt = X_gt.toarray()
    gt_log1p = np.asarray(X_gt, dtype=np.float64)
    gt_genes  = gt_log1p_adata.var.index.tolist()
    gt_cells  = gt_log1p_adata.obs.index.tolist()
    
    # 공통 gene 교집합
    common_genes = set(gt_genes)
    for adata in sim_data.values():
        common_genes &= set(adata.var.index.tolist())
    common_genes = [g for g in gt_genes if g in common_genes]
    gt_col_idx   = [gt_genes.index(g) for g in common_genes]
    gt_sub       = gt_log1p[:, gt_col_idx]
    print(f"Common genes: {len(common_genes)}")
    
    # GT 이진화
    gt_binary = (gt_sub > 0).astype(np.int8)  # (n_cells, n_genes)
    
    # cell sampling (모든 tool에 동일하게 적용)
    np.random.seed(42)
    n_cells = gt_binary.shape[0]
    if n_cells > n_sample_cells:
        sampled_idx = np.random.choice(n_cells, n_sample_cells, replace=False)
        sampled_idx = np.sort(sampled_idx)
    else:
        sampled_idx = np.arange(n_cells)
    print(f"Sampled cells: {len(sampled_idx)}")
    
    gt_binary_sub = gt_binary[sampled_idx, :]  # (n_sample, n_genes)
    
    auprc_records = []
    
    for tool_name, adata in sim_data.items():
        imp_genes     = adata.var.index.tolist()
        imp_gene_dict = {g: i for i, g in enumerate(imp_genes)}
        imp_col_idx   = [imp_gene_dict[g] for g in common_genes]
    
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        imp_sub = np.asarray(X, dtype=np.float64)[sampled_idx, :][:, imp_col_idx]
    
        # cell-wise AUPRC
        scores = []
        for r in range(gt_binary_sub.shape[0]):
            y_true  = gt_binary_sub[r]
            y_score = imp_sub[r]
            y_score = np.nan_to_num(y_score, nan=0.0, posinf=0.0, neginf=0.0)
            if y_true.sum() == 0:
                scores.append(np.nan)
            else:
                scores.append(average_precision_score(y_true, y_score))
    
        scores = np.array(scores)
        auprc_records.append(pd.DataFrame({"AUPRC": scores, "Tool": tool_name}))
        print(f"  {tool_name}: median AUPRC = {np.nanmedian(scores):.4f}")
    
        del imp_sub
        gc.collect()
    
    final_auprc_df = pd.concat(auprc_records, ignore_index=True)
    tool_order = list(sim_data.keys())
    
    raw_med = None
    if "RAW" in final_auprc_df["Tool"].values:
        raw_med = final_auprc_df.loc[final_auprc_df["Tool"] == "RAW", "AUPRC"].median()
    
    fig, ax = plt.subplots(figsize=(max(6, len(tool_order) * 0.9), 5.8))
    sns.boxplot(
        data=final_auprc_df,
        x="Tool", y="AUPRC",
        order=tool_order,
        hue="Tool", hue_order=tool_order,
        dodge=False, palette="tab20",
        showfliers=False, ax=ax
    )
    if raw_med is not None:
        ax.axhline(raw_med, color="black", linestyle="--",
                   linewidth=1, alpha=0.8, label=f"RAW median ({raw_med:.3f})")
        ax.legend(fontsize=8)
    
    ax.set_title(f"Cell-wise AUPRC - {sim_name} (n={len(sampled_idx)} cells)")
    ax.set_ylabel("AUPRC")
    ax.set_xlabel("")
    ax.set_xticks(range(len(tool_order)))
    ax.set_xticklabels(tool_order, rotation=45, ha="right")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    
    out_fn = os.path.join(SAVE_DIR, f"{sim_name}_Cell_AUPRC_Boxplot.png")
    fig.savefig(out_fn, dpi=300)
    plt.close(fig)
    print(f"[saved] {out_fn}")
    
    del gt_binary_sub, auprc_records
    gc.collect()
    
    return final_auprc_df

# 실행
auprc_results = {}
for sim_name, sim_data in sims.items():
    auprc_df = run_cell_auprc_boxplot(
        sim_name       = sim_name,
        sim_data       = sim_data,
        gt_log1p_adata = gt,
        n_sample_cells = 500
    )
    auprc_results[sim_name] = auprc_df

# 랭킹
auprc_rows = []
for sim_name, auprc_df in auprc_results.items():
    for tool, val in auprc_df.groupby("Tool")["AUPRC"].median().items():
        auprc_rows.append({"sim": sim_name, "tool": tool, "median_AUPRC": val})

auprc_ranking = (
    pd.DataFrame(auprc_rows)
    .groupby("tool")["median_AUPRC"]
    .mean()
    .reset_index()
    .rename(columns={"median_AUPRC": "overall_median_AUPRC"})
    .sort_values("overall_median_AUPRC", ascending=False)
    .reset_index(drop=True)
)
auprc_ranking["rank"] = auprc_ranking.index + 1

print("\n=== Tool Ranking (Cell-wise AUPRC) ===")
print(auprc_ranking.to_string(index=False))

# 8) PCC-AUPRC 분석 함수 정의

def run_cell_metric_scatter_facets(sim_name, sim_data, gt_log1p_adata,
                                   x_metric="AUPRC", y_metric=" PCC",
                                   n_sample_cells=500, ncols=4):
    """
    Cell-wise 2D scatter: x=AUPRC, y=SCC (또는 PCC)
    한 점 = cell 하나
    RAW 대비 두 지표 모두 개선된 cell 비율을 title에 표시
    """
    print(f"\n--- Cell-wise {x_metric} vs {y_metric} Scatter: {sim_name} ---")
    
    X_gt = gt_log1p_adata.X
    if hasattr(X_gt, "toarray"):
        X_gt = X_gt.toarray()
    gt_log1p = np.asarray(X_gt, dtype=np.float64)
    gt_genes  = gt_log1p_adata.var.index.tolist()
    groups    = np.asarray(gt_log1p_adata.obs["Group"].astype(str))
    
    # 공통 gene 교집합
    common_genes = set(gt_genes)
    for adata in sim_data.values():
        common_genes &= set(adata.var.index.tolist())
    common_genes = [g for g in gt_genes if g in common_genes]
    gt_col_idx   = [gt_genes.index(g) for g in common_genes]
    gt_sub       = gt_log1p[:, gt_col_idx]
    print(f"Common genes: {len(common_genes)}")
    
    # GT 이진화 (AUPRC용)
    gt_binary = (gt_sub > 0).astype(np.int8)
    
    # cell sampling (모든 tool 동일)
    np.random.seed(42)
    n_cells = gt_sub.shape[0]
    if n_cells > n_sample_cells:
        sampled_idx = np.sort(np.random.choice(n_cells, n_sample_cells, replace=False))
    else:
        sampled_idx = np.arange(n_cells)
    
    gt_val_sub = gt_sub[sampled_idx, :]       # (n_sample, n_genes) - SCC/PCC용
    gt_bin_sub = gt_binary[sampled_idx, :]    # (n_sample, n_genes) - AUPRC용
    print(f"Sampled cells: {len(sampled_idx)}")
    
    # 지표 계산 함수
    def cell_metric(gt_val, gt_bin, imp_vec, which):
        imp_vec = np.nan_to_num(imp_vec, nan=0.0, posinf=0.0, neginf=0.0)
        if which == "AUPRC":
            if gt_bin.sum() == 0:
                return np.nan
            return average_precision_score(gt_bin, imp_vec)
        elif which == "SCC":
            if np.all(gt_val == gt_val[0]) or np.all(imp_vec == imp_vec[0]):
                return np.nan
            return spearmanr(gt_val, imp_vec).statistic
        elif which == "PCC":
            if np.std(gt_val) == 0 or np.std(imp_vec) == 0:
                return np.nan
            return np.corrcoef(gt_val, imp_vec)[0, 1]
        raise ValueError(f"Unknown metric: {which}")
    
    # tool별 계산
    tool_points = {}
    x_all, y_all = [], []
    
    for tool_name, adata in sim_data.items():
        imp_genes     = adata.var.index.tolist()
        imp_gene_dict = {g: i for i, g in enumerate(imp_genes)}
        imp_col_idx   = [imp_gene_dict[g] for g in common_genes]
    
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        imp_sub = np.asarray(X, dtype=np.float64)[sampled_idx, :][:, imp_col_idx]
    
        xs, ys = [], []
        for i in range(len(sampled_idx)):
            xs.append(cell_metric(gt_val_sub[i], gt_bin_sub[i], imp_sub[i], x_metric))
            ys.append(cell_metric(gt_val_sub[i], gt_bin_sub[i], imp_sub[i], y_metric))
    
        pts = pd.DataFrame({
            "cell_idx": sampled_idx,
            x_metric: xs,
            y_metric: ys
        }).dropna()
        tool_points[tool_name] = pts
    
        x_all.extend(pts[x_metric].tolist())
        y_all.extend(pts[y_metric].tolist())
    
        print(f"  {tool_name}: median {x_metric}={np.nanmedian(xs):.3f}, "
              f"median {y_metric}={np.nanmedian(ys):.3f}")
    
        del imp_sub
        gc.collect()
    
    # 축 범위 공유
    pad = lambda v: 0.02 * (np.nanmax(v) - np.nanmin(v) + 1e-9)
    x_lim = (np.nanmin(x_all) - pad(x_all), np.nanmax(x_all) + pad(x_all))
    y_lim = (np.nanmin(y_all) - pad(y_all), np.nanmax(y_all) + pad(y_all))
    
    # RAW baseline
    if "RAW" not in tool_points or len(tool_points["RAW"]) == 0:
        print("[Warn] RAW not found.")
        return
    raw_df    = tool_points["RAW"][["cell_idx", x_metric, y_metric]].copy()
    raw_x_med = raw_df[x_metric].median()
    raw_y_med = raw_df[y_metric].median()
    
    # subplot
    tools = list(sim_data.keys())
    nrows = int(np.ceil(len(tools) / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(4.2 * ncols, 3.8 * nrows),
                              sharex=True, sharey=True)
    axes = np.array(axes).flatten()
    
    for ax, tool_name in zip(axes, tools):
        pts = tool_points.get(tool_name, pd.DataFrame())
    
        if len(pts):
            ax.scatter(pts[x_metric], pts[y_metric],
                       s=8, alpha=0.5, rasterized=True)
    
        # RAW median 기준선 (검은 점선)
        ax.axvline(raw_x_med, color="black", linestyle="--", linewidth=1.0, alpha=0.8)
        ax.axhline(raw_y_med, color="black", linestyle="--", linewidth=1.0, alpha=0.8)
    
        # tool median (빨간 점선)
        if len(pts):
            ax.axvline(pts[x_metric].median(), color="red", linestyle="--", linewidth=1.0)
            ax.axhline(pts[y_metric].median(), color="red", linestyle="--", linewidth=1.0)
    
            # per-cell improvement: RAW보다 두 지표 모두 높은 cell 비율
            merged = pts.merge(
                raw_df.rename(columns={x_metric: f"{x_metric}_RAW",
                                       y_metric: f"{y_metric}_RAW"}),
                on="cell_idx", how="inner"
            ).dropna()
    
            if len(merged):
                improved = (
                    (merged[x_metric] > merged[f"{x_metric}_RAW"]) &
                    (merged[y_metric] > merged[f"{y_metric}_RAW"])
                ).mean() * 100
                ax.set_title(f"{tool_name}\nImproved: {improved:.1f}% (n={len(merged)})",
                             fontsize=8)
            else:
                ax.set_title(f"{tool_name}\nImproved: n/a", fontsize=8)
        else:
            ax.set_title(f"{tool_name}\nno points", fontsize=8)
    
        ax.set_xlim(*x_lim)
        ax.set_ylim(*y_lim)
        ax.grid(True, linestyle="--", alpha=0.3)
    
    for j in range(len(tools), len(axes)):
        axes[j].axis("off")
    
    fig.suptitle(
        f"{sim_name}: Cell-wise {x_metric} vs {y_metric} (n={len(sampled_idx)} cells)\n"
        f"black dashed = RAW median / red dashed = tool median",
        fontsize=11
    )
    fig.supxlabel(x_metric, fontsize=10)
    fig.supylabel(y_metric, fontsize=10, x=0.01)
    plt.tight_layout(rect=[0.03, 0.0, 1, 0.95])
    
    out_fn = os.path.join(SAVE_DIR, f"{sim_name}_Cell_{x_metric}_vs_{y_metric}_Scatter.png")
    fig.savefig(out_fn, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out_fn}")
    
    del gt_val_sub, gt_bin_sub
    gc.collect()

# 실행
for sim_name, sim_data in sims.items():
    run_cell_metric_scatter_facets(
        sim_name       = sim_name,
        sim_data       = sim_data,
        gt_log1p_adata = gt,
        x_metric       = "AUPRC",
        y_metric       = "PCC",
        n_sample_cells = 500,
        ncols          = 4
    )

# 9) Cell-wise logRMSE 분석 함수 정의
def run_cell_logrmse_boxplot(sim_name, sim_data, gt_log1p_adata, n_sample_cells=500):
    """
    Cell-wise logRMSE:
    - GT log1p vs Imputed log1p 간 RMSE를 cell마다 계산
    - 낮을수록 좋음
    - RAW baseline (검은 점선)
    """
    print(f"\n--- Cell-wise logRMSE: {sim_name} ---")
    
    X_gt = gt_log1p_adata.X
    if hasattr(X_gt, "toarray"):
        X_gt = X_gt.toarray()
    gt_log1p = np.asarray(X_gt, dtype=np.float64)
    gt_genes  = gt_log1p_adata.var.index.tolist()
    
    # 공통 gene 교집합
    common_genes = set(gt_genes)
    for adata in sim_data.values():
        common_genes &= set(adata.var.index.tolist())
    common_genes = [g for g in gt_genes if g in common_genes]
    gt_col_idx   = [gt_genes.index(g) for g in common_genes]
    gt_sub       = gt_log1p[:, gt_col_idx]
    print(f"Common genes: {len(common_genes)}")
    
    # cell sampling
    np.random.seed(42)
    n_cells = gt_sub.shape[0]
    if n_cells > n_sample_cells:
        sampled_idx = np.sort(np.random.choice(n_cells, n_sample_cells, replace=False))
    else:
        sampled_idx = np.arange(n_cells)
    gt_sub_sampled = gt_sub[sampled_idx, :]
    print(f"Sampled cells: {len(sampled_idx)}")
    
    def cell_logrmse(gt_vec, imp_vec):
        imp_vec = np.nan_to_num(imp_vec, nan=0.0, posinf=0.0, neginf=0.0)
        return np.sqrt(np.mean((gt_vec - imp_vec) ** 2))
    
    rmse_records = []
    
    for tool_name, adata in sim_data.items():
        imp_genes     = adata.var.index.tolist()
        imp_gene_dict = {g: i for i, g in enumerate(imp_genes)}
        imp_col_idx   = [imp_gene_dict[g] for g in common_genes]
    
        X = adata.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        imp_sub = np.asarray(X, dtype=np.float64)[sampled_idx, :][:, imp_col_idx]
    
        scores = np.array([
            cell_logrmse(gt_sub_sampled[i], imp_sub[i])
            for i in range(len(sampled_idx))
        ])
    
        rmse_records.append(pd.DataFrame({"logRMSE": scores, "Tool": tool_name}))
        print(f"  {tool_name}: median logRMSE = {np.median(scores):.4f}")
    
        del imp_sub
        gc.collect()
    
    final_df = pd.concat(rmse_records, ignore_index=True)
    tool_order = list(sim_data.keys())
    
    raw_med = None
    if "RAW" in final_df["Tool"].values:
        raw_med = final_df.loc[final_df["Tool"] == "RAW", "logRMSE"].median()
    
    fig, ax = plt.subplots(figsize=(max(6, len(tool_order) * 0.9), 5.8))
    sns.boxplot(
        data=final_df,
        x="Tool", y="logRMSE",
        order=tool_order,
        hue="Tool", hue_order=tool_order,
        dodge=False, palette="tab20",
        showfliers=False, ax=ax
    )
    if raw_med is not None:
        ax.axhline(raw_med, color="black", linestyle="--",
                   linewidth=1, alpha=0.8, label=f"RAW median ({raw_med:.3f})")
        ax.legend(fontsize=8)
    
    ax.set_title(f"Cell-wise logRMSE vs GT - {sim_name}  (lower is better)")
    ax.set_ylabel("logRMSE (log1p scale)")
    ax.set_xlabel("")
    ax.set_xticks(range(len(tool_order)))
    ax.set_xticklabels(tool_order, rotation=45, ha="right")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    
    out_fn = os.path.join(SAVE_DIR, f"{sim_name}_Cell_logRMSE_Boxplot.png")
    fig.savefig(out_fn, dpi=300)
    plt.close(fig)
    print(f"[saved] {out_fn}")
    
    del gt_sub_sampled, rmse_records
    gc.collect()
    
    return final_df


# 실행
logrmse_results = {}
for sim_name, sim_data in sims.items():
    rmse_df = run_cell_logrmse_boxplot(
        sim_name       = sim_name,
        sim_data       = sim_data,
        gt_log1p_adata = gt,
        n_sample_cells = 500
    )
    logrmse_results[sim_name] = rmse_df

# 랭킹 (낮을수록 좋으므로 ascending=True)
rmse_rows = []
for sim_name, rmse_df in logrmse_results.items():
    for tool, val in rmse_df.groupby("Tool")["logRMSE"].median().items():
        rmse_rows.append({"sim": sim_name, "tool": tool, "median_logRMSE": val})

rmse_ranking = (
    pd.DataFrame(rmse_rows)
    .groupby("tool")["median_logRMSE"]
    .mean()
    .reset_index()
    .rename(columns={"median_logRMSE": "overall_median_logRMSE"})
    .sort_values("overall_median_logRMSE", ascending=True)   # 낮을수록 좋음
    .reset_index(drop=True)
)
rmse_ranking["rank"] = rmse_ranking.index + 1

print("\n=== Tool Ranking (Cell-wise logRMSE, lower is better) ===")
print(rmse_ranking.to_string(index=False))

# ===================================================================
# Export results to CSV for summary_visualization.py
# ===================================================================
sim_results = pd.DataFrame({
    "tool": pcc_ranking["tool"],
}).merge(
    pcc_ranking[["tool", "overall_median_PCC"]].rename(columns={"overall_median_PCC": "cell_PCC"}), on="tool"
).merge(
    ranking[["tool", "mean_SCC"]].rename(columns={"mean_SCC": "Pseudobulk_AUPRC"}), on="tool"
).merge(
    auprc_ranking[["tool", "overall_median_AUPRC"]].rename(columns={"overall_median_AUPRC": "AUPRC (cell-level)"}), on="tool"
).merge(
    rmse_ranking[["tool", "overall_median_logRMSE"]].rename(columns={"overall_median_logRMSE": "logRMSE"}), on="tool"
).merge(
    overall_mean[["tool", "overall_mean_SCC"]].rename(columns={"overall_mean_SCC": "SCC_FC"}), on="tool"
)

sim_results.to_csv(os.path.join(SAVE_DIR, "../results/eval_simulation.csv"), index=False)
print("\n✅ Saved: results/eval_simulation.csv")
