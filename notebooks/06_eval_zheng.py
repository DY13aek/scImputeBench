# Set the base directory for datasets
# Change this path to point to your local data directory
BASE_DIR = "../data"

# 필요한 라이브러리 임포트
import os
import scanpy as sc
import numpy as np
import pandas as pd

os.chdir(f"{BASE_DIR}/Zheng/Raw")
raw = sc.read_h5ad("zheng.h5ad")
sc.pp.normalize_total(raw, target_sum=1e4)
sc.pp.log1p(raw)
cell_idx = raw.obs.index.tolist()
gene_idx = raw.var.index.tolist()    

# imputed datasets
os.chdir(f"{BASE_DIR}/Zheng/imputed")

# === 1. DCA (latent.tsv) 말고 일단 mean으로 가져옴: Done ===
#dca_df = pd.read_csv("./DCA/mean.tsv", sep="\t", header=0, index_col=0, engine="pyarrow").T
#sc.AnnData(dca_df.values.astype('float32')).write_h5ad("./DCA/dca.h5ad")# 2. X만 떼서 저장 (인덱스/메타데이터 문제 원천 차단)
dca = sc.read_h5ad("./DCA/dca.h5ad") # 3. 다시 로드
dca.obs, dca.var = raw.obs.copy(), raw.var.copy() # 4. raw의 정보 그대로 이식 (Shape과 순서가 같으므로)
dca.obs_names, dca.var_names = raw.obs_names, raw.var_names
#dca.write_h5ad("./DCA/dca.h5ad")# 최종 저장 및 불러오기
print("Shape:", dca.shape)
sc.pp.normalize_total(dca, target_sum=1e4)
sc.pp.log1p(dca)
sc.pp.highly_variable_genes(dca, n_top_genes=2000)

# === 2. DrImpute: Done ===
#drimpute_df = pd.read_csv("./DrImpute/zheng_DrImpute.csv", sep=",", index_col=0, engine="pyarrow").T
#temp_adata = sc.AnnData(drimpute_df.values.astype('float32'))
#saved_var_names = drimpute_df.columns # var(유전자 이름) 따로 빼두기
#sc.AnnData(temp_adata.X).write_h5ad("./DrImpute/drimpute.h5ad") # X만 저장 후 다시 불러오기 (메타데이터 에러 원천 차단)
drimpute = sc.read_h5ad("./DrImpute/drimpute.h5ad")
drimpute.obs, drimpute.var = raw.obs.copy(), raw.var.copy() # obs는 raw에서 가져와서 붙이기
drimpute.obs_names, drimpute.var_names = cell_idx, gene_idx  # 순서 같으니 이름 덮어씌우기
#drimpute.write_h5ad("./DrImpute/drimpute.h5ad") # 최종 저장
# 확인
print("Shape:", drimpute.shape)

# === 3. MAGIC: Done ===
magic = sc.read_h5ad("./MAGIC/zheng_MAGIC.h5ad")
print("magic shape:", magic.shape)

# === 4. SAVER: Done ===
#saver_df = pd.read_csv("./SAVER/zheng_saver_imputed.csv", sep=",", index_col=0, engine="pyarrow").T
#sc.AnnData(saver_df.values.astype('float32')).write_h5ad("./SAVER/saver.h5ad")
saver = sc.read_h5ad("./SAVER/saver.h5ad")
saver.obs = raw.obs.copy() 
saver.var = raw.var.copy()
saver.obs_names, saver.var_names = cell_idx, gene_idx
#saver.write_h5ad("./SAVER/saver.h5ad")
print("saver shape:", saver.shape)
sc.pp.normalize_total(saver, target_sum=1e4)
sc.pp.log1p(saver)
sc.pp.highly_variable_genes(saver, n_top_genes=2000)

# === 5. scIDPMs: Done ===
#scidpms_df = pd.read_csv("./scIDPMs/imputed.csv", sep=",", index_col=None, header=None)
#sc.AnnData(scidpms_df.values.astype('float32')).write_h5ad("./scIDPMs/scidpms.h5ad")
scidpms = sc.read_h5ad("./scIDPMs/scidpms.h5ad")
scidpms.obs, scidpms.var = raw.obs.copy(), raw.var.copy()
scidpms.obs_names, scidpms.var_names = cell_idx, gene_idx
print("scidpms shape:", scidpms.shape)
#scidpms.write_h5ad("./scIDPMs/scidpms.h5ad")
sc.pp.normalize_total(scidpms, target_sum=1e4)
sc.pp.log1p(scidpms)
sc.pp.highly_variable_genes(scidpms, n_top_genes=2000)

# === 5. scIDPMs_norm: Done ===
#scidpms_norm_df = pd.read_csv("./scIDPMs_norm/imputed.csv", sep=",", index_col=None, header=None)
#sc.AnnData(scidpms_norm_df.values.astype('float32')).write_h5ad("./scIDPMs_norm/scidpms_norm.h5ad")
scidpms_norm = sc.read_h5ad("./scIDPMs_norm/scidpms_norm.h5ad")
scidpms_norm.obs, scidpms_norm.var = raw.obs.copy(), raw.var.copy()
scidpms_norm.obs_names, scidpms_norm.var_names = cell_idx, gene_idx
print("scidpms_norm shape:", scidpms_norm.shape)
#scidpms_norm.write_h5ad("./scIDPMs_norm/scidpms_norm.h5ad")

# === 6. scIGANs: Done ===
#scigans_df = pd.read_csv("./scIGANs/scIGANs_zheng.txt", sep="\t", index_col=0, engine="pyarrow").T
#sc.AnnData(scigans_df.values.astype('float32')).write_h5ad("./scIGANs/scigans.h5ad")
scigans = sc.read_h5ad("./scIGANs/scigans.h5ad")
scigans.obs = raw.obs.copy()
scigans.var = raw.var.copy()
scigans.obs_names, scigans.var_names = cell_idx, gene_idx
#scigans.write_h5ad("./scIGANs/scigans.h5ad")
print("scigans shape:", scigans.shape)
sc.pp.normalize_total(scigans, target_sum=1e4)
sc.pp.log1p(scigans)
sc.pp.highly_variable_genes(scigans, n_top_genes=2000)

# === 8. scMultiGAN ===
#scmultigan_df = pd.read_csv("./scMultiGAN/Restored/zheng_scMultiGAN.tsv", sep="\t", header=None, index_col=None, engine="pyarrow").T
#sc.AnnData(scmultigan_df.values.astype('float32')).write_h5ad("./scMultiGAN/scmultigan.h5ad")
scmultigan = sc.read_h5ad("./scMultiGAN/scmultigan.h5ad")
scmultigan.obs, scmultigan.var = raw.obs.copy(), raw.var.copy()
scmultigan.obs_names, scmultigan.var_names = cell_idx, gene_idx
#scmultigan.write_h5ad("./scMultiGAN/scmultigan.h5ad")
print("scmultigan shape:", scmultigan.shape)
sc.pp.normalize_total(scmultigan, target_sum=1e4)
sc.pp.log1p(scmultigan)
sc.pp.highly_variable_genes(scmultigan, n_top_genes=2000)

# === 9. scSTD: Done ===
#scstd_df = pd.read_csv("./scSTD/Result/zheng_imputed.txt", sep=",", index_col=None, header=None, engine="pyarrow")
#sc.AnnData(scstd_df.values.astype('float32')).write_h5ad("./scSTD/scstd.h5ad")
scstd = sc.read_h5ad("./scSTD/scstd.h5ad")
scstd.obs, scstd.var = raw.obs.copy(), raw.var.copy()
scstd.obs_names, scstd.var_names = raw.obs_names, raw.var_names
#scstd.write_h5ad("./scSTD/scstd.h5ad")
print("scstd shape:", scstd.shape)
sc.pp.normalize_total(scstd, target_sum=1e4)
sc.pp.log1p(scstd)
sc.pp.highly_variable_genes(scstd, n_top_genes=2000)

# === 10. scVI: Done ===
scvi = sc.read_h5ad("./scVI/scvi_imputed_zheng.h5ad")
scvi.obs_names, scvi.var_names = cell_idx, gene_idx
print("scvi shape:", scvi.shape)

# === 11. ALRA: Done ===
#alra_df = pd.read_csv("./ALRA/zheng_ALRA_imputed.csv", index_col=0, engine="pyarrow")
#sc.AnnData(alra_df.values.astype('float32')).write_h5ad("./ALRA/alra.h5ad")
alra = sc.read_h5ad("./ALRA/alra.h5ad")
alra.obs, alra.var = raw.obs.copy(), raw.var.copy()
alra.obs_names, alra.var_names = cell_idx, gene_idx
#alra.write_h5ad("./ALRA/alra.h5ad")
print("alra shape:", alra.shape)


zheng_data = {
    "RAW":          raw,
    "ALRA":         alra,
    "DCA":          dca,
    "DrImpute":     drimpute,
    "MAGIC":        magic,
    "SAVER":        saver,
    "scIDPMs":      scidpms,
    "scIGANs":      scigans,
    "scMultiGAN":   scmultigan,
    "scSTD":        scstd,
    "scVI":         scvi,
}

SAVE_DIR_ZHENG = f"{BASE_DIR}/Zheng/figures"
os.makedirs(SAVE_DIR_ZHENG, exist_ok=True)

#── 1. Marker gene specificity (Fine grain) ─────────────────────────────────────

# ================================================================
# 0. 설정
# ================================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import LabelEncoder

SAVE_DIR = f"{BASE_DIR}/Zheng/figures"
os.makedirs(SAVE_DIR, exist_ok=True)

CELLTYPE_COL = "cell_type"  # raw.obs.columns 확인 후 수정

# cell type 그룹핑 두 가지
cell_type_broad = {
    "CD14+ Monocyte":               "Monocyte",
    "CD19+ B":                      "B cell",
    "CD4+/CD25 T Reg":              "CD4+ T",
    "CD4+/CD45RA+/CD25- Naive T":   "CD4+ T",
    "CD4+/CD45RO+ Memory":          "CD4+ T",
    "CD56+ NK":                     "NK cell",
    "CD8+ Cytotoxic T":             "CD8+ T",
    "CD8+/CD45RA+ Naive Cytotoxic": "CD8+ T",
    "Dendritic":                    "Dendritic"
}

# fine: 원본 그대로 사용
cell_type_fine = {ct: ct for ct in cell_type_broad.keys()}

# marker gene도 fine grain으로 확장
marker_genes_broad = {
    "Monocyte":  ["CD14", "LYZ", "CST3"],
    "B cell":    ["CD19", "MS4A1", "CD79A"],
    "CD4+ T":    ["CD4", "IL7R"],
    "CD8+ T":    ["CD8A", "CD8B"],
    "NK cell":   ["GNLY", "NKG7", "NCAM1"],
    "Dendritic": ["FCER1A", "CST7"]
}

marker_genes_fine = {
    "CD14+ Monocyte":               ["CD14", "LYZ", "CST3"],
    "CD19+ B":                      ["CD19", "MS4A1", "CD79A"],
    "CD4+/CD25 T Reg":              ["CD4", "FOXP3", "IL2RA"],
    "CD4+/CD45RA+/CD25- Naive T":   ["CD4", "CCR7", "SELL"],
    "CD4+/CD45RO+ Memory":          ["CD4", "IL7R", "S100A4"],
    "CD56+ NK":                     ["GNLY", "NKG7", "NCAM1"],
    "CD8+ Cytotoxic T":             ["CD8A", "GZMB", "PRF1"],
    "CD8+/CD45RA+ Naive Cytotoxic": ["CD8A", "CCR7", "SELL"],
    "Dendritic":                    ["FCER1A", "CST7", "HLA-DRA"]
}

data_dict = {
    "RAW":        raw,
    "ALRA":       alra,
    "DCA":        dca,
    "DrImpute":   drimpute,
    "MAGIC":      magic,
    "SAVER":      saver,
    "scIDPMs":    scidpms,
    "scIGANs":    scigans,
    "scMultiGAN": scmultigan,
    "scSTD":      scstd,
    "scVI":       scvi,
}

# ================================================================
# 1. Marker Gene Specificity 함수
# ================================================================
def get_marker_specificity(adata, ct_map, marker_dict, tool_name):
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)
    
    genes     = adata.var.index.tolist()
    gene_dict = {g: i for i, g in enumerate(genes)}
    
    broad_labels = np.array([
        ct_map.get(ct, ct)
        for ct in adata.obs[CELLTYPE_COL].astype(str)
    ])
    
    rows = []
    for cell_type, markers in marker_dict.items():
        target_idx = np.where(broad_labels == cell_type)[0]
        other_idx  = np.where(broad_labels != cell_type)[0]
    
        for gene in markers:
            if gene not in gene_dict:
                print(f"  [skip] {gene} not in {tool_name}")
                continue
            vec        = X[:, gene_dict[gene]]
            med_target = np.median(vec[target_idx])
            med_other  = np.median(vec[other_idx])
            denom      = med_target + med_other
            specificity = med_target / denom if denom > 0 else np.nan
    
            rows.append({
                "tool": tool_name, "cell_type": cell_type,
                "gene": gene, "med_target": med_target,
                "med_other": med_other, "specificity": specificity
            })
    return pd.DataFrame(rows)


def run_marker_specificity(marker_dict, ct_map, label="broad"):
    """Marker specificity 계산 → heatmap + boxplot 저장"""
    all_rows = []
    for tool_name, adata in data_dict.items():
        print(f"  [{label}] Processing {tool_name}...")
        df = get_marker_specificity(adata, ct_map, marker_dict, tool_name)
        nan_genes = df[df["specificity"].isna()]["gene"].tolist()
        if nan_genes:
            print(f"    [NaN→0] {tool_name}: {nan_genes}")
        all_rows.append(df)
    
    spec_df = pd.concat(all_rows, ignore_index=True)
    
    # 랭킹 (NaN=0 페널티)
    summary = (
        spec_df.fillna({"specificity": 0.0})
        .groupby("tool")["specificity"].mean()
        .reset_index()
        .rename(columns={"specificity": "mean_specificity"})
        .sort_values("mean_specificity", ascending=False)
        .reset_index(drop=True)
    )
    summary["rank"] = summary.index + 1
    summary["nan_count"] = summary["tool"].map(
        spec_df.groupby("tool")["specificity"].apply(lambda x: x.isna().sum())
    )
    print(f"\n=== Marker Specificity Ranking [{label}] ===")
    print(summary.to_string(index=False))
    
    # Heatmap
    spec_df["ct_gene"] = spec_df["cell_type"] + "  |  " + spec_df["gene"]
    ct_order = [f"{ct}  |  {g}" for ct, gs in marker_dict.items() for g in gs]
    
    pivot = (
        spec_df.fillna({"specificity": 0.0})
        .pivot_table(index="ct_gene", columns="tool",
                     values="specificity", aggfunc="mean")
    )
    pivot = pivot.loc[[c for c in ct_order if c in pivot.index]]
    pivot = pivot[summary["tool"].tolist()]
    
    fig, ax = plt.subplots(figsize=(max(10, len(pivot.columns) * 1.0),
                                    max(6,  len(pivot) * 0.55 + 2)))
    sns.heatmap(pivot, annot=True, fmt=".2f",
                cmap="YlOrRd", vmin=0.0, vmax=1.0,
                linewidths=0.3, linecolor="gray", ax=ax)
    
    # cell type 구분선
    cumsum = 0
    for ct, gs in list(marker_dict.items())[:-1]:
        cumsum += len(gs)
        ax.axhline(cumsum, color="black", linewidth=1.5)
    
    ax.set_title(f"Marker Gene Specificity [{label}]  |  NaN→0 (dropout penalty)",
                 fontsize=10)
    ax.set_ylabel("Cell Type  |  Marker Gene")
    ax.set_xlabel("Tool")
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)
    plt.tight_layout()
    fig.savefig(os.path.join(SAVE_DIR, f"Marker_Specificity_Heatmap_{label}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] Marker_Specificity_Heatmap_{label}.png")
    
    # Boxplot
    tool_order = summary["tool"].tolist()
    fig, ax = plt.subplots(figsize=(max(8, len(tool_order) * 0.9), 5.5))
    sns.boxplot(
        data=spec_df.fillna({"specificity": 0.0}),
        x="tool", y="specificity",
        order=tool_order, hue="tool", hue_order=tool_order,
        dodge=False, palette="tab20", showfliers=False, ax=ax
    )
    if "RAW" in spec_df["tool"].values:
        raw_med = spec_df.fillna({"specificity": 0.0}).loc[
            spec_df["tool"] == "RAW", "specificity"].median()
        ax.axhline(raw_med, color="black", linestyle="--", linewidth=1,
                   alpha=0.8, label=f"RAW median ({raw_med:.3f})")
        ax.legend(fontsize=8)
    ax.set_title(f"Marker Gene Specificity [{label}]  (NaN→0)")
    ax.set_ylabel("Specificity (higher is better)")
    ax.set_xlabel("")
    ax.set_xticks(range(len(tool_order)))
    ax.set_xticklabels(tool_order, rotation=45, ha="right")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    fig.savefig(os.path.join(SAVE_DIR, f"Marker_Specificity_Boxplot_{label}.png"),
                dpi=300)
    plt.close(fig)
    print(f"[saved] Marker_Specificity_Boxplot_{label}.png")
    return spec_df, summary


# ================================================================
# 2. Silhouette Score 함수
# ================================================================
def run_silhouette(ct_map, label="broad", n_hvg=2000):
    """
    각 tool의 log1p norm + HVG 기반 PCA embedding에서 silhouette score 계산
    """
    print(f"\n=== Silhouette Score [{label}] ===")
    rows = []
    
    for tool_name, adata in data_dict.items():
        print(f"  Processing {tool_name}...")
        adata_tmp = adata.copy()
    
        # log1p norm이 안 된 경우 처리 (RAW)
        if tool_name == "RAW":
            sc.pp.normalize_total(adata_tmp, target_sum=1e4)
            sc.pp.log1p(adata_tmp)
    
        # broad/fine label 추가
        adata_tmp.obs["label"] = [
            ct_map.get(ct, ct)
            for ct in adata_tmp.obs[CELLTYPE_COL].astype(str)
        ]
    
        # HVG → PCA
        sc.pp.highly_variable_genes(adata_tmp, n_top_genes=n_hvg)
        adata_tmp = adata_tmp[:, adata_tmp.var["highly_variable"]]
        sc.pp.scale(adata_tmp, max_value=10)
        sc.tl.pca(adata_tmp, n_comps=50)
    
        X_pca = adata_tmp.obsm["X_pca"]
        labels = adata_tmp.obs["label"].values
        le = LabelEncoder()
        labels_enc = le.fit_transform(labels)
    
        score = silhouette_score(X_pca, labels_enc, sample_size=3000,
                                 random_state=42)
        rows.append({"tool": tool_name, "silhouette": score})
        print(f"    {tool_name}: silhouette = {score:.4f}")
    
    sil_df = pd.DataFrame(rows).sort_values("silhouette", ascending=False).reset_index(drop=True)
    sil_df["rank"] = sil_df.index + 1
    
    print(sil_df.to_string(index=False))
    
    # Barplot
    fig, ax = plt.subplots(figsize=(max(8, len(sil_df) * 0.9), 5))
    sns.barplot(data=sil_df, x="tool", y="silhouette",
                order=sil_df["tool"], palette="tab20", ax=ax)
    if "RAW" in sil_df["tool"].values:
        raw_val = sil_df.loc[sil_df["tool"] == "RAW", "silhouette"].values[0]
        ax.axhline(raw_val, color="black", linestyle="--", linewidth=1,
                   alpha=0.8, label=f"RAW ({raw_val:.3f})")
        ax.legend(fontsize=8)
    ax.set_title(f"Silhouette Score [{label}]  (higher is better)")
    ax.set_ylabel("Silhouette Score")
    ax.set_xlabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    fig.savefig(os.path.join(SAVE_DIR, f"Silhouette_{label}.png"), dpi=300)
    plt.close(fig)
    print(f"[saved] Silhouette_{label}.png")
    
    return sil_df


# ================================================================
# 3. UMAP 함수
# ================================================================
def run_umap(ct_map, label="broad", n_hvg=2000):
    """
    tool별 UMAP subplot - 한 figure에 모든 tool
    cell type은 ct_map 기준으로 색상 구분
    """
    print(f"\n=== UMAP [{label}] ===")
    
    tools     = list(data_dict.keys())
    n_tools   = len(tools)
    ncols     = 4
    nrows     = int(np.ceil(n_tools / ncols))
    all_labels = sorted(set(ct_map.values()))
    palette    = dict(zip(all_labels,
                          sns.color_palette("tab10", n_colors=len(all_labels))))
    
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(ncols * 4.5, nrows * 4.0))
    axes = np.array(axes).flatten()
    
    for ax, tool_name in zip(axes, tools):
        print(f"  Processing {tool_name}...")
        adata_tmp = data_dict[tool_name].copy()
    
        adata_tmp.obs["label"] = [
            ct_map.get(ct, ct)
            for ct in adata_tmp.obs[CELLTYPE_COL].astype(str)
        ]
    
        sc.pp.highly_variable_genes(adata_tmp, n_top_genes=n_hvg)
        adata_tmp = adata_tmp[:, adata_tmp.var["highly_variable"]]
        sc.pp.scale(adata_tmp, max_value=10)
        sc.tl.pca(adata_tmp, n_comps=50)
        sc.pp.neighbors(adata_tmp, n_neighbors=15, n_pcs=50)
        sc.tl.umap(adata_tmp)
    
        umap_coords = adata_tmp.obsm["X_umap"]
        labels      = adata_tmp.obs["label"].values
    
        for ct in all_labels:
            idx = np.where(labels == ct)[0]
            ax.scatter(umap_coords[idx, 0], umap_coords[idx, 1],
                       c=[palette[ct]], s=2, alpha=0.5,
                       label=ct, rasterized=True)
    
        ax.set_title(tool_name, fontsize=10)
        ax.set_xlabel("UMAP1", fontsize=7)
        ax.set_ylabel("UMAP2", fontsize=7)
        ax.tick_params(labelsize=6)
    
    # 공통 legend
    handles = [plt.Line2D([0], [0], marker='o', color='w',
                           markerfacecolor=palette[ct], markersize=7, label=ct)
               for ct in all_labels]
    fig.legend(handles=handles, loc="lower center",
               ncol=min(len(all_labels), 5),
               bbox_to_anchor=(0.5, -0.02), fontsize=8,
               title="Cell Type", title_fontsize=9)
    
    # 남는 subplot 숨기기
    for j in range(n_tools, len(axes)):
        axes[j].axis("off")
    
    fig.suptitle(f"UMAP [{label}]  —  tool별 비교", fontsize=13, y=1.01)
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    
    out_fn = os.path.join(SAVE_DIR, f"UMAP_{label}.png")
    fig.savefig(out_fn, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] UMAP_{label}.png")


# ================================================================
# 실행
# ================================================================
print("=" * 60)
print("1. Marker Gene Specificity")
print("=" * 60)
spec_broad, summary_broad = run_marker_specificity(
    marker_genes_broad, cell_type_broad, label="broad")
spec_fine, summary_fine = run_marker_specificity(
    marker_genes_fine, cell_type_fine, label="fine")

print("=" * 60)
print("2. Silhouette Score")
print("=" * 60)
sil_broad = run_silhouette(cell_type_broad, label="broad")
sil_fine  = run_silhouette(cell_type_fine,  label="fine")

print("=" * 60)
print("3. UMAP")
print("=" * 60)
run_umap(cell_type_broad, label="broad")
run_umap(cell_type_fine,  label="fine")

# ── 2. cell_cell_correlation (Heatmap) ─────────────────────────────────────

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import os

SEED        = 42
SAMPLE_FRAC = 0.1
SAVE_DIR    = f"{BASE_DIR}/Zheng/figures"
CELLTYPE_COL = "cell_type"  # raw.obs.columns 확인 후 수정
os.makedirs(SAVE_DIR, exist_ok=True)

data_dict = {
    "RAW":        raw,
    "ALRA":       alra,
    "DCA":        dca,
    "DrImpute":   drimpute,
    "MAGIC":      magic,
    "SAVER":      saver,
    "scIDPMs":    scidpms,
    "scIGANs":    scigans,
    "scMultiGAN": scmultigan,
    "scSTD":      scstd,
    "scVI":       scvi,
}
tools_ordered = list(data_dict.keys())

# broad/fine 라벨 정의
cell_type_broad = {
    "CD14+ Monocyte":               "Monocyte",
    "CD19+ B":                      "B cell",
    "CD4+/CD25 T Reg":              "CD4+ T",
    "CD4+/CD45RA+/CD25- Naive T":   "CD4+ T",
    "CD4+/CD45RO+ Memory":          "CD4+ T",
    "CD56+ NK":                     "NK cell",
    "CD8+ Cytotoxic T":             "CD8+ T",
    "CD8+/CD45RA+ Naive Cytotoxic": "CD8+ T",
    "Dendritic":                    "Dendritic"
}

broad_order = ["Monocyte", "B cell", "CD4+ T", "CD8+ T", "NK cell", "Dendritic"]
fine_order  = [
    "CD14+ Monocyte", "CD19+ B",
    "CD4+/CD25 T Reg", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory",
    "CD8+ Cytotoxic T", "CD8+/CD45RA+ Naive Cytotoxic",
    "CD56+ NK", "Dendritic"
]

# ================================================================
def run_cell_cell_correlation_plot(ct_map, order_list, file_suffix):
    """
    ct_map    : None이면 원본 label 그대로, dict이면 broad로 합쳐서 사용
    order_list: cell type 순서
    """
    print(f"\n--- Cell-Cell Correlation [{file_suffix}] ---")
    np.random.seed(SEED)
    
    # 샘플링: raw.obs 기준으로 인덱스 고정
    obs_labels = raw.obs[CELLTYPE_COL].astype(str)
    if ct_map is not None:
        obs_labels = obs_labels.map(ct_map)
    
    sampled_indices  = []
    cell_type_counts = {}
    
    for ct in order_list:
        indices = np.where(obs_labels == ct)[0]
        if len(indices) > 0:
            n_sample = max(int(len(indices) * SAMPLE_FRAC), min(len(indices), 10))
            selected = np.random.choice(indices, n_sample, replace=False)
            sampled_indices.extend(selected.tolist())
            cell_type_counts[ct] = len(selected)
        else:
            cell_type_counts[ct] = 0
    
    print(f"  Total sampled cells: {len(sampled_indices)}")
    
    palette      = sns.color_palette("husl", len(order_list))
    group_colors = {ct: c for ct, c in zip(order_list, palette)}
    
    n_cols = 4
    n_rows = int(np.ceil(len(tools_ordered) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(7 * n_cols, 7 * n_rows))
    axes = axes.flatten()
    last_im = None
    
    for i, tool_name in enumerate(tools_ordered):
        ax = axes[i]
        print(f"  -> {tool_name}...")
    
        adata = data_dict[tool_name]
        X = adata.X[sampled_indices, :]
        if hasattr(X, "toarray"):
            X = X.toarray()
        X = np.asarray(X, dtype=np.float64)
    
        corr_mat = np.corrcoef(X)
    
        last_im = sns.heatmap(
            corr_mat,
            ax=ax, cmap="RdBu_r", center=0,
            vmin=-0.2, vmax=1.0,
            cbar=False,
            xticklabels=False, yticklabels=False,
            square=True, rasterized=True
        )
    
        # cell type 색상 띠
        n_total       = len(sampled_indices)
        bar_thickness = n_total * 0.04
        current_idx   = 0
    
        for ct in order_list:
            n = cell_type_counts[ct]
            if n > 0:
                color = group_colors[ct]
                ax.add_patch(mpatches.Rectangle(
                    (current_idx, -bar_thickness), n, bar_thickness,
                    facecolor=color, edgecolor="none",
                    transform=ax.transData, clip_on=False
                ))
                ax.add_patch(mpatches.Rectangle(
                    (-bar_thickness, current_idx), bar_thickness, n,
                    facecolor=color, edgecolor="none",
                    transform=ax.transData, clip_on=False
                ))
                current_idx += n
    
        ax.set_title(tool_name, fontsize=20, fontweight="bold", pad=30)
    
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")
    
    # 범례
    legend_handles = [
        mpatches.Patch(color=c, label=ct)
        for ct, c in group_colors.items()
    ]
    fig.legend(handles=legend_handles, loc="center left",
               bbox_to_anchor=(0.91, 0.6),
               title="Cell Type", fontsize=14, title_fontsize=16,
               frameon=False)
    
    # 공통 컬러바
    cbar_ax = fig.add_axes([0.92, 0.12, 0.015, 0.25])
    cbar = fig.colorbar(
        last_im.get_children()[0], cax=cbar_ax, orientation="vertical"
    )
    cbar.set_label("Pearson r", fontsize=12)
    
    plt.suptitle(f"Cell-Cell Correlation [{file_suffix}]",
                 fontsize=28, y=0.98, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 0.90, 0.95])
    
    out_fn = os.path.join(SAVE_DIR, f"Cell_Cell_Correlation_{file_suffix}.png")
    fig.savefig(out_fn, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out_fn}")


# ================================================================
# 실행
# ================================================================
# Broad (합쳐진 cell type)
run_cell_cell_correlation_plot(
    ct_map      = cell_type_broad,
    order_list  = broad_order,
    file_suffix = "broad"
)

# Fine (세분화된 cell type)
run_cell_cell_correlation_plot(
    ct_map      = None,
    order_list  = fine_order,
    file_suffix = "fine"
)

# ── 3. gene-gene co-expression contrast (Barplot) ─────────────────────────────────────

def run_gene_gene_contrast(file_suffix="zheng"):
    """
    Gene-gene co-expression contrast:
    - intra-group corr: 같은 cell type marker 간 상관 평균
    - inter-group corr: 다른 cell type marker 간 상관 평균
    - contrast = intra - inter (높을수록 co-expression 구조 보존)
    broad / fine 두 가지 marker set으로 계산
    """
    print(f"\n--- Gene-Gene Co-expression Contrast: {file_suffix} ---")
    
    marker_sets = {
        "broad": {
            "Monocyte":  ["CD14", "LYZ", "CST3"],
            "B cell":    ["CD19", "MS4A1", "CD79A"],
            "CD4+ T":    ["CD4", "IL7R"],
            "CD8+ T":    ["CD8A", "CD8B"],
            "NK cell":   ["GNLY", "NKG7", "NCAM1"],
            "Dendritic": ["FCER1A", "CST7"]
        },
        "fine": {
            "CD14+ Monocyte":               ["CD14", "LYZ", "CST3", "FCGR3A"],
            "CD19+ B":                      ["CD19", "MS4A1", "CD79A", "CD79B"],
            "CD4+/CD25 T Reg":              ["CD4", "FOXP3", "IL2RA"],
            "CD4+/CD45RA+/CD25- Naive T":   ["CD4", "CCR7", "SELL"],
            "CD4+/CD45RO+ Memory":          ["CD4", "IL7R", "S100A4"],
            "CD56+ NK":                     ["GNLY", "NKG7", "NCAM1"],
            "CD8+ Cytotoxic T":             ["CD8A", "CD8B", "GZMB", "PRF1"],
            "CD8+/CD45RA+ Naive Cytotoxic": ["CD8A", "CD8B", "CCR7", "SELL"],
            "Dendritic":                    ["FCER1A", "CST7", "HLA-DRA"]
        }
    }
    
    all_summary = {}
    
    for level, marker_dict in marker_sets.items():
        print(f"\n  === {level} ===")
    
        all_markers = [g for genes in marker_dict.values() for g in genes]
    
        # intra/inter 마스크 사전 계산
        intra_mask = np.zeros((len(all_markers), len(all_markers)), dtype=bool)
        current = 0
        for ct, genes in marker_dict.items():
            n = len(genes)
            for ii in range(current, current + n):
                for jj in range(current, current + n):
                    if ii != jj:
                        intra_mask[ii, jj] = True
            current += n
        inter_mask = ~intra_mask & ~np.eye(len(all_markers), dtype=bool)
    
        rows = []
    
        for tool_name, adata in data_dict.items():
            avail_genes    = adata.var.index.tolist()
            valid_markers  = [g for g in all_markers if g in avail_genes]
            missing        = set(all_markers) - set(valid_markers)
            if missing:
                print(f"    [skip genes] {tool_name}: {missing}")
    
            valid_idx  = [all_markers.index(g) for g in valid_markers]
            intra_sub  = intra_mask[np.ix_(valid_idx, valid_idx)]
            inter_sub  = inter_mask[np.ix_(valid_idx, valid_idx)]
    
            X = adata[:, valid_markers].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            corr_mat = np.corrcoef(np.asarray(X, dtype=np.float64).T)
    
            intra_med = np.nanmedian(corr_mat[intra_sub])
            inter_med = np.nanmedian(corr_mat[inter_sub])
            contrast  = intra_med - inter_med
    
            rows.append({
                "tool":      tool_name,
                "intra_med": intra_med,
                "inter_med": inter_med,
                "contrast":  contrast,
            })
            print(f"    {tool_name}: intra={intra_med:.3f}  "
                  f"inter={inter_med:.3f}  contrast={contrast:.3f}")
    
        summary_df = (
            pd.DataFrame(rows)
            .sort_values("contrast", ascending=False)
            .reset_index(drop=True)
        )
        summary_df["rank"] = summary_df.index + 1
        all_summary[level]  = summary_df
    
        print(f"\n  === Ranking [{level}] ===")
        print(summary_df.to_string(index=False))
    
        # ── Figure 1: barplot (intra / inter / contrast) ──────
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        raw_row   = summary_df[summary_df["tool"] == "RAW"].iloc[0]
    
        for ax, col, title, color in zip(
            axes,
            ["intra_med", "inter_med", "contrast"],
            ["Intra-group Correlation\n(same cell type markers)",
             "Inter-group Correlation\n(different cell type markers)",
             "Contrast (intra - inter)\n(higher = better co-expression structure)"],
            ["tab:red", "tab:blue", "tab:green"]
        ):
            order = summary_df["tool"].tolist()
            sns.barplot(data=summary_df, x="tool", y=col,
                        order=order, color=color, alpha=0.8, ax=ax)
            raw_val = raw_row[col]
            ax.axhline(raw_val, color="black", linestyle="--",
                       linewidth=1.2, label=f"RAW ({raw_val:.3f})")
            ax.set_title(title, fontsize=10)
            ax.set_ylabel("")
            ax.set_xlabel("")
            ax.set_xticklabels(order, rotation=45, ha="right", fontsize=8)
            ax.legend(fontsize=8)
            ax.grid(axis="y", linestyle="--", alpha=0.4)
    
        plt.suptitle(
            f"Gene-Gene Co-expression Contrast [{level}]",
            fontsize=13
        )
        plt.tight_layout()
        out_fn = os.path.join(SAVE_DIR,
                              f"GeneGene_Contrast_Barplot_{level}_{file_suffix}.png")
        fig.savefig(out_fn, dpi=300)
        plt.close(fig)
        print(f"  [saved] {out_fn}")
    
        # ── Figure 2: heatmap subplot (tool별) ────────────────
        tools_ordered = ["RAW"] + [t for t in data_dict if t != "RAW"]
        ncols = 4
        nrows = int(np.ceil(len(tools_ordered) / ncols))
    
        fig, axes = plt.subplots(nrows, ncols,
                                  figsize=(6 * ncols, 6 * nrows))
        axes      = np.array(axes).flatten()
        palette   = sns.color_palette("husl", len(marker_dict))
        group_colors = {ct: c for ct, c in zip(marker_dict.keys(), palette)}
    
        for i, tool_name in enumerate(tools_ordered):
            ax    = axes[i]
            adata = data_dict[tool_name]
    
            avail_genes   = adata.var.index.tolist()
            valid_markers = [g for g in all_markers if g in avail_genes]
    
            X = adata[:, valid_markers].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            corr_mat = np.corrcoef(np.asarray(X, dtype=np.float64).T)
            corr_df  = pd.DataFrame(corr_mat,
                                     index=valid_markers,
                                     columns=valid_markers)
    
            row_data     = summary_df[summary_df["tool"] == tool_name]
            contrast_val = row_data["contrast"].values[0] if len(row_data) else np.nan
    
            sns.heatmap(
                corr_df, ax=ax,
                cmap="RdBu_r", center=0, vmin=-0.6, vmax=1.0,
                cbar=False, square=True,
                xticklabels=True, yticklabels=True,
                linewidths=0.3, linecolor="lightgray"
            )
    
            # cell type 구분선 + 색상 띠
            current = 0
            n_total = len(valid_markers)
            bar_w   = n_total * 0.04
    
            for ct, genes in marker_dict.items():
                n = sum(1 for g in genes if g in valid_markers)
                if n > 0:
                    color = group_colors[ct]
                    # top bar
                    ax.add_patch(mpatches.Rectangle(
                        (current, -bar_w), n, bar_w,
                        facecolor=color, edgecolor="none",
                        transform=ax.transData, clip_on=False
                    ))
                    # right bar
                    ax.add_patch(mpatches.Rectangle(
                        (n_total, current), bar_w, n,
                        facecolor=color, edgecolor="none",
                        transform=ax.transData, clip_on=False
                    ))
                    # 구분선
                    ax.axhline(current + n, color="black",
                               linewidth=1.0, linestyle="-")
                    ax.axvline(current + n, color="black",
                               linewidth=1.0, linestyle="-")
                    current += n
    
            ax.tick_params(axis="x", labelsize=7, rotation=90)
            ax.tick_params(axis="y", labelsize=7, rotation=0)
            ax.set_title(f"{tool_name}  |  contrast={contrast_val:.3f}",
                         fontsize=9, fontweight="bold", pad=20)
    
        for j in range(len(tools_ordered), len(axes)):
            axes[j].axis("off")
    
        legend_handles = [
            mpatches.Patch(color=c, label=ct)
            for ct, c in group_colors.items()
        ]
        fig.legend(handles=legend_handles, loc="center left",
                   bbox_to_anchor=(0.91, 0.5),
                   title="Cell Type", fontsize=10,
                   title_fontsize=12, frameon=False)
    
        plt.suptitle(
            f"Gene-Gene Correlation Heatmap [{level}]  |  "
            f"contrast = intra - inter  (higher is better)",
            fontsize=13, y=1.01
        )
        plt.tight_layout(rect=[0, 0, 0.90, 0.98])
        out_fn = os.path.join(SAVE_DIR,
                              f"GeneGene_Heatmap_{level}_{file_suffix}.png")
        fig.savefig(out_fn, dpi=200, bbox_inches="tight")
        plt.close(fig)
        print(f"  [saved] {out_fn}")
    
    return all_summary


# 실행
gene_gene_summary = run_gene_gene_contrast(file_suffix="zheng")

#7. Gene-gene Correlation: 산점도
 
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import os

# ===================================================================
# 1. 설정: 분석할 유전자 쌍
# ===================================================================
GENE_X = "CD8A"
GENE_Y = "CD8B"

# 전처리 기준
RAW_SCALE = ["DCA", "SAVER", "scIGANs", "scSTD", "scMASKGAN", "scMultiGAN", "RAW"]
SAVE_DIR = f"{BASE_DIR}/Zheng/figures"
os.makedirs(SAVE_DIR, exist_ok=True)

print(f"--- Scatter Plot Analysis: {GENE_X} vs {GENE_Y} ---")

# ===================================================================
# 2. 툴 순서 정의 (RAW + imputed_ads 순서)
# ===================================================================
# imputed_ads에 있는 키 순서대로 가져오되, 혹시 RAW가 들어있으면 제외
imputed_names = [k for k in imputed_ads.keys() if k != "RAW"]

# RAW를 맨 앞에 고정
tools_ordered = ['RAW'] + imputed_names

print(f"Plot Order: {tools_ordered}")

# ===================================================================
# 3. 시각화 루프
# ===================================================================
n_cols = 4
n_rows = int(np.ceil(len(tools_ordered) / n_cols))

fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
axes = axes.flatten()

for i, tool_name in enumerate(tools_ordered):
    ax = axes[i]
    
    # 1. 데이터 로드 (RAW 별도 처리)
    if tool_name == 'RAW':
        ad = raw.copy()
    else:
        ad = imputed_ads[tool_name].copy()
    
    # 유전자 존재 여부 확인
    if GENE_X not in ad.var_names or GENE_Y not in ad.var_names:
        ax.text(0.5, 0.5, "Gene not found", ha='center', va='center')
        continue
    
    # 2. 전처리 (Log1p)
    # 이미 Log되어 있는지 확인 (max < 20 이면 보통 Log된 상태)
    is_log_scale = (ad.X.max() < 20) or ("log1p" in ad.uns)
    
    if tool_name in RAW_SCALE and not is_log_scale:
        if tool_name == "afMF":
            # afMF 음수 클리핑
            X_mat = ad.X.toarray() if hasattr(ad.X, "toarray") else ad.X
            X_mat[X_mat < 0] = 0
            ad.X = X_mat
            
        sc.pp.normalize_total(ad, target_sum=1e4)
        sc.pp.log1p(ad)
        
    # 3. 발현량 추출
    # toarray()로 밀집 행렬 변환 후 1차원 벡터로 만듦
    x_val = ad[:, GENE_X].X.toarray().flatten() if hasattr(ad.X, 'toarray') else ad[:, GENE_X].X.flatten()
    y_val = ad[:, GENE_Y].X.toarray().flatten() if hasattr(ad.X, 'toarray') else ad[:, GENE_Y].X.flatten()
    
    # 4. 상관계수 계산
    if len(x_val) > 0 and len(y_val) > 0:
        corr, _ = pearsonr(x_val, y_val)
    else:
        corr = 0
    
    # 5. 산점도 그리기
    # s: 점 크기, alpha: 투명도 (데이터가 많으므로 투명하게)
    ax.scatter(x_val, y_val, s=3, alpha=0.3, c='black', edgecolors='none', rasterized=True)
    
    # 6. 스타일링
    ax.set_title(f"{tool_name}\nR = {corr:.3f}", fontsize=16, fontweight='bold')
    ax.set_xlabel(f"{GENE_X} Expression")
    ax.set_ylabel(f"{GENE_Y} Expression")
    
    # 0이 아닌 데이터의 비율(Sparsity Overlap) 표시
    non_zero_overlap = np.count_nonzero((x_val > 0) & (y_val > 0))
    non_zero_ratio = non_zero_overlap / len(x_val) * 100
    ax.text(0.05, 0.95, f"Co-expr: {non_zero_ratio:.1f}%", 
            transform=ax.transAxes, fontsize=11, color='red', va='top')
    
    # 가이드라인 (대각선)
    max_val = max(x_val.max(), y_val.max())
    if max_val > 0:
        ax.plot([0, max_val], [0, max_val], ls='--', c='gray', alpha=0.5)

# 빈 subplot 정리
for j in range(i + 1, len(axes)):
    axes[j].axis('off')

plt.suptitle(f"Gene-Gene Expression Scatter: {GENE_X} vs {GENE_Y}", fontsize=24, y=1.02, fontweight='bold')
plt.tight_layout()

# 저장
save_path = os.path.join(SAVE_DIR, f"Scatter_{GENE_X}_vs_{GENE_Y}.png")
plt.savefig(save_path, dpi=300, bbox_inches='tight')

print(f"[저장 완료] {save_path}")
