# Set the base directory for datasets
# Change this path to point to your local data directory
BASE_DIR = "../data"

import os
import scanpy as sc
import pandas as pd

os.chdir(f"{BASE_DIR}/Zheng")
zheng_raw = sc.read_10x_mtx(f"{BASE_DIR}/Zheng",var_names="gene_symbols", make_unique=True)
df_anno = pd.read_csv('./68k_pbmc_barcodes_annotation.tsv', sep='\t')
# 2. 'barcodes' 열을 인덱스로 설정 (중복이 있다면 제거해야 함)
# 이렇게 하면 나중에 zheng_raw의 바코드와 1:1로 매칭하기 쉬워집니다.
df_anno = df_anno.set_index('barcodes')
# 3. zheng_raw의 바코드(obs_names)를 기준으로 celltype을 찾아 매핑
# df_anno['celltype']은 바코드가 인덱스인 Series여야 합니다.
zheng_raw.obs['cell_type'] = zheng_raw.obs_names.map(df_anno['celltype'])

import scanpy as sc
import pandas as pd
import numpy as np

# 0. 데이터 복사 (원본 보존)
adata = zheng_raw.copy()

# ====================================================
# 1. 희귀 세포 집단 제거 (Subpopulation Filtering)
# 전체의 1% 미만인 세포 유형 제거 
# ====================================================
# 각 세포 유형의 비율 계산
type_counts = adata.obs['cell_type'].value_counts(normalize=True)

# 1% 이상인 세포 유형만 필터링
major_types = type_counts[type_counts >= 0.01].index
adata_filtered = adata[adata.obs['cell_type'].isin(major_types)].copy()

print(f"제거된 희귀 세포 유형: {list(type_counts[type_counts < 0.01].index)}")
print(f"남은 세포 수: {adata_filtered.n_obs}")

# ====================================================
# 2. 세포 품질 기반 순위 매기기 (Cell Ranking)
# 발현된 유전자 비율(fraction of expressed genes) 계산 
# ====================================================
# QC 메트릭 계산 (n_genes_by_counts 등을 얻음)
sc.pp.calculate_qc_metrics(adata_filtered, inplace=True)

# 전체 유전자 수 대비 발현된 유전자 비율 계산
adata_filtered.obs['frac_genes'] = adata_filtered.obs['n_genes_by_counts'] / adata_filtered.n_vars

# ====================================================
# 3. 비율 유지 샘플링 (Proportionate Sampling)
# 상위 품질 세포를 선택하여 총 10,000개 추출 
# ====================================================
target_total = 10000
final_indices = []

# 필터링된 데이터에서 비율 재계산
current_props = adata_filtered.obs['cell_type'].value_counts(normalize=True)

for cell_type in current_props.index:
    # 해당 세포 유형의 데이터만 추출
    subset = adata_filtered.obs[adata_filtered.obs['cell_type'] == cell_type]
    
    # 목표 할당량 계산 (비율 * 10,000)
    n_target = int(current_props[cell_type] * target_total)
    
    # 'frac_genes' 기준으로 내림차순 정렬 (품질 좋은 세포 우선)
    subset_sorted = subset.sort_values(by='frac_genes', ascending=False)
    
    # 상위 n_target 개만큼 인덱스 선택
    top_cells = subset_sorted.index[:n_target]
    final_indices.extend(top_cells)

# 최종 선택된 세포들로 AnnData 생성
adata_10k = adata_filtered[final_indices].copy()

print(f"최종 다운샘플링된 세포 수: {adata_10k.n_obs}")
print(adata_10k.obs['cell_type'].value_counts())
# 부족한 개수 계산 (여기서는 5개)
diff = target_total - len(final_indices)

if diff > 0:
    # 1. 가장 개수가 많은 세포 타입 찾기 (여기서는 CD8+ Cytotoxic T)
    largest_group = adata_10k.obs['cell_type'].value_counts().idxmax()
    
    # 2. 필터링된 원본 데이터(adata_filtered)에서 아직 선택되지 않은 해당 타입의 세포들 찾기
    already_selected = set(final_indices)
    remaining_in_largest = [
        idx for idx in adata_filtered.obs_names 
        if adata_filtered.obs.loc[idx, 'cell_type'] == largest_group 
        and idx not in already_selected
    ]
    
    # 3. 부족한 만큼 더 선택 (품질 순서대로 정렬된 상태를 원한다면 정렬 후 추출)
    extra_indices = remaining_in_largest[:diff]
    final_indices.extend(extra_indices)

# 다시 AnnData 생성
adata_10k = adata_filtered[final_indices].copy()

print(f"보정 후 최종 세포 수: {adata_10k.n_obs}")
print(adata_10k.obs['cell_type'].value_counts())
# ====================================================
# 4. 고변동성 유전자 선별 (HVG Selection)
# 4가지 버전 생성 (예: 2000개 버전) [cite: 2899]
# ====================================================
# 정규화 및 로그 변환 (HVG 계산 전 필수)
#sc.pp.normalize_total(adata_10k, target_sum=1e4)
#sc.pp.log1p(adata_10k)

# HVG 선별 (예: 상위 2000개)
#sc.pp.highly_variable_genes(adata_10k, n_top_genes=2000, subset=True)

#print(f"최종 데이터 형태 (Cells x Genes): {adata_10k.shape}")

import scanpy as sc
import os

# 0. 데이터 복사 (원본 보존)
adata = adata_10k.copy()  # 또는 zheng_raw.copy()

# ====================================================
# 1. 기본 필터링 (품질 낮은 세포/유전자 제거)
# 참고: Huang et al. 2025 기준 (min_genes=200, min_cells=10)
# 1. 저장할 경로 설정
save_dir = f'{BASE_DIR}/Zheng'
# 2. Scanpy의 그림 저장 경로 지정
sc.settings.figdir = save_dir

# ====================================================

# 세포 필터링: 유전자가 200개 미만으로 발현된 죽은 세포 제거
sc.pp.filter_cells(adata, min_genes=200)

# 유전자 필터링: 3개 미만의 세포에서만 발견된 희귀 유전자 제거
sc.pp.filter_genes(adata, min_cells=3)

print(f"기본 필터링 후 크기: {adata.shape}")


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


# ====================================================
# 3. 추가 아웃라이어 제거 (논문 기준 적용)
# ====================================================

# 기준 1: 미토콘드리아 비율이 너무 높은 세포 제거 (보통 5~20% 기준)
# PBMC는 깨끗한 편이라 5%~10%를 많이 씁니다. (Huang et al.은 outlier filtering 언급)
adata = adata[adata.obs['pct_counts_mt'] < 10, :]

# 기준 2: 유전자 수가 비정상적으로 많은 세포(Multiplet 가능성) 제거
# Cheng et al.은 75th percentile + 3*IQR 방식을 썼으나, 보통 2500~5000개로 자르기도 합니다.
# 여기서는 상위 1% 정도를 자르는 예시입니다.
#max_genes = adata.obs['n_genes_by_counts'].quantile(0.99)
#adata = adata[adata.obs['n_genes_by_counts'] < max_genes, :]
#print(f"QC 후 최종 세포 수: {adata.n_obs}")


# ====================================================
# 4. 정규화 및 로그 변환 (필수)
# Huang et al. 등 모든 논문의 Baseline
# ====================================================
# 라이브러리 사이즈 맞추기 (Target sum = 10,000)
#sc.pp.normalize_total(adata, target_sum=1e4)
# 로그 변환 (log(x+1))
#sc.pp.log1p(adata)
# ====================================================
# 5. (옵션) HVG 선별
# Cheng et al.의 대용량 데이터 처리 방식 또는 scSTD 방식
# ====================================================
# 상위 2,000개 고변동성 유전자만 남기기
#sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
#print(f"최종 분석용 데이터 형태: {adata.shape}")

# Save the processed data
#adata.write(os.path.join(save_dir, 'zheng.h5ad'))

save_dir = f'{BASE_DIR}/Zheng/Raw'
zheng = sc.read_h5ad(os.path.join(save_dir,"zheng.h5ad"))


# 2. 매핑 딕셔너리 정의 (CD4 vs CD8 구분)
mapping_dict = {
    # --- CD8+ T cells 계열 ---
    'CD8+ Cytotoxic T': 'CD8+ T cells',
    'CD8+/CD45RA+ Naive Cytotoxic': 'CD8+ T cells',
    
    # --- CD4+ T cells 계열 ---
    'CD4+/CD25 T Reg': 'CD4+ T cells',
    'CD4+/CD45RO+ Memory': 'CD4+ T cells',
    'CD4+/CD45RA+/CD25- Naive T': 'CD4+ T cells',
    
    # --- 기타 주요 세포 ---
    'CD56+ NK': 'NK cells',
    'CD19+ B': 'B cells',
    'CD14+ Monocyte': 'Monocytes',
    'Dendritic': 'Dendritic cells'
}

# 3. 새로운 컬럼 'major_cell_type' 생성
# (혹시 모를 에러 방지를 위해 string으로 변환 후 매핑)
zheng.obs['major_cell_type'] = zheng.obs['cell_type'].astype(str).map(mapping_dict)

# 4. 결과 확인
print("=== 분류 전 (Original) ===")
print(zheng.obs['cell_type'].value_counts())
print("\n=== 분류 후 (Major Cell Type) ===")
print(zheng.obs['major_cell_type'].value_counts())
