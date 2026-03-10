# Set the base directory for datasets
# Change this path to point to your local data directory
BASE_DIR = "../data"

import scanpy as sc
import pandas as pd
import os
import glob
import scipy.sparse as sp

# ================== 설정 ==================
INPUT_DIR = f"{BASE_DIR}/Chu_time/Raw"
OUTPUT_DIR = f"{BASE_DIR}/Chu_time/Raw/scIDPMs_input"
os.makedirs(OUTPUT_DIR, exist_ok=True)
# ==========================================

file_list = sorted(glob.glob(os.path.join(INPUT_DIR, "*.h5ad")))
print(f"📂 총 {len(file_list)}개 h5ad 파일 발견\n")
for file_path in file_list:
    file_name = os.path.basename(file_path)
    name_base = os.path.splitext(file_name)[0]
    print(f"🔄 처리 중: {file_name}")
    try:
        adata = sc.read_h5ad(file_path)
        # ====== [중요] feature / cell 축 고정 ======
        adata = adata[:, adata.var_names].copy()
        adata.obs_names = adata.obs_names.astype(str)
        adata.var_names = adata.var_names.astype(str)
        # ====== count matrix ======
        X = adata.X
        if sp.issparse(X):
            X = X.toarray()
        df_counts = pd.DataFrame(
            X,
            index=adata.obs_names,      # cell barcode
            columns=adata.var_names     # ENSG gene id
        )
        count_path = os.path.join(OUTPUT_DIR, f"{name_base}_counts.csv")
        # 🔴 핵심: header=True, index=True (명시)
        df_counts.to_csv(
            count_path,
            index=True,
            header=True
        )
        print(f"   ✅ Counts 저장 완료: {df_counts.shape}")
        # ====== label ======
        if "cell_type" in adata.obs.columns:
            df_labels = adata.obs.loc[:, ["cell_type"]]
            label_path = os.path.join(OUTPUT_DIR, f"{name_base}_label.csv")
            df_labels.to_csv(
                label_path,
                index=True,
                header=True
            )
            print("   ✅ Label 저장 완료")
        else:
            print("   ⚠️ cell_type 없음 → label 생략")
    except Exception as e:
        print(f"   ❌ 오류 발생: {e}")
    print("-" * 50)

print("✨ 모든 작업 완료")



# Normalized Data 생성 코드
import scanpy as sc
import pandas as pd
import os
import glob
import scipy.sparse as sp

# ================== 설정 ==================
INPUT_DIR = f"{BASE_DIR}/Chu/Raw"

# ⚠️ 저장 경로를 변경했습니다 (_norm 추가) -> 섞이지 않게 하기 위함
OUTPUT_DIR = f"{BASE_DIR}/Chu/Raw/scIDPMs_input_norm"
os.makedirs(OUTPUT_DIR, exist_ok=True)
# ==========================================

file_list = sorted(glob.glob(os.path.join(INPUT_DIR, "*.h5ad")))
print(f"📂 총 {len(file_list)}개 h5ad 파일 발견\n")

for file_path in file_list:
    file_name = os.path.basename(file_path)
    name_base = os.path.splitext(file_name)[0]
    print(f"🔄 처리 중: {file_name}")
    
    try:
        adata = sc.read_h5ad(file_path)
        
        # ====== [중요] feature / cell 축 고정 ======
        adata = adata[:, adata.var_names].copy()
        adata.obs_names = adata.obs_names.astype(str)
        adata.var_names = adata.var_names.astype(str)
        
        # ====== [★ 핵심 추가] Normalization Process ======
        print("   ⚡ Normalizing (10k) & Log1p Transformation 적용 중...")
        
        # 1. Library Size Normalization
        # 세포마다 읽힌 총 read 수를 10,000(1e4)으로 맞춰줍니다. (Depth 차이 보정)
        sc.pp.normalize_total(adata)
        # 2. Log Transformation
        # 데이터 분포를 펴줍니다 (log(x+1)). 폭죽(Starburst) 현상 방지 핵심
        sc.pp.log1p(adata)
        # ===============================================
    
        # ====== count matrix ======
        X = adata.X
        
        # Sparse Matrix인 경우 Dense로 변환 (CSV 저장을 위해)
        if sp.issparse(X):
            X = X.toarray()
            
        df_counts = pd.DataFrame(
            X,
            index=adata.obs_names,      # cell barcode
            columns=adata.var_names     # ENSG gene id
        )
        
        count_path = os.path.join(OUTPUT_DIR, f"{name_base}_counts.csv")
        
        # 🔴 핵심: header=True, index=True (명시)
        df_counts.to_csv(
            count_path,
            index=True,
            header=True
        )
        print(f"   ✅ Processed Counts 저장 완료: {df_counts.shape}")
        
        # ====== label ======
        if "cell_type" in adata.obs.columns:
            df_labels = adata.obs.loc[:, ["cell_type"]]
            label_path = os.path.join(OUTPUT_DIR, f"{name_base}_label.csv")
            df_labels.to_csv(
                label_path,
                index=True,
                header=True
            )
            print("   ✅ Label 저장 완료")
        else:
            print("   ⚠️ cell_type 없음 → label 생략")
            
    except Exception as e:
        print(f"   ❌ 오류 발생: {e}")
        
    print("-" * 50)

print("✨ 모든 작업 완료 (Normalized Data 생성됨)")

# 연구실 서버에서 KBDS로 옮기기
rsync -avzP --progress -e "ssh -p 22" \
  $DATA_DIR/Chu_time/Raw/scIDPMs_input \
  user@server:/home01/user/Datasets/Chu_time/Raw
rsync -avzP --progress -e "ssh -p 22" \
  $BENCHMARK_DIR/scIDPMs/src/utils_table.py \
  user@server:/home01/user/scIDPMs/src/
rsync -avzP --progress -e "ssh -p 22" \
  $BENCHMARK_DIR/scIDPMs/exe.py \
  user@server:/home01/user/scIDPMs/
rsync -avzP --progress -e "ssh -p 22" \
  $BENCHMARK_DIR/scIDPMs/dataset.py \
  user@server:/home01/user/scIDPMs/


# KBDS에서 서버로 옮기기
scp user@server:/home01/user/scIDPMs/save/scidpms_zheng_20260220_194419/imputed.csv $DATA_DIR/Zheng/imputed/scIDPMs/
scp user@server:/home01/user/scIDPMs/save/scidpms_zheng_norm_20260129_205129/imputed.csv $DATA_DIR/Zheng/imputed/scIDPMs_norm/
scp user@server:/home01/user/scIDPMs/save/scidpms_chu_20260131_015958/imputed.csv $DATA_DIR/Chu/imputed/scIDPMs/
scp user@server:/home01/user/scIDPMs/save/scidpms_chu_time_20260220_194413/imputed.csv $DATA_DIR/Chu_time/imputed/scIDPMs/

# KBDS에서 서버로 옮기기
scp user@server:/home01/user/scIDPMs/save/sim1_counts_LSH_20251119_195420/imputed.csv $DATA_DIR/splatter/imputed/scIDPMs/sim1/
scp user@server:/home01/user/5K_PBMC_10X/Raw/scIDPMs/save/sim2_counts_LSH_20251119_201145/imputed.csv $DATA_DIR/5K_PBMC_10X/imputed/scIDPMs/sim2/
scp user@server:/home01/user/5K_PBMC_10X/Raw/scIDPMs/save/sim3_counts_LSH_20251119_210607/imputed.csv $DATA_DIR/5K_PBMC_10X/imputed/scIDPMs/sim3/


scp user@server:/home01/user/scIDPMs/save/sim1_counts_LSH_20251201_141042/imputed.csv $DATA_DIR/splatter/imputed/scIDPMs/sim1/

# 로컬에서 KBDS 서버로 옮기기
scp -r ~/다운로드/scidpm_edited/src/utils_table.py  user@server:/home01/user/scIDPMs/src/
rsync -avzP --progress -e "ssh -p 22" \
  ~/다운로드/scidpm_edited/scIDPMs \
  user@server:/home01/user/




#!/bin/bash -l
#SBATCH --job-name=scidpms_chu_time
#SBATCH --partition=4gpu
#SBATCH --gres=gpu:1
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --output=/home01/user/scIDPMs/logs/%x-%j.out
cd /home01/user/scIDPMs
module load compilers/cuda/12.4
module load libraries/nccl/2.21.5
module load applications/python/3.10.4

# conda 환경 활성화
conda activate scIDPMs

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True,max_split_size_mb:128
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

unset SCIDPMS_DEBUG_GEN

# 1. 시작 시점의 타이머를 0으로 초기화 (또는 현재 시간 기록)
SECONDS=0
echo "Imputation Start: $(date)"

python /home01/user/scIDPMs/exe.py \
  --config default.yaml \
  --device cuda \
  --att LSH \
  --nsample 10 \
  --file_path /home01/user/Datasets/Chu_time/Raw/scIDPMs_input/chu_counts.csv \
  --label_path /home01/user/Datasets/Chu_time/Raw/scIDPMs_input/chu_label.csv \

# 2. 걸린 시간 계산 및 포맷팅 출력
duration=$SECONDS
hours=$(($duration / 3600))
minutes=$(($duration % 3600 / 60))
seconds=$(($duration % 60))

echo "------------------------------------------------"
echo "Job Finished."
echo "Total Execution Time: ${hours}h ${minutes}m ${seconds}s" 
echo "------------------------------------------------"
echo "Imputation End: $(date)"

#실행 명령어
python exe.py \
    --config lsh_safe.yaml \
    --device cuda:5 \
    --att LSH \
    --nsample 10 \
    --file_path='$DATA_DIR/5K_PBMC_10X/Raw/scIDPMs_input/pbmc_counts.csv' \
    --label_path='$DATA_DIR/5K_PBMC_10X/Raw/scIDPMs_input/pbmc_label.csv'

import scanpy as sc
import numpy as np
import pandas as pd
import os

# 데이터 로드
in_dir = f"{BASE_DIR}/Chu_time/Raw"
chu = sc.read_h5ad(in_dir + "/chu.h5ad")

# 요청하신 저장 경로 설정
output_dir = f"{BASE_DIR}/Chu_time/Raw/scSTD_input"
# 훈련 및 Imputation에 사용될 데이터셋 목록
adata_list = [
    (chu, "chu.txt")
]

print("✅ AnnData 객체에서 TXT 파일로 데이터 추출 및 저장 시작...")

for i, (chu, filename) in enumerate(adata_list):
    # 1. 데이터 추출: .X 필드 사용 (셀 x 유전자 행렬)
    data = chu.X
    
    # 2. 희소 행렬인 경우 조밀 행렬로 변환
    if not isinstance(data, np.ndarray):
        data = data.toarray()
    
    # 3. 데이터 저장: autoencoder.py에서 기대하는 쉼표(,) 구분자로 저장
    output_path = os.path.join(output_dir, filename)
    pd.DataFrame(data).to_csv(output_path, sep=',', header=False, index=False)
    print(f"   -> [{i+1}/1] {filename} ({data.shape[0]} cells, {data.shape[1]} genes) 저장 완료: {output_path}")

print("\n✅ 모든 TXT 파일 저장 완료.")

cd $BENCHMARK_DIR/scstd
conda activate scstd

# chu_time 데이터로 훈련 시작
time python autoencoder.py \
    --file_path $DATA_DIR/Chu_time/Raw/scSTD_input/chu.txt \
    --gpu_id 1 \
    --backup_dir $DATA_DIR/Chu_time/imputed/scSTD

# chu 데이터로 훈련 시작
time python autoencoder.py \
    --file_path $DATA_DIR/Chu/Raw/scSTD_input/chu.txt \
    --gpu_id 1 \
    --backup_dir $DATA_DIR/Chu/imputed/scSTD


#zheng 데이터로 훈련 시작
time python autoencoder.py \
    --file_path $DATA_DIR/Zheng/Raw/scSTD_input/zheng.txt \
    --gpu_id 3 \
    --backup_dir $DATA_DIR/Zheng/imputed/scSTD
#pbmc 데이터로 훈련 시작
python autoencoder.py \
    --file_path $DATA_DIR/5K_PBMC_10X/Raw/scSTD_input/pbmc.txt \
    --gpu_id 7 \
    --backup_dir $DATA_DIR/5K_PBMC_10X/imputed/scSTD
# Sim1 데이터로 훈련 시작
python autoencoder.py \
    --file_path $DATA_DIR/splatter/Raw/scSTD_input/sim1.txt \
    --gpu_id 5 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD
# Sim2 데이터로 훈련 시작
python autoencoder.py \
    --file_path $DATA_DIR/splatter/Raw/scSTD_input/sim2.txt \
    --gpu_id 6 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD
# Sim3 데이터로 훈련 시작
python autoencoder.py \
    --file_path $DATA_DIR/splatter/Raw/scSTD_input/sim3.txt \
    --gpu_id 7 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD

import numpy as np
import scanpy as sc
import anndata as ad
import os
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics import adjusted_rand_score

# Matplotlib 백엔드 설정 (서버 환경을 위해 유지)
matplotlib.use('Agg')

# --- 1. 경로 및 상수 설정 ---
BASE_DIR = Path(f"{BASE_DIR}/splatter")
EMBEDDINGS_SUBDIR = BASE_DIR / "imputed" / "scSTD"
ORIGINAL_SUBDIR = BASE_DIR / "Raw"
DATASETS = ["sim1", "sim2", "sim3"]
# 원래 기대하는 obs 컬럼 후보
COLOR_BY_CAND = ['Group', 'group', 'cell_type', 'sub_cell_type']
EMBEDDINGS_SUBDIR.mkdir(parents=True, exist_ok=True)
(EMBEDDINGS_SUBDIR / "Debug").mkdir(parents=True, exist_ok=True) # Debug 디렉토리 생성
print(f"현재 작업 디렉토리: {os.getcwd()}")
print(f"원본 데이터 디렉토리: {ORIGINAL_SUBDIR}")
print(f"임베딩 디렉토리: {EMBEDDINGS_SUBDIR}")
# --- 2. 원본 AnnData 객체 로드 (메모리에 유지) ---
print("✅ 원본 AnnData 객체 로드 중...")
try:
    sim1 = sc.read_h5ad(ORIGINAL_SUBDIR / "sim1.h5ad")
    sim2 = sc.read_h5ad(ORIGINAL_SUBDIR / "sim2.h5ad")
    sim3 = sc.read_h5ad(ORIGINAL_SUBDIR / "sim3.h5ad")
    adata_objects = {'sim1': sim1, 'sim2': sim2, 'sim3': sim3}
    print("✅ AnnData 객체 로드 완료.")
except FileNotFoundError as e:
    print(f"❌ AnnData 파일 로드 오류: {e}. 경로를 확인하세요.")
    exit()
# --- 3. 수정된 시각화 및 평가 함수 ---
def visualize_and_save_embeddings_single(dataset_name: str, original_adata: ad.AnnData):
    print(f"\n--- 데이터셋: {dataset_name} 처리 시작 ---")
    embeddings_path = EMBEDDINGS_SUBDIR / f"{dataset_name}_embeddings.npy"
    # 3.1. 임베딩 로드
    if not embeddings_path.exists():
        print(f"   ❌ 오류: 임베딩 파일이 {embeddings_path} 에 없습니다. 건너뜀.")
        return
    embeddings = np.load(embeddings_path)
    print(f"   -> 임베딩 로드 완료. 형태: {embeddings.shape}")
    # 임베딩 데이터가 원본 세포 수와 일치하는지 확인
    if embeddings.shape[0] != original_adata.n_obs:
        print(f"   ❌ 오류: 임베딩 수 ({embeddings.shape[0]})가 원본 세포 수 ({original_adata.n_obs})와 불일치합니다. 건너뜀.")
        return
    # 3.2. obs 컬럼 선택 및 AnnData 생성
    available_cols = list(original_adata.obs.columns)
    # ARI 계산에 사용할 수 있는 Group Label 찾기
    group_key = next((k for k in COLOR_BY_CAND if k in available_cols), None) 
    # obs 데이터 복사 (원본 라벨 유지)
    obs_data = original_adata.obs.copy()
    embeddings_adata = ad.AnnData(X=embeddings, obs=obs_data)
    print("   -> 임베딩 AnnData 생성 완료.")
    # 3.3. neighbors/UMAP/Leiden 계산 (t-SNE 제거)
    print("   -> neighbors / UMAP / Leiden 계산 시작...")
    sc.pp.neighbors(embeddings_adata, use_rep='X', n_neighbors=15)
    sc.tl.umap(embeddings_adata)
    # 클러스터링 수행 (ARI 계산을 위해 필요)
    sc.tl.leiden(embeddings_adata, key_added="leiden") 
    print("   -> 계산 완료.")
    # 3.4. ARI 계산
    ari_score = None
    color_key_to_use = None
    if group_key:
        # Group Label과 Leiden 클러스터링 결과 비교
        y_true = embeddings_adata.obs[group_key].astype("category")
        y_pred = embeddings_adata.obs["leiden"].astype("category")
        ari_score = adjusted_rand_score(y_true, y_pred)
        color_key_to_use = group_key
        print(f"   ✅ ARI ({group_key} vs Leiden): {ari_score:.4f}")
    else:
        print(f"   ⚠ {COLOR_BY_CAND} 중 유효한 그룹 키를 찾지 못했습니다. ARI 계산 생략.")
        # 만약 그룹 키가 없다면, Leiden 결과를 시각화 키로 사용
        color_key_to_use = "leiden"
    # 3.5. UMAP 시각화 & 저장 (t-SNE 제거)
    output_filename = f"{dataset_name}_embeddings_umap.png"
    save_path = EMBEDDINGS_SUBDIR / "Debug" / output_filename
    # UMAP 플롯 생성
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    try:
        title_str = f"UMAP of {dataset_name} Autoencoder Embeddings"
        if ari_score is not None:
            title_str += f" (ARI={ari_score:.4f})"
        sc.pl.embedding(
            embeddings_adata,
            basis='umap',
            color=color_key_to_use, # Group이 있으면 Group으로, 없으면 Leiden으로 색상 지정
            title=title_str,
            ax=ax,
            show=False,
            save=False,
            legend_loc='right margin' # 범례를 플롯 바깥쪽에 배치
        )
        plt.tight_layout()
        fig.savefig(save_path, dpi=300)
        plt.close(fig)
        print(f"   ✅ UMAP 이미지 저장 완료: {save_path}")
    except Exception as e:
        plt.close(fig)
        print(f"   ❌ 오류: UMAP 플롯 생성 중 오류 발생: {e}")
# --- 4. sim1, sim2, sim3 시각화 실행 ---
if __name__ == "__main__":
    for dataset in DATASETS:
        print("\n" + "=" * 50)
        print(f"--- 🔬 {dataset} 임베딩 시각화 및 평가 시작 ---")
        visualize_and_save_embeddings_single(dataset, adata_objects[dataset])
        print("=" * 50)

cd $BENCHMARK_DIR/scstd


# chu_time 데이터로 훈련 시작
time CUDA_VISIBLE_DEVICES=0 python sc_train.py \
  --embeddings_file $DATA_DIR/Chu_time/imputed/scSTD/chu_embeddings.npy \
  --gpu_id 0 \
  --batch_size 64 \
  --microbatch 16 \
  --use_fp16 False \
  --lr 0.001 \
  --lr_anneal_steps 100000 \
  --save_interval 10000 \
  --backup_dir $DATA_DIR/Chu_time/imputed/scSTD



# chu 데이터로 훈련 시작
time CUDA_VISIBLE_DEVICES=0 python sc_train.py \
  --embeddings_file $DATA_DIR/Chu/imputed/scSTD/chu_embeddings.npy \
  --gpu_id 0 \
  --batch_size 64 \
  --microbatch 16 \
  --use_fp16 False \
  --lr 0.001 \
  --lr_anneal_steps 100000 \
  --save_interval 10000 \
  --backup_dir $DATA_DIR/Chu/imputed/scSTD



# zheng 데이터로 훈련 시작
time CUDA_VISIBLE_DEVICES=3 python sc_train.py \
  --embeddings_file $DATA_DIR/Zheng/imputed/scSTD/zheng_embeddings.npy \
  --gpu_id 3 \
  --batch_size 64 \
  --microbatch 16 \
  --use_fp16 False \
  --lr 0.001 \
  --lr_anneal_steps 100000 \
  --save_interval 10000 \
  --backup_dir $DATA_DIR/Zheng/imputed/scSTD

# Sim1 데이터로 훈련 시작
CUDA_VISIBLE_DEVICES=5 python sc_train.py \
  --embeddings_file $DATA_DIR/splatter/imputed/scSTD/sim1_embeddings.npy \
  --gpu_id 5 \
  --batch_size 64 \
  --microbatch 16 \
  --use_fp16 False \
  --lr 0.001 \
  --lr_anneal_steps 100000 \
  --save_interval 10000 \
  --backup_dir $DATA_DIR/splatter/imputed/scSTD

# Sim2 데이터로 훈련 시작
CUDA_VISIBLE_DEVICES=6 python sc_train.py \
  --embeddings_file $DATA_DIR/splatter/imputed/scSTD/sim2_embeddings.npy \
  --gpu_id 6 \
  --batch_size 64 \
  --microbatch 16 \
  --use_fp16 False \
  --lr 0.001 \
  --lr_anneal_steps 100000 \
  --save_interval 10000 \
  --backup_dir $DATA_DIR/splatter/imputed/scSTD

# Sim3 데이터로 훈련 시작
CUDA_VISIBLE_DEVICES=7 python sc_train.py \
  --embeddings_file $DATA_DIR/splatter/imputed/scSTD/sim3_embeddings.npy \
  --gpu_id 7 \
  --batch_size 64 \
  --microbatch 16 \
  --use_fp16 False \
  --lr 0.001 \
  --lr_anneal_steps 100000 \
  --save_interval 10000 \
  --backup_dir $DATA_DIR/splatter/imputed/scSTD

# chu_time 데이터로 샘플링
time CUDA_VISIBLE_DEVICES=0 python sc_sample.py \
    --model_path $DATA_DIR/Chu_time/imputed/scSTD/Diffusion_models/chu/ema_0.9999_100000.pt \
    --num_samples 10000 \
    --gpu_id 0 \
    --batch_size 100 \
    --backup_dir $DATA_DIR/Chu_time/imputed/scSTD


# chu 데이터로 샘플링
time CUDA_VISIBLE_DEVICES=7 python sc_sample.py \
    --model_path $DATA_DIR/Chu/imputed/scSTD/Diffusion_models/chu/ema_0.9999_100000.pt \
    --num_samples 10000 \
    --gpu_id 7 \
    --batch_size 100 \
    --backup_dir $DATA_DIR/Chu/imputed/scSTD

# zheng 데이터로 샘플링
time CUDA_VISIBLE_DEVICES=7 python sc_sample.py \
    --model_path $DATA_DIR/Zheng/imputed/scSTD/Diffusion_models/zheng/ema_0.9999_100000.pt \
    --num_samples 10000 \
    --gpu_id 7 \
    --batch_size 100 \
    --backup_dir $DATA_DIR/Zheng/imputed/scSTD

# Sim1 데이터로 샘플링
CUDA_VISIBLE_DEVICES=5 python sc_sample.py \
    --model_path $DATA_DIR/splatter/imputed/scSTD/Diffusion_models/sim1/ema_0.9999_100000.pt \
    --num_samples 10000 \
    --gpu_id 5 \
    --batch_size 100 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD
# Sim2 데이터로 샘플링
CUDA_VISIBLE_DEVICES=6 python sc_sample.py \
    --model_path $DATA_DIR/splatter/imputed/scSTD/Diffusion_models/sim2/ema_0.9999_100000.pt \
    --num_samples 10000 \
    --gpu_id 6 \
    --batch_size 100 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD
# Sim3 데이터로 샘플링
CUDA_VISIBLE_DEVICES=7 python sc_sample.py \
    --model_path $DATA_DIR/splatter/imputed/scSTD/Diffusion_models/sim3/ema_0.9999_100000.pt \
    --num_samples 10000 \
    --gpu_id 7 \
    --batch_size 100 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD

# Sim1 임퓨테이션
CUDA_VISIBLE_DEVICES=7 python imputation.py \
    --original_txt_path $DATA_DIR/splatter/Raw/scSTD_input/sim1.txt \
    --decoder_path $DATA_DIR/splatter/imputed/scSTD/AE_models/sim1_decoder.pth \
    --sampled_codes_path $DATA_DIR/splatter/imputed/scSTD/Sample/sim1/sim1_sampled_10000x256.npy \
    --query_embeddings_path $DATA_DIR/splatter/imputed/scSTD/sim1_embeddings.npy \
    --gpu_id 7 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD
# Sim2 임퓨테이션
CUDA_VISIBLE_DEVICES=6 python imputation.py \
    --original_txt_path $DATA_DIR/splatter/Raw/scSTD_input/sim2.txt \
    --decoder_path $DATA_DIR/splatter/imputed/scSTD/AE_models/sim2_decoder.pth \
    --sampled_codes_path $DATA_DIR/splatter/imputed/scSTD/Sample/sim2/sim2_sampled_10000x256.npy \
    --query_embeddings_path $DATA_DIR/splatter/imputed/scSTD/sim2_embeddings.npy \
    --gpu_id 6 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD
# Sim3 임퓨테이션
CUDA_VISIBLE_DEVICES=7 python imputation.py \
    --original_txt_path $DATA_DIR/splatter/Raw/scSTD_input/sim3.txt \
    --decoder_path $DATA_DIR/splatter/imputed/scSTD/AE_models/sim3_decoder.pth \
    --sampled_codes_path $DATA_DIR/splatter/imputed/scSTD/Sample/sim3/sim3_sampled_10000x256.npy \
    --query_embeddings_path $DATA_DIR/splatter/imputed/scSTD/sim3_embeddings.npy \
    --gpu_id 7 \
    --backup_dir $DATA_DIR/splatter/imputed/scSTD

# Zheng 임퓨테이션
CUDA_VISIBLE_DEVICES=7 python imputation.py \
    --original_txt_path $DATA_DIR/Zheng/Raw/scSTD_input/zheng.txt \
    --decoder_path $DATA_DIR/Zheng/imputed/scSTD/AE_models/zheng_decoder.pth \
    --sampled_codes_path $DATA_DIR/Zheng/imputed/scSTD/Sample/zheng/zheng_sampled_10000x256.npy \
    --query_embeddings_path $DATA_DIR/Zheng/imputed/scSTD/zheng_embeddings.npy \
    --gpu_id 7 \
    --backup_dir $DATA_DIR/Zheng/imputed/scSTD

# Chu_time 임퓨테이션
CUDA_VISIBLE_DEVICES=0 python imputation.py \
    --original_txt_path $DATA_DIR/Chu_time/Raw/scSTD_input/chu.txt \
    --decoder_path $DATA_DIR/Chu_time/imputed/scSTD/AE_models/chu_decoder.pth \
    --sampled_codes_path $DATA_DIR/Chu_time/imputed/scSTD/Sample/chu/chu_sampled_10000x256.npy \
    --query_embeddings_path $DATA_DIR/Chu_time/imputed/scSTD/chu_embeddings.npy \
    --gpu_id 0 \
    --backup_dir $DATA_DIR/Chu_time/imputed/scSTD

# scSTD 폴더 KBDS 서버로 옮기기
```rsync -avzP --progress -e "ssh -p 22" \
  $BENCHMARK_DIR/scstd \
  user@server:/home01/user/```

# 폴더 KBDS 서버로 옮기기: input
rsync -avzP --progress -e "ssh -p 22" \
  $DATA_DIR/Chu_time/Raw/scSTD_input \
  user@server:/home01/user/Datasets/Chu_time/Raw
# 폴더 KBDS 서버로 옮기기: imputed
rsync -avzP --progress -e "ssh -p 22" \
  $DATA_DIR/Chu_time/imputed/scSTD \
  user@server:/home01/user/Datasets/Chu_time/imputed  
# ========================================== KBDS 접속 ==========================================
#!/bin/bash -l
#SBATCH --job-name=scstd_impute_chu_time
#SBATCH --partition=debug-1gpu
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/home01/user/scstd/logs/%x-%j.out
cd /home01/user/scstd
module load compilers/cuda/12.4

# conda 환경 활성화
conda activate scSTD

# 1. 시작 시점의 타이머를 0으로 초기화 (또는 현재 시간 기록)
SECONDS=0
echo "Imputation Start: $(date)"

#pbmc 임퓨테이션
python imputation.py \
    --original_txt_path /home01/user/Datasets/Chu_time/Raw/scSTD_input/chu.txt \
    --decoder_path /home01/user/Datasets/Chu_time/imputed/scSTD/AE_models/chu_decoder.pth \
    --sampled_codes_path /home01/user/Datasets/Chu_time/imputed/scSTD/Sample/chu/chu_sampled_10000x256.npy \
    --query_embeddings_path /home01/user/Datasets/Chu_time/imputed/scSTD/chu_embeddings.npy \
    --backup_dir /home01/user/Datasets/Chu_time/imputed/scSTD
# 2. 걸린 시간 계산 및 포맷팅 출력
duration=$SECONDS
hours=$(($duration / 3600))
minutes=$(($duration % 3600 / 60))
seconds=$(($duration % 60))

echo "------------------------------------------------"
echo "Job Finished."
echo "Total Execution Time: ${hours}h ${minutes}m ${seconds}s" 
echo "------------------------------------------------"


# KBDS에서 서버로 옮기기
scp user@server:/home01/user/Datasets/Chu/imputed/scSTD/Result/chu_imputed.txt $DATA_DIR/Chu/imputed/scSTD/Result/

#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
scSTD 파이프라인 디버깅용 (torch 없이)

체크하는 것:
1) latent(real AE embedding) vs latent(sampled) UMAP
2) Raw vs scSTD-imputed UMAP + ARI(Group vs Leiden)

쓰기 전에 경로만 한 번 확인해봐라.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.metrics import adjusted_rand_score
import matplotlib.pyplot as plt

# --- 0. 경로 설정 ---

BASE_DIR = Path(f"{BASE_DIR}/splatter")
RAW_H5AD = BASE_DIR / "Raw" / "sim1.h5ad"  # AnnData (Group 정보 포함)
SCSTD_DIR = BASE_DIR / "imputed" / "scSTD"

SCSTD_INPUT_TXT = BASE_DIR / "Raw" / "scSTD_input" / "sim1.txt"  # scSTD용 입력
AE_EMBED_PATH = SCSTD_DIR / "sim1_embeddings.npy"                # AE latent
SAMPLED_LATENT_PATH = SCSTD_DIR / "Sample" / "sim1" / "sim1_sampled_10000x256.npy"

RESULT_DIR = SCSTD_DIR / "Result"
SCSTD_IMPUTED_PATH = RESULT_DIR / "sim1_imputed.txt"             # imputation.py 결과

DEBUG_DIR = SCSTD_DIR / "Debug"
DEBUG_DIR.mkdir(parents=True, exist_ok=True)


# --- 1. 공통: matrix → AnnData → UMAP + Leiden + ARI ---

def build_adata_and_umap(matrix: np.ndarray,
                         raw_adata: ad.AnnData,
                         label_key: str = "Group",
                         title_prefix: str = "",
                         out_png: Path = None):
    """
    matrix: (n_cells, n_genes) count-like matrix
    raw_adata.obs에서 label_key 사용 (예: 'Group')
    """
    assert matrix.shape[0] == raw_adata.n_obs, \
        f"Cell 수 mismatch: matrix {matrix.shape[0]} vs raw_adata {raw_adata.n_obs}"

    adata = ad.AnnData(X=matrix.copy(), obs=raw_adata.obs.copy())

    # 공정 비교: 전부 같은 방식으로 처리
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, use_rep="X_pca")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="leiden")

    if label_key in adata.obs.columns:
        y_true = adata.obs[label_key].astype("category")
        y_pred = adata.obs["leiden"].astype("category")
        ari = adjusted_rand_score(y_true, y_pred)
    else:
        ari = None

    if out_png is not None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        if label_key in adata.obs.columns:
            sc.pl.umap(
                adata,
                color=label_key,
                ax=ax,
                show=False
            )
            ax.set_title(f"{title_prefix} (color={label_key}, ARI={ari:.3f})")
        else:
            sc.pl.umap(
                adata,
                color=None,
                ax=ax,
                show=False
            )
            ax.set_title(f"{title_prefix} (no {label_key} in obs)")
        fig.tight_layout()
        fig.savefig(out_png, dpi=300)
        plt.close(fig)

    return adata, ari


def main():
    # --- 2. Raw AnnData / input txt 로드 ---
    if not RAW_H5AD.exists():
        raise FileNotFoundError(f"RAW_H5AD not found: {RAW_H5AD}")
    adata_raw = sc.read_h5ad(RAW_H5AD)
    print(f"[info] Loaded RAW_H5AD: {RAW_H5AD} | shape={adata_raw.shape}")
    print(f"[info] obs columns: {list(adata_raw.obs.columns)}")

    if not SCSTD_INPUT_TXT.exists():
        raise FileNotFoundError(f"SCSTD_INPUT_TXT not found: {SCSTD_INPUT_TXT}")
    raw_matrix = pd.read_csv(SCSTD_INPUT_TXT, header=None, sep=",").values
    print(f"[info] Loaded scSTD_input txt: {SCSTD_INPUT_TXT} | shape={raw_matrix.shape}")

    if raw_matrix.shape[0] != adata_raw.n_obs:
        print(f"[warn] Cell count mismatch: sim1.h5ad n_obs={adata_raw.n_obs}, "
              f"sim1.txt rows={raw_matrix.shape[0]}")

    # --- 3. AE latent / sampled latent 로드 + 통계 ---
    if not AE_EMBED_PATH.exists():
        raise FileNotFoundError(f"AE_EMBED_PATH not found: {AE_EMBED_PATH}")
    ae_emb = np.load(AE_EMBED_PATH)
    print(f"[info] Loaded AE embeddings: {AE_EMBED_PATH} | shape={ae_emb.shape}")
    print(f"       AE latent mean={ae_emb.mean():.4f}, std={ae_emb.std():.4f}")

    if not SAMPLED_LATENT_PATH.exists():
        raise FileNotFoundError(f"SAMPLED_LATENT_PATH not found: {SAMPLED_LATENT_PATH}")
    sampled_latent = np.load(SAMPLED_LATENT_PATH)
    print(f"[info] Loaded sampled latent: {SAMPLED_LATENT_PATH} | shape={sampled_latent.shape}")
    print(f"       sampled latent mean={sampled_latent.mean():.4f}, std={sampled_latent.std():.4f}")

    # --- 4. latent real vs sample UMAP ---
    print("[step] UMAP: real latent vs sampled latent")

    adata_q = ad.AnnData(ae_emb.copy())
    adata_s = ad.AnnData(sampled_latent.copy())
    adata_q.obs["type"] = "real"
    adata_s.obs["type"] = "sample"
    adata_lat = ad.concat([adata_q, adata_s])

    sc.pp.neighbors(adata_lat, use_rep="X")
    sc.tl.umap(adata_lat)

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    sc.pl.umap(adata_lat, color="type", ax=ax, show=False)
    fig.tight_layout()
    latent_umap_path = DEBUG_DIR / "sim1_latent_real_vs_sample_umap.png"
    fig.savefig(latent_umap_path, dpi=300)
    plt.close(fig)
    print(f"[save] latent UMAP saved to: {latent_umap_path}")

    # --- 5. 최종 scSTD-imputed 결과 로드 ---
    if not SCSTD_IMPUTED_PATH.exists():
        raise FileNotFoundError(f"SCSTD_IMPUTED_PATH not found: {SCSTD_IMPUTED_PATH}")
    scstd_imputed = np.loadtxt(SCSTD_IMPUTED_PATH, delimiter=",")
    print(f"[info] Loaded scSTD-imputed result: {SCSTD_IMPUTED_PATH} | shape={scstd_imputed.shape}")

    if scstd_imputed.shape[0] != adata_raw.n_obs:
        print(f"[warn] Cell count mismatch: scSTD_imputed rows={scstd_imputed.shape[0]}, "
              f"adata_raw n_obs={adata_raw.n_obs}")

    # --- 6. Raw vs scSTD-imputed UMAP + ARI 비교 ---
    print("[step] Raw vs scSTD-imputed: UMAP + ARI(Group vs Leiden)")

    raw_adata_proc, ari_raw = build_adata_and_umap(
        matrix=raw_matrix,
        raw_adata=adata_raw,
        label_key="Group",
        title_prefix="Raw (sim1)",
        out_png=DEBUG_DIR / "sim1_raw_umap_Group.png"
    )
    print(f"[metric] Raw ARI(Group vs Leiden) = {ari_raw:.4f}" if ari_raw is not None else "[metric] Raw ARI = None")

    scstd_adata_proc, ari_scstd = build_adata_and_umap(
        matrix=scstd_imputed,
        raw_adata=adata_raw,
        label_key="Group",
        title_prefix="scSTD full pipeline (sim1)",
        out_png=DEBUG_DIR / "sim1_scSTD_umap_Group.png"
    )
    print(f"[metric] scSTD-full ARI(Group vs Leiden) = {ari_scstd:.4f}" if ari_scstd is not None else "[metric] scSTD-full ARI = None")

    print("\n[summary]")
    print(f"  Raw      ARI: {ari_raw}")
    print(f"  scSTD    ARI: {ari_scstd}")
    print(f"\nUMAP png들은 {DEBUG_DIR} 밑에 저장됨.")


if __name__ == "__main__":
    main()


import pandas as pd
import numpy as np
import scanpy as sc
import os

# --- 1. 경로 설정 ---
# 원본 AnnData 객체 파일 경로 (이름을 가져오기 위해 필요)
ORIGINAL_DATA_DIR = f"{BASE_DIR}/splatter/Raw"
# Imputed TXT 파일이 저장된 경로
IMPUTED_RESULT_DIR = f"{BASE_DIR}/splatter/imputed/scSTD/Result"
# 처리할 데이터셋 목록
DATASETS = ["sim1", "sim2", "sim3"]
# 모든 결과를 저장할 딕셔너리
imputed_dataframes = {}

# --- 2. 데이터 처리 루프 ---
print("✅ Imputation 결과 파일 로드 및 이름 복원 시작...")
for dataset_name in DATASETS:
    print(f"\n--- 데이터셋: {dataset_name} 처리 중 ---")
    # 2.1. 경로 조합
    original_h5ad_path = os.path.join(ORIGINAL_DATA_DIR, f"{dataset_name}.h5ad")
    imputed_txt_path = os.path.join(IMPUTED_RESULT_DIR, f"{dataset_name}_imputed.txt")
    # 2.2. 원본 AnnData 로드 (이름 정보 획득용)
    try:
        # scanpy를 사용하여 AnnData 객체 로드
        original_adata = sc.read_h5ad(original_h5ad_path)
        gene_names = original_adata.var_names.tolist() # 유전자 이름 (열 이름)
        cell_names = original_adata.obs_names.tolist() # 세포 이름 (행 인덱스)
        print(f"   -> 원본 데이터 로드 완료. (Cells: {len(cell_names)}, Genes: {len(gene_names)})")
    except FileNotFoundError:
        print(f"   ❌ 오류: 원본 AnnData 파일 '{original_h5ad_path}'을(를) 찾을 수 없습니다. 건너뜁니다.")
        continue
    except Exception as e:
        print(f"   ❌ 오류: 원본 AnnData 로드 중 예상치 못한 오류 발생: {e}")
        continue
    # 2.3. Imputed TXT 파일 로드
    try:
        imputed_df = pd.read_csv(imputed_txt_path, header=None, index_col=None, sep=',')
        # 2.4. 이름 연결 (행 인덱스, 열 이름 설정)
        imputed_df.index = cell_names
        imputed_df.columns = gene_names
        # 2.5. 결과 저장
        imputed_dataframes[dataset_name] = imputed_df
        print(f"   -> Imputed TXT 로드 및 이름 복원 완료. (결과 DataFrame 저장)")
    except FileNotFoundError:
        print(f"   ❌ 오류: Imputed TXT 파일 '{imputed_txt_path}'을(를) 찾을 수 없습니다. 건너뜁니다.")
        continue

print("\n=======================================================")
print("✅ 모든 데이터셋 처리가 완료되었습니다.")
print("   결과는 'imputed_dataframes' 딕셔너리에 저장되어 있습니다.")

# --- 3. 결과 확인 (선택 사항) ---
for name, df in imputed_dataframes.items():
    print(f"\n[결과 확인] {name} DataFrame:")
    print(f"   - 형태: {df.shape}")
    print(f"   - 행 인덱스 (세포 이름) 미리보기: {df.index[:5].tolist()}")
    print(f"   - 열 이름 (유전자 이름) 미리보기: {df.columns[:5].tolist()}")
    print("   - 상위 5x5 데이터:")
    print(df.iloc[:5, :5])


import os
import pandas as pd
import scanpy as sc
import numpy as np
import scipy.sparse  # 👈 [추가] 희소 행렬 확인을 위한 라이브러리

# --- 설정 ---
# os.getcwd() 
os.chdir(f"{BASE_DIR}/Chu_time/Raw")
chu = sc.read_h5ad("chu.h5ad")

adata_list = [chu]
file_names = ['chu']
output_dir = f"{BASE_DIR}/Chu_time/Raw/scIGANs_input"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# --- 변환 및 저장 루프 ---
for adata, name in zip(adata_list, file_names):
    print(f"Processing {name}...")
    
    # 1. Sparse 여부 확인 후 Dense 변환
    # sc.issparse -> scipy.sparse.issparse 로 변경 
    if scipy.sparse.issparse(adata.X):  # 👈 [수정됨]
        data_genes_x_cells = adata.X.T.toarray()
    else:
        data_genes_x_cells = adata.X.T
    
    # 2. DataFrame 생성
    df = pd.DataFrame(
        data_genes_x_cells,
        index=adata.var_names,
        columns=adata.obs_names
    )
    
    # 3. 데이터 저장
    output_path = os.path.join(output_dir, f"{name}.T.tsv")
    df.to_csv(output_path, sep='\t', index=True, index_label='Gene_ID')
    
    # 5. 라벨 변환 (문자열 -> 숫자)
    labels = adata.obs['cell_type'].astype('category').cat.codes
    
    # 6. 라벨 저장
    label_output_path = os.path.join(output_dir, f"{name}_labels.tsv")
    labels.to_csv(
        label_output_path,
        sep='\t',
        index=False,  
        header=False  
    )
    print(f"-> Saved labels to {label_output_path}")

print("\nAll files and labels saved successfully.")

conda activate scIGANs
cd $BENCHMARK_DIR/scIGANs

time CUDA_VISIBLE_DEVICES=0 scIGANs \
    $DATA_DIR/Chu_time/Raw/scIGANs_input/chu.T.tsv \
    -l $DATA_DIR/Chu_time/Raw/scIGANs_input/chu_labels.tsv \
    -o $DATA_DIR/Chu_time/imputed/scIGANs \
    -j chu \
    -e 200


import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

BASE = f"{BASE_DIR}/splatter/Raw"
OUT  = os.path.join(BASE, "scMultiGAN_input")
os.makedirs(OUT, exist_ok=True)

def export_scMultiGAN_inputs(adata: sc.AnnData,
                             basename: str,
                             out_dir: str,
                             label_keys=("group","Group","cell_type","label")):
    """
    저장 결과:
      - {basename}_matrix.tsv  : genes x cells, 숫자만, 헤더=cell IDs, index 없음
      - {basename}_label.txt   : 단일 컬럼 헤더 'x', 1-based 정수 라벨
    """
    # 1) counts 우선 선택
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    M = X.T.toarray() if issparse(X) else X.T  # genes x cells
    
    if (M < 0).any():
        raise ValueError("음수 값이 있습니다. 원본이 카운트인지 확인하세요.")
    M = np.rint(M).astype(np.int64)
    
    # 2) matrix 저장 (헤더=cell IDs, 유전자명 컬럼 절대 추가하지 않음)
    df = pd.DataFrame(M, columns=adata.obs_names)
    mat_path = os.path.join(out_dir, f"{basename}_matrix.tsv")
    df.to_csv(mat_path, sep="\t", index=False, header=True)
    
    # 3) 라벨 저장 (단일 컬럼, 헤더 'x', 1-based)
    for k in label_keys:
        if k in adata.obs:
            lab_key = k; break
    else:
        raise KeyError(f"라벨 컬럼을 찾을 수 없습니다. candidates={label_keys}")
    
    labels_1based = adata.obs[lab_key].astype("category").cat.codes + 1
    lab_path = os.path.join(out_dir, f"{basename}_label.txt")
    pd.DataFrame({"x": labels_1based}).to_csv(lab_path, sep="\t", index=False)
    
    print(f"✅ {basename}: matrix={df.shape} (genes,cells) | "
          f"labels={len(labels_1based)} cells, classes={labels_1based.nunique()} (key='{lab_key}')")
    return mat_path, lab_path


# --- 실행 예시 ---
sim1 = sc.read_h5ad(os.path.join(BASE, "sim1.h5ad"))
export_scMultiGAN_inputs(sim1, "sim1", OUT)
sim2 = sc.read_h5ad(os.path.join(BASE, "sim2.h5ad"))
export_scMultiGAN_inputs(sim2, "sim2", OUT)
sim3 = sc.read_h5ad(os.path.join(BASE, "sim3.h5ad"))
export_scMultiGAN_inputs(sim3, "sim3", OUT)


BASE_RAW=$DATA_DIR/splatter/Raw/scMultiGAN_input
BASE_OUT=$DATA_DIR/splatter/imputed/scMultiGAN

# 데이터셋별 런 디렉토리
RUN1=$BASE_OUT
# 환경 접속
conda activate scMultiGAN_R
cd $BENCHMARK_DIR/scMultiGAN/code
# Data preprocessing process
Rscript generate.data.R \
  --expression_matrix_path "$DATA_DIR/splatter/Raw/scMultiGAN_input/sim2_matrix.tsv" \
  --file_suffix "tsv" \
  --label_file_path "$DATA_DIR/splatter/Raw/scMultiGAN_input/sim2_label.txt"
mv scMultiGAN.csv $RUN1/scMultiGAN.csv
cp $BASE_RAW/sim1_label.txt $RUN1/sim1_label.txt

Rscript generate.data.R \
  --expression_matrix_path "$DATA_DIR/splatter/Raw/scMultiGAN_input/sim2_matrix.tsv" \
  --file_suffix "tsv" \
  --label_file_path "$DATA_DIR/splatter/Raw/scMultiGAN_input/sim2_label.txt"
mv scMultiGAN.csv $RUN1/scMultiGAN.csv
cp $BASE_RAW/sim2_label.txt $RUN1/sim2_label.txt

Rscript generate.data.R \
  --expression_matrix_path "$DATA_DIR/splatter/Raw/scMultiGAN_input/sim1_matrix.tsv" \
  --file_suffix "tsv" \
  --label_file_path "$DATA_DIR/splatter/Raw/scMultiGAN_input/sim3_label.txt"
mv scMultiGAN.csv $RUN1/scMultiGAN.csv
cp $BASE_RAW/sim3_label.txt $RUN1/sim3_label.txt

# ===================================== Train scMultiGAN ===================================== 
conda deactivate
conda activate scMultiGAN
cd $BENCHMARK_DIR/scMultiGAN/code

SECONDS=0
RUN1=$DATA_DIR/splatter/imputed/scMultiGAN/sim1

CUDA_VISIBLE_DEVICES=0 python train_scMultiGAN.py \
  --epoch 200 --batch-size 8 --save-interval 10 \
  --d_file "$RUN1/scMultiGAN.csv" \
  --c_file "$RUN1/sim1_label.txt" \
  --img_size 100 --ncls 8 \
  --output_dir "$RUN1/result" \
  --lr 1e-4 \
  --num-workers 0

elapsed=$SECONDS
printf "⏱ Training time: %02d:%02d:%02d (hh:mm:ss)\n" \
  $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

conda activate scMultiGAN
cd $BENCHMARK_DIR/scMultiGAN/code

SECONDS=0
RUN=$DATA_DIR/splatter/imputed/scMultiGAN/sim2

CUDA_VISIBLE_DEVICES=2 python train_scMultiGAN.py \
  --epoch 150 --batch-size 8 --save-interval 10 \
  --d_file "$RUN/scMultiGAN.csv" \
  --c_file "$RUN/sim2_label.txt" \
  --img_size 100 --ncls 8 \
  --output_dir "$RUN/result" \
  --lr 1e-4 \
  --num-workers 0

elapsed=$SECONDS
printf "⏱ Training time: %02d:%02d:%02d (hh:mm:ss)\n" \
  $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

conda activate scMultiGAN
cd $BENCHMARK_DIR/scMultiGAN/code

SECONDS=0
RUN1=$DATA_DIR/splatter/imputed/scMultiGAN/sim3

CUDA_VISIBLE_DEVICES=3 python train_scMultiGAN.py \
  --epoch 150 --batch-size 8 --save-interval 10 \
  --d_file "$RUN1/scMultiGAN.csv" \
  --c_file "$RUN1/sim3_label.txt" \
  --img_size 100 --ncls 8 \
  --output_dir "$RUN1/result" \
  --lr 1e-4 \
  --num-workers 0

elapsed=$SECONDS
printf "⏱ Training time: %02d:%02d:%02d (hh:mm:ss)\n" \
  $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))    

# ===================================== Train scMultiGAN_impute ===================================== 
SECONDS=0

CUDA_VISIBLE_DEVICES=7 python train_scMultiGAN_impute.py \
  --d_file "$RUN1/scMultiGAN.csv" \
  --c_file "$RUN1/chu_label.txt" \
  --img_size 130 --ncls 6 \
  --output_dir "$RUN1/result" \
  --pretrain "$RUN1/result/model/0099.pth" \
  --imputer_model "$RUN1/result/model/0099.pth"
  --batch-size 8

elapsed=$SECONDS
printf "⏱ Training time: %02d:%02d:%02d (hh:mm:ss)\n" \
  $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

# ===================================== KBDS로 이동 (1) ===================================== 
rsync -avzP --progress -e "ssh -p 22" \
  $DATA_DIR/splatter/Raw/scMultiGAN_input \
  user@server:/home01/user/Datasets/splatter/Raw
rsync -avzP --progress -e "ssh -p 22" \
  $DATA_DIR/splatter/imputed/scMultiGAN \
  user@server:/home01/user/Datasets/splatter/imputed

# ===================================== 또는: KBDS에서 실행 ===================================== 
#!/bin/bash -l
#SBATCH --job-name=sim1_scmultigan
#SBATCH --partition=debug-1gpu
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/home01/user/scMultiGAN/logs/%x-%j.out

cd /home01/user/scMultiGAN
module load compilers/cuda/12.4
module load libraries/nccl/2.21.5

# conda 환경 활성화
conda activate scMultiGAN
BASE_OUT=/home01/user/Datasets/splatter/imputed/scMultiGAN/sim1
RUN=$BASE_OUT

echo "=============================="
echo "Job started at: $(date)"
echo "=============================="

python /home01/user/scMultiGAN/code/train_scMultiGAN.py \
  --epoch 100 --batch-size 32 --save-interval 10 \
  --d_file "$RUN/scMultiGAN.csv" \
  --c_file "$RUN/sim1_label.txt" \
  --img_size 100 --ncls 8 \
  --output_dir "$RUN/result" \
  --lr 1e-4 \
  --num-workers 0

elapsed=$SECONDS
printf "⏱ Training time: %02d:%02d:%02d (hh:mm:ss)\n" \
  $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))

# ===================================== KBDS 에서 받아옴 ===================================== 
rsync -avzP --progress -e "ssh -p 22" \
  $DATA_DIR/Zheng/Raw/scMultiGAN_input \
  user@server:/home01/user/Datasets/Zheng/Raw
rsync -avzP --progress -e "ssh -p 22" \
  user@server:Datasets/Chu_time/imputed/scMultiGAN/imputed_data/scMultiGAN.csv \
  $DATA_DIR/Chu_time/imputed/scMultiGAN/result/

# ===================================== Train scMultiGAN_impute (in KBDS) ===================================== 
#!/bin/bash -l
#SBATCH --job-name=scmultigan_impute
#SBATCH --partition=debug-1gpu
#SBATCH --gres=gpu:1
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/home01/user/scMultiGAN/logs/%x-%j.out

set -euo pipefail

cd /home01/user/scMultiGAN
module load compilers/cuda/12.4
module load libraries/nccl/2.21.5

# conda 환경 활성화
conda activate scMultiGAN
BASE_OUT=/home01/user/Datasets/splatter/imputed/scMultiGAN/sim1
RUN=$BASE_OUT

echo "=============================="
echo "Job started at: $(date)"
echo "=============================="

# training and impute
# --pretrain 경로에 /result/model/ 추가됨 (ls 확인 결과 반영)
python /home01/user/scMultiGAN/code/train_scMultiGAN_impute.py \
  --epoch 100 \
  --batch-size 12 \
  --save-interval 10 \
  --d_file "$RUN/scMultiGAN.csv" \
  --c_file "$RUN/sim1_label.txt" \
  --img_size 100 \
  --ncls 8 \
  --lr 1e-4 \
  --output_dir "$RUN" \
  --pretrain "$RUN/result/model/0099.pth" \
  --imputer_model "$RUN/model_impute/0099.pth" \
  --num_workers 0

elapsed=$SECONDS

echo "=============================="
echo "Job finished at: $(date)"
printf "Total elapsed time: %02d:%02d:%02d (hh:mm:ss)\n" \
  $((elapsed/3600)) $(((elapsed%3600)/60)) $((elapsed%60))
echo "=============================="

# ===================================== KBDS에서 서버로 이동 =====================================
scp user@server:/home01/user/Datasets/splatter/imputed/scMultiGAN/sim1/imputed_data/scMultiGAN.csv $DATA_DIR/splatter/imputed/scMultiGAN/sim1/result/
scp user@server:/home01/user/Datasets/splatter/imputed/scMultiGAN/sim2/imputed_data/scMultiGAN.csv $DATA_DIR/splatter/imputed/scMultiGAN/sim2/result/
scp user@server:/home01/user/Datasets/splatter/imputed/scMultiGAN/sim3/imputed_data/scMultiGAN.csv $DATA_DIR/splatter/imputed/scMultiGAN/sim3/result/

import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import sparse

# ---------------------------------------------------------
# 공통 원본 데이터 경로 설정
# ---------------------------------------------------------
RAW_DIR = f"{BASE_DIR}/splatter/Raw"
sim_list = ["sim1", "sim2", "sim3"]

for sim in sim_list:
    print(f"\n{'='*50}")
    print(f"🚀 Processing dataset: {sim}")
    print(f"{'='*50}")
    
    # ---------------------------------------------------------
    # 동적 경로 설정 (Dynamic Path Generation)
    # ---------------------------------------------------------
    RAW_PATH = os.path.join(RAW_DIR, f"{sim}.h5ad")
    CSV_PATH = ff"{BASE_DIR}/splatter/imputed/scMultiGAN/{sim}/result/scMultiGAN.csv"
    OUTDIR   = ff"{BASE_DIR}/splatter/imputed/scMultiGAN/{sim}/Restored"
    
    os.makedirs(OUTDIR, exist_ok=True)
    
    # ---------------------------------------------------------
    # 1. 원본 데이터 로드 (Scanpy)
    # ---------------------------------------------------------
    print(f"Loading raw data from {RAW_PATH} ...")
    raw_adata = sc.read_h5ad(RAW_PATH)
    
    # ---------------------------------------------------------
    # 2. 복원용 Factor 계산 (Library Size & Scaling Factor)
    # ---------------------------------------------------------
    print("Calculating scaling factors from raw data...")
    
    # (1) Library Size (Total Counts per Cell) - CPM 복원용
    if sparse.issparse(raw_adata.X):
        lib_sizes = raw_adata.X.sum(axis=1).A1  # (n_cells,)
    else:
        lib_sizes = raw_adata.X.sum(axis=1)
    
    # (2) Scaling Factor (Max Log CPM) - 0~1 스케일 복원용
    raw_norm = raw_adata.copy()
    sc.pp.normalize_total(raw_norm, target_sum=1e6) # CPM
    sc.pp.log1p(raw_norm, base=2)                   # Log2(x+1)
    
    if sparse.issparse(raw_norm.X):
        scaling_factors = raw_norm.X.max(axis=1).toarray().flatten() # (n_cells,)
    else:
        scaling_factors = raw_norm.X.max(axis=1)
    
    # ---------------------------------------------------------
    # 3. Imputed 데이터 로드 및 Truncation
    # ---------------------------------------------------------
    print(f"Loading imputed data from {CSV_PATH} ...")
    
    # load (Genes x Cells -> .T makes it Cells x Genes)
    imputed_df = pd.read_csv(CSV_PATH, sep=",", index_col=None, dtype=np.float32).T
    
    # Truncation (genes 방향: columns, padding 제거)
    imputed_df = imputed_df.iloc[:, :raw_adata.n_vars]
    print("Truncated shape:", imputed_df.shape)  # (cells, genes)
    
    # ---------------------------------------------------------
    # 4. 역변환 수행 (0~1 -> LogCPM -> CPM -> Counts)
    # ---------------------------------------------------------
    print("Restoring data to counts...")
    imputed_data = imputed_df.to_numpy(copy=False)
    
    # (1) 0~1 -> Log Scale 복원 (Multiply by Max Log CPM)
    imputed_log = imputed_data * scaling_factors[:, np.newaxis]
    
    # (2) Log Scale -> CPM 복원 (Inverse Log2: 2^x - 1)
    imputed_cpm = np.power(2, imputed_log) - 1
    imputed_cpm = np.maximum(imputed_cpm, 0) # 음수 방지
    
    # (3) CPM -> Counts 복원 (Inverse CPM using Library Size)
    imputed_counts = (imputed_cpm / 1e6) * lib_sizes[:, np.newaxis]
    
    # ---------------------------------------------------------
    # 5. 저장 (Genes x Cells) TSV, no header
    # ---------------------------------------------------------
    out_tsv = os.path.join(OUTDIR, f"{sim}_scMultiGAN.tsv")
    print(f"Saving restored count data to {out_tsv} ...")
    
    # 저장할 때는 다시 Transpose 하여 (Genes, Cells) 형태로 저장
    np.savetxt(out_tsv, imputed_counts.T, delimiter="\t", fmt="%.6g")
    print(f"✅ {sim} saved → {out_tsv}")

print("\n🎉 모든 시뮬레이션 데이터 복원 및 저장이 완료되었습니다!")

import scanpy as sc
import numpy as np
import pandas as pd
import os
import scipy.sparse as sp

# 데이터 로드
in_dir  = f"{BASE_DIR}/splatter/Raw"
out_dir = f"{BASE_DIR}/splatter/Raw/scGNN_input"
os.makedirs(out_dir, exist_ok=True)

sim1 = sc.read_h5ad(os.path.join(in_dir, "sim1.h5ad"))
sim2 = sc.read_h5ad(os.path.join(in_dir, "sim2.h5ad"))
sim3 = sc.read_h5ad(os.path.join(in_dir, "sim3.h5ad"))

def write_scGNN_input_csv(adata, out_path, layer=None):
    """
    adata -> genes x cells CSV (scGNN PreprocessingscGNN.py CSV 모드용)
    - layer=None 이면 .X 사용
    - index = gene (var_names)
    - columns = cell (obs_names)
    """
    if layer is None:
        X = adata.X
    else:
        X = adata.layers[layer]
    # sparse면 dense로
    if sp.issparse(X):
        X = X.toarray()
    # AnnData는 (cells x genes)이므로 전치
    df = pd.DataFrame(
        X.T,
        index=adata.var_names,   # genes
        columns=adata.obs_names  # cells
    )
    df.to_csv(out_path)

# 여기서는 dropout이 들어간 "관측값"을 쓰고 싶다고 가정해서 .X 사용
write_scGNN_input_csv(sim1, os.path.join(out_dir, "sim1_raw_counts.csv"))
write_scGNN_input_csv(sim2, os.path.join(out_dir, "sim2_raw_counts.csv"))
write_scGNN_input_csv(sim3, os.path.join(out_dir, "sim3_raw_counts.csv"))

print("CSV 저장 완료")


python -W ignore PreprocessingscGNN.py \
  --datasetName sim1_raw_counts.csv \
  --datasetDir $DATA_DIR/splatter/Raw/scGNN_input/ \
  --LTMGDir    $DATA_DIR/splatter/Raw/scGNN_input/sim1/ \
  --filetype   CSV \
  --delim      comma \
  --geneSelectnum 2000 \
  --inferLTMGTag

python -W ignore PreprocessingscGNN.py \
  --datasetName sim2_raw_counts.csv \
  --datasetDir $DATA_DIR/splatter/Raw/scGNN_input/ \
  --LTMGDir    $DATA_DIR/splatter/Raw/scGNN_input/sim2/ \
  --filetype   CSV \
  --delim      comma \
  --geneSelectnum 2000 \
  --inferLTMGTag

python -W ignore PreprocessingscGNN.py \
  --datasetName sim3_raw_counts.csv \
  --datasetDir $DATA_DIR/splatter/Raw/scGNN_input/ \
  --LTMGDir    $DATA_DIR/splatter/Raw/scGNN_input/sim3/ \
  --filetype   CSV \
  --delim      comma \
  --geneSelectnum 2000 \
  --inferLTMGTag

# sim1 impute
CUDA_VISIBLE_DEVICES=5 python -W ignore scGNN.py \
  --datasetName sim1 \
  --datasetDir $DATA_DIR/splatter/Raw/scGNN_input/ \
  --LTMGDir    $DATA_DIR/splatter/Raw/scGNN_input/ \
  --outputDir  $DATA_DIR/splatter/imputed/scGNN/sim1/ \
  --ltmgExpressionFile Use_expression.csv \
  --ltmgFile   LTMG_sparse.mtx \
  --nonsparseMode \
  --regulized-type LTMG

CUDA_VISIBLE_DEVICES=6 python -W ignore scGNN.py \
  --datasetName sim2 \
  --datasetDir $DATA_DIR/splatter/Raw/scGNN_input/ \
  --LTMGDir    $DATA_DIR/splatter/Raw/scGNN_input/ \
  --outputDir  $DATA_DIR/splatter/imputed/scGNN/sim2/ \
  --ltmgExpressionFile Use_expression.csv \
  --ltmgFile   LTMG_sparse.mtx \
  --nonsparseMode \
  --regulized-type LTMG

CUDA_VISIBLE_DEVICES=7 python -W ignore scGNN.py \
  --datasetName sim3 \
  --datasetDir $DATA_DIR/splatter/Raw/scGNN_input/ \
  --LTMGDir    $DATA_DIR/splatter/Raw/scGNN_input/ \
  --outputDir  $DATA_DIR/splatter/imputed/scGNN/sim3/ \
  --ltmgExpressionFile Use_expression.csv \
  --ltmgFile   LTMG_sparse.mtx \
  --nonsparseMode \
  --regulized-type LTMG

import scanpy as sc
import pandas as pd

# 1. 파일 읽기
chu = sc.read_h5ad(f"{BASE_DIR}/Chu_time/Raw/chu.h5ad")

# Dense로 변환
if hasattr(chu.X, 'toarray'):
    print("chu.X를 dense로 변환합니다.")
    chu.X = chu.X.toarray()

print("\nDataFrame 생성 시도...")

try:
    df1 = pd.DataFrame(chu.X.T, index=chu.var_names, columns=chu.obs_names)
    print("df1 생성 성공!")
   
except ValueError as e:
    print(f"[에러 발생] DataFrame 생성 실패: {e}")
    print("--- chu 상세 정보 ---")
    print(chu)
    print(f"chu.X.shape: {chu.X.shape}")
    print(f"len(chu.var_names): {len(chu.var_names)}")
    print(f"len(chu.obs_names): {len(chu.obs_names)}")

# 2. (gene x cell) 형태로 변환
#    .X.toarray()를 사용하여 희소 행렬을 조밀한 배열로 변환한 후 .T로 전치
df1 = pd.DataFrame(chu.X.T, index=chu.var_names, columns=chu.obs_names)
# 3. CSV로 저장
df1.to_csv(f"{BASE_DIR}/Chu_time/Raw/DCA_input/chu.csv")

time CUDA_VISIBLE_DEVICES=7 dca \
{BASE_DIR}/Chu_time/Raw/DCA_input/chu.csv \
{BASE_DIR}/Chu_time/imputed/DCA/ \
-s 64,50,64 \
--nocheckcounts

import os
import torch
import scanpy as sc
import scvi
import numpy as np
import pandas as pd
import scipy.sparse as sp
import time  # 👈 [추가] 시간 측정을 위한 모듈

# --- ⏱️ 전체 시작 시간 기록 ---
start_time = time.time()

# --- 환경 설정 및 시드 고정 ---
os.environ["CUDA_VISIBLE_DEVICES"] = "7"
print(f"CUDA available: {torch.cuda.is_available()}")

scvi.settings.seed = 0
torch.manual_seed(0)

# --- 💡 1. scVI 처리 함수 (h5ad 저장으로 변경됨) ---

def run_scvi_imputation(
    adata: sc.AnnData,
    output_path: str,
    layer_key: str = "counts",
    n_latent: int = 10,
    max_epochs: int = 200,
    lr: float = 1e-3,
    library_size: float = 1e4
):
    """
    AnnData 객체를 받아 scVI 모델을 학습시키고,
    정규화된 발현값을 계산하여 h5ad 파일로 저장합니다.
    """
    
    print(f"--- Processing {output_path} ---")
    # 1. AnnData 복사
    adata_copy = adata.copy()
    # 2. [핵심] counts layer 다시 보장
    adata_copy.layers[layer_key] = adata_copy.X.copy()
    # 3. scVI 모델 설정
    scvi.model.SCVI.setup_anndata(
        adata_copy,
        layer=layer_key
    )
    # 3. 모델 초기화 및 학습
    model = scvi.model.SCVI(adata_copy, n_latent=n_latent)
    print(f"Starting model training...")
    model.train(
        max_epochs=max_epochs,
        early_stopping=True,
        plan_kwargs={"lr": lr},
        accelerator="gpu" if torch.cuda.is_available() else "cpu",
        devices=1,
    )
    print("Training complete.")
    
    # 4. 정규화/Imputation된 발현값 추출
    # library_size=1e4로 정규화된 값을 가져옴 (DataFrame 형태 반환됨)
    denoised = model.get_normalized_expression(library_size=library_size)
    denoised_log1p = np.log1p(denoised)  # Log 변환
    
    # 5. [수정됨] AnnData 생성 및 메타데이터 이식
    print("Converting to AnnData...")
    
    # DataFrame을 AnnData로 변환
    imputed_adata = sc.AnnData(denoised_log1p)
    
    # [중요] 원본의 메타데이터(obs, var) 복사
    # 이렇게 해야 나중에 Cell Type이나 Gene Name을 그대로 쓸 수 있음
    imputed_adata.obs = adata.obs.copy()
    imputed_adata.var = adata.var.copy()
    
    # 6. [수정됨] h5ad로 저장
    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    # compression="gzip"을 쓰면 용량이 획기적으로 줄어듭니다.
    imputed_adata.write_h5ad(output_path, compression="gzip")
    print(f"Successfully saved imputed data to {output_path}\n")

# --- 2. 데이터 로드 ---
os.chdir(f"{BASE_DIR}/splatter/Raw")
datasets = {
    "sim1": sc.read_h5ad("sim1.h5ad"),
    "sim2": sc.read_h5ad("sim2.h5ad"),
    "sim3": sc.read_h5ad("sim3.h5ad")
}

# --- 3. 정의된 함수 실행 ---
base_output_dir = f"{BASE_DIR}/splatter/imputed/scVI"
os.makedirs(base_output_dir, exist_ok=True)
# 확장자를 .tsv에서 .h5ad로 변경하여 호출
# --- 3. 반복 실행 (수정됨) ---
for name, adata in datasets.items():
    run_scvi_imputation(
        adata=adata,
        output_path=os.path.join(base_output_dir, f"{name}_scVI.h5ad")
    )

print("--- All processes complete. ---")

# --- ⏱️ 시간 계산 및 출력 ---
end_time = time.time()
elapsed_time = end_time - start_time

hours = int(elapsed_time // 3600)
minutes = int((elapsed_time % 3600) // 60)
seconds = elapsed_time % 60

print(f"\n🚀 총 소요 시간: {hours}시간 {minutes}분 {seconds:.2f}초")

import scanpy as sc
import magic
import os
import time

# 1. 전체 시작 시간 기록
start_time = time.time()

# --- 설정 ---
os.chdir(f'{BASE_DIR}/Chu/Raw')

# 2. 데이터 로드
print("Loading data...")
chu = sc.read_h5ad("chu.h5ad")

# 3. 전처리
print("Preprocessing...")
chu_magic_log = chu.copy()
sc.pp.normalize_total(chu_magic_log, target_sum = 1e4)
sc.pp.log1p(chu_magic_log)

# 4. MAGIC 실행
print("Running MAGIC imputation...")
magic_op_log = magic.MAGIC()
chu_magic_log.X = magic_op_log.fit_transform(chu_magic_log.X)

# 5. 결과 저장 (압축 없음)
output_path = f"{BASE_DIR}/Chu/imputed/MAGIC/chu_MAGIC.h5ad"
output_dir = os.path.dirname(output_path)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print(f"Saving results to {output_path}...")

# 수정됨: compression 옵션 제거 (기본값으로 저장)
chu_magic_log.write(output_path)

print("✅ MAGIC Imputation completed.")

# 6. 소요 시간 계산
end_time = time.time()
elapsed_time = end_time - start_time

hours = int(elapsed_time // 3600)
minutes = int((elapsed_time % 3600) // 60)
seconds = int(elapsed_time % 60)

# 7. 결과 출력
print(f"🚀 총 소요 시간: {hours}시간 {minutes}분 {seconds}초")

import os
import scanpy as sc
import pandas as pd
from scipy.sparse import issparse

# 현재 경로 확인 및 이동 (이미 실행하셨으나, 명시를 위해 포함)
print(f"Current working directory: {os.getcwd()}")
os.chdir(f"{BASE_DIR}/Chu_time/Raw")

# AnnData 파일 로드
chu = sc.read_h5ad("chu.h5ad")
output_dir = f'{BASE_DIR}/Chu_time/Raw/SAVER_input'
os.makedirs(output_dir, exist_ok=True)
# --- 데이터 변환 및 저장 함수 ---
def save_saver_input(adata, name, output_path):
    # 1. Raw Count Data 추출: .X는 (세포 x 유전자) 형태
    # SAVER는 Raw Count를 입력받는 것이 일반적입니다.
    # AnnData의 .X가 이미 정규화(normalized)된 데이터일 수 있으므로,
    # 가능하다면 .raw.X를 사용하는 것이 더 정확하지만, 여기서는 .X를 사용하겠습니다.
    # (splatter 데이터셋의 경우 .X가 일반적으로 Raw Count임)
    X_data = adata.X
    # 2. 형태 변환: (유전자 x 세포)로 전치(Transpose)
    # R에서 Matrix로 쉽게 로드하기 위해 pandas DataFrame으로 변환
    if issparse(X_data):
        # 희소 행렬(Sparse matrix)인 경우 밀집 행렬(Dense matrix)로 변환 후 전치
        X_T = X_data.transpose().toarray()
    else:
        # 밀집 행렬인 경우 바로 전치
        X_T = X_data.T
    # DataFrame 생성 (R에서 Gene/Cell 이름 활용을 위해)
    df = pd.DataFrame(
        X_T,
        index=adata.var_names,  # 유전자 이름 (Rows)
        columns=adata.obs_names  # 세포 이름 (Columns)
    )
    # 3. CSV 파일로 저장
    file_path = os.path.join(output_path, f"{name}.csv")
    df.to_csv(file_path, sep=",", header=True, index=True)
    print(f"Saved {name} data (Genes x Cells) to: {file_path}")
    return file_path

# 데이터셋별로 변환 및 저장
chu_file = save_saver_input(chu, "chu", output_dir)

suppressPackageStartupMessages({
  library(SAVER)
})

# 1. 전체 시작 시간 기록
start_time <- Sys.time()

input_dir   <- f"{BASE_DIR}/Chu_time/Raw/SAVER_input"
imputed_dir <- f"{BASE_DIR}/Chu_time/imputed/SAVER"

impute_saver <- function(data_name, input_dir, output_dir, ncores = 4) {
  message(sprintf("Starting imputation for: %s", data_name))

  file_path <- file.path(input_dir, paste0(data_name, ".csv"))
  if (!file.exists(file_path)) {
    stop(sprintf("Error: File not found at %s", file_path))
  }

  raw_counts <- read.csv(
    file_path,
    header      = TRUE,
    row.names   = 1,
    check.names = FALSE
  )

  # ---- 여기부터 중요 ----
  input_matrix <- as.matrix(raw_counts)
  storage.mode(input_matrix) <- "double"
  # R 4.5 + SAVER class 체크 버그 회피용 강제 패치
  attr(input_matrix, "class") <- "matrix"
  # ---- 여기까지 ----

  message(sprintf("Data dimensions (Genes x Cells): %d x %d",
                  nrow(input_matrix), ncol(input_matrix)))
  message(sprintf("Input matrix class: %s",
                  paste(class(input_matrix), collapse = ", ")))

  saver_result <- SAVER::saver(
    input_matrix,
    do.fast = TRUE,
    ncores  = ncores
  )

  imputed_estimate <- saver_result$estimate

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, paste0(data_name, "_saver_imputed.csv"))
  write.csv(imputed_estimate, output_file, quote = FALSE)

  message(sprintf("Completed and saved results to: %s", output_file))
}

# 2. 작업 실행
impute_saver("chu", input_dir, imputed_dir, ncores = 20)

message("All SAVER imputations completed.")

# 3. 전체 종료 시간 기록 및 소요 시간 계산
end_time <- Sys.time()
time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))

hours <- floor(time_diff / 3600)
minutes <- floor((time_diff %% 3600) / 60)
seconds <- round(time_diff %% 60, 2)

# 4. 결과 출력
message(sprintf("🚀 총 소요 시간: %d시간 %d분 %.2f초", hours, minutes, seconds))

import os
import scanpy as sc
import pandas as pd
from scipy.sparse import issparse

os.chdir(f"{BASE_DIR}/splatter/Raw")
#os.chdir(f"{BASE_DIR}/Zheng/Raw")
#os.chdir(f"{BASE_DIR}/Chu_time/Raw")
# AnnData 파일 로드
sim1, sim2, sim3 = sc.read_h5ad("sim1.h5ad"), sc.read_h5ad("sim2.h5ad"), sc.read_h5ad("sim3.h5ad")
#chu = sc.read_h5ad("chu.h5ad")
#output_dir = f'{BASE_DIR}/Chu_time/Raw/DrImpute_input'
output_dir = f'{BASE_DIR}/splatter/Raw/DrImpute_input'
os.makedirs(output_dir, exist_ok=True)

def save_drimpute_input(adata, name, output_folder):
    """
    AnnData 객체를 로그 정규화한 후, 
    DrImpute 입력 형식(행: 유전자, 열: 세포)의 CSV로 저장합니다.
    """
    # 1. 원본 데이터 보호를 위해 복사
    ad = adata.copy()
    
    # 2. Scanpy 표준 로그 정규화 (10k Scaling + log1p)
    # 이렇게 하면 다른 툴들과 데이터 스케일이 일치하게 됩니다.
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    print(f"[{name}] Log-normalization completed (target_sum=1e4).")
    
    # 3. DrImpute 형식으로 전치 (Gene x Cell)
    X = ad.X.T
   
    # 4. 희소 행렬인 경우 밀집 행렬로 변환
    if issparse(X):
        X = X.toarray()
        
    # 5. 데이터프레임 생성
    df = pd.DataFrame(
        X, 
        index=ad.var_names, 
        columns=ad.obs_names
    )
    
    # 6. CSV 파일로 저장
    file_path = os.path.join(output_folder, f'{name}.csv')
    df.to_csv(file_path, index=True, header=True)
    print(f"Successfully saved normalized {name} to: {file_path}")

# 딕셔너리에 데이터셋을 담아 반복 적용
datasets = {
    "sim1": sim1,
    "sim2": sim2,
    "sim3": sim3
}

for name, data in datasets.items():
    save_drimpute_input(data, name, output_dir)

print("\n🎉 모든 데이터셋의 DrImpute 입력 변환 및 저장이 완료되었습니다.")

# 1. 스크립트 파일 로드
library(DrImpute)

# 전체 작업 시작 시간 기록
total_start_time <- Sys.time()

# 2. 입출력 경로 설정
input_dir <- f'{BASE_DIR}/splatter/Raw/DrImpute_input'
output_dir <- f'{BASE_DIR}/splatter/imputed/DrImpute'

# 3. 분석할 데이터셋 목록
dataset_names <- c("sim1.csv","sim2.csv","sim3.csv")

# 4. DrImpute 실행 함수 정의
run_drimpute <- function(file_name) {
    cat('--------------------------------------------------\n')
    cat(sprintf("Processing file: %s\n", file_name))
    
    # 개별 파일 시작 시간 기록
    file_start_time <- Sys.time()
    
    # 1. 데이터 로드
    input_path <- file.path(input_dir, file_name)
    X_input <- read.csv(input_path, row.names = 1, check.names = FALSE)
    X_matrix <- as.matrix(X_input)
    
    # 2. DrImpute 전처리 적용
    # 🌟 min.expressed.cell = 0 으로 수정 (유전자 누락 방지)
    # 🌟 normalize.by.size.effect = TRUE (로그 변환 수행)
    X_preprocessed <- preprocessSC(
        X = X_matrix,
        min.expressed.gene = 0,
        min.expressed.cell = 0, 
        normalize.by.size.effect = FALSE
    )
    
    # 3. DrImpute 실행
    imputed_data <- DrImpute(X = X_preprocessed)
    
    # 4. 결과 저장
    output_file_name <- gsub(".csv", "_DrImpute.csv", file_name)
    output_path <- file.path(output_dir, output_file_name)
    write.csv(imputed_data, output_path, row.names = TRUE)
    
    # 개별 파일 종료 시간 및 소요 시간 계산
    file_end_time <- Sys.time()
    file_diff <- as.numeric(difftime(file_end_time, file_start_time, units = "secs"))
    
    cat(sprintf("Result saved to: %s\n", output_path))
    cat(sprintf("File processing time: %.2f seconds\n", file_diff))
    cat('--------------------------------------------------\n')
}

# 5. 모든 데이터셋에 함수 적용
for (ds in dataset_names) {
    run_drimpute(ds)
}

# 전체 종료 시간 기록 및 최종 소요 시간 계산
total_end_time <- Sys.time()
total_diff <- as.numeric(difftime(total_end_time, total_start_time, units = "secs"))

# 초 단위를 시/분/초로 변환
hours <- floor(total_diff / 3600)
minutes <- floor((total_diff %% 3600) / 60)
seconds <- round(total_diff %% 60, 2)

cat("\n✅ 모든 DrImpute 보정 작업이 완료되었습니다.\n")
cat(sprintf("🚀 총 소요 시간: %d시 %d분 %.2f초\n", hours, minutes, seconds))

import os
import scanpy as sc
import pandas as pd
from scipy.sparse import issparse

#os.chdir(f"{BASE_DIR}/splatter/Raw")
#os.chdir(f"{BASE_DIR}/5K_PBMC_10X/Raw")
#os.chdir(f"{BASE_DIR}/Zheng/Raw")
os.chdir(f"{BASE_DIR}/Chu_time/Raw")
#output_dir = f"{BASE_DIR}/splatter/Raw/ALRA_input"
output_dir = f"{BASE_DIR}/Chu_time/Raw/ALRA_input"
os.makedirs(output_dir, exist_ok=True)

#sim1, sim2, sim3 = sc.read_h5ad("sim1.h5ad"), sc.read_h5ad("sim2.h5ad"), sc.read_h5ad("sim3.h5ad")
#pbmc = sc.read_h5ad("pbmc.h5ad")
#zheng = sc.read_h5ad("zheng.h5ad")
chu = sc.read_h5ad("chu.h5ad")
def export_for_alra(adata, filepath, layer=None):
    """AnnData → dense matrix → CSV (cells x genes)"""
    if layer is None:
        X = adata.X
    else:
        X = adata.layers[layer]
    
    if issparse(X):
        X = X.toarray()
    
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    df.to_csv(filepath)

# 1) 관측된 drop-out 포함 counts (X 사용) → ALRA에 넣어서 impute 할 입력
#export_for_alra(sim1, os.path.join(output_dir, "sim1.csv"))
#export_for_alra(sim2, os.path.join(output_dir, "sim2.csv"))
#export_for_alra(sim3, os.path.join(output_dir, "sim3.csv"))
#export_for_alra(pbmc, os.path.join(output_dir, "pbmc.csv"))
#export_for_alra(zheng, os.path.join(output_dir, "zheng.csv"))
export_for_alra(chu, os.path.join(output_dir, "chu.csv"))

library(ALRA)   # 여기 안에 alra(), normalize_data() 있음

input_dir  <- f"{BASE_DIR}/Chu/Raw/ALRA_input"
output_dir <- f"{BASE_DIR}/Chu/imputed/ALRA"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(0)  # randomized SVD 재현성용(원하면 숫자 바꿔도 됨)

run_alra_one <- function(name) {
  cat("==========", name, "==========\n")
  in_path  <- file.path(input_dir,  paste0(name, ".csv"))
  out_path <- file.path(output_dir, paste0(name, "_ALRA_imputed.csv"))
  
  rawdata <- as.matrix(read.csv(in_path, row.names = 1, check.names = FALSE))
  cat(sprintf("Loaded %s: %d cells x %d genes\n",
              name, nrow(rawdata), ncol(rawdata)))
  t_start <- Sys.time()
  A_norm  <- normalize_data(rawdata,library.size = 1e4)
  res     <- alra(A_norm)
  t_end   <- Sys.time()
  # ---- 여기부터 깔끔하게 포맷 ----
  elapsed_sec <- as.numeric(difftime(t_end, t_start, units = "secs"))
  h <- floor(elapsed_sec / 3600)
  m <- floor((elapsed_sec %% 3600) / 60)
  s <- round(elapsed_sec %% 60)

  cat(sprintf(
    "Total ALRA time for %s: %02dh %02dm %02ds (%.1f min)\n\n",
    name, h, m, s, elapsed_sec / 60
  ))
  # -------------------------------
  imputed <- res$A_norm_rank_k_cor_sc
  write.csv(imputed, out_path, quote = FALSE)
  cat("Saved imputed matrix to:\n", out_path, "\n\n")
}

# 실제 실행
for (nm in c("chu")) {
  run_alra_one(nm)
}

