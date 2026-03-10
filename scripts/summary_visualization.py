import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb


# 실제 성능 값 데이터
performance_data = {
    'tool': ['RAW', 'ALRA', 'DrImpute', 'MAGIC', 'SAVER', 'scVI', 'DCA', 'scIGANs', 'scMultiGAN', 'scSTD', 'scIDPMs'],
    # Accuracy - Value Reconstruction
    'cell_PCC': [0.853169, 0.865347, 0.940459, 0.92637, 0.885462, 0.90796, 0.923269, 0.841268, 0.861064, 0.920518, 0.864485],
    'Pseudobulk_AUPRC': [0.569347, 0.942346, 0.956207, 0.95637, 0.957048, 0.93844, 0.955999, 0.580193, 0.62341, 0.955355, 0.838376],
    # Accuracy - Signal Recovery
    'AUPRC (cell-level)': [0.721498, 0.853048, 0.936961, 0.923502, 0.938693, 0.907618, 0.91944, 0.770415, 0.762126, 0.920289, 0.829273],
    'logRMSE': [0.356246, 0.42077, 0.22687, 0.396688, 0.31091, 0.26065, 0.250991, 0.336821, 0.322051, 0.256507, 0.320118],
    # Structural Preservation
    'ARI': [0.8415833333, 0.80123, 0.9410843333, 0.6766606667, 0.9003293333, 0.66623, 0.3890303333, 0.969693, 0.8745813333, 0.37414, 0.888385],
    'NMI': [0.8479960833, 0.8296315, 0.9567903333, 0.7274011667, 0.8936123333, 0.72928475, 0.4059188333, 0.969864, 0.9013020833, 0.4501505, 0.849436],
    'Purity': [0.9198958333, 0.905007, 0.9767170833, 0.8277151667, 0.9415823333, 0.8280217, 0.5963575833, 0.98512325, 0.9205453333, 0.617835, 0.91329625],
    'SCC_FC': [0.225814, 0.504338, 0.536459, 0.513973, 0.456138, 0.432284, 0.348096, 0.153445, 0.193296, 0.346732, 0.427714],
    # Biological Consistency - Marker Gene
    'Specificity': [0.4473275, 0.769098, 0.827087, 0.8313805, 0.7913305, 0.806215, 0.7596595, 0.739585, 0.5928635, 0.792473, 0.4852495],
    # Biological Consistency - DEG-based
    'Coexpression': [0.2964735, 0.8695215, 0.4978385, 0.8923175, 0.7490145, 0.789516, 0.79273, 0.296407, 0.0457885, 0.8000175, -0.0120835],
    'DEG F1': [0.4763328351, 0.6519130112, 0.6732632857, 0.6749319913, 0.6682266699, 0.5253251701, 0.6602438627, 0.8895641204, 0.6251975754, 0.4492403551, 0.5046133886],
    'DEG SCC':[0.7078948054,0.6519130112,0.6732632857,0.6749319913,0.6682266699,0.5253251701,0.6602438627,0.08956412043,0.6251975754,0.4492403551,0.504613886],
    'DEG AUPRC': [0.5092558364, 0.4304653905, 0.5022200229, 0.4938014007, 0.4896726608, 0.4778683108, 0.4852531208, 0.4565090876, 0.4638295955, 0.4608688284, 0.4522938148]
}

# 2. DataFrame 생성
df = pd.DataFrame(performance_data)

# 3. 순위 데이터 생성을 위한 함수
def generate_rank_data(perf_df):
    rank_df = pd.DataFrame()
    rank_df['tool'] = perf_df['tool']
    
    # 순위를 매길 컬럼들 (tool 제외)
    metrics = [col for col in perf_df.columns if col != 'tool']
    
    for metric in metrics:
        # logRMSE는 낮을수록 순위가 높음 (ascending=True)
        # 그 외 대부분의 지표는 높을수록 순위가 높음 (ascending=False)
        is_ascending = True if 'RMSE' in metric else False
        
        # rank(method='min')을 사용하여 동일 값일 경우 최소 순위 부여
        # .astype(int)로 정수형 변환
        rank_df[metric] = perf_df[metric].rank(ascending=is_ascending, method='min', descending=not is_ascending).astype(int)
        
        # 참고: pandas의 rank(ascending=False)는 값이 클수록 1등
        # logRMSE처럼 값이 작아야 1등인 경우 ascending=True 사용
        if is_ascending:
            rank_df[metric] = perf_df[metric].rank(ascending=True, method='min').astype(int)
        else:
            rank_df[metric] = perf_df[metric].rank(ascending=False, method='min').astype(int)
            
    return rank_df.to_dict(orient='list')

# 4. 결과 출력
rank_data = generate_rank_data(df)

# 확인을 위해 예쁘게 출력
import pprint

# Set the base directory for datasets
# Change this path to point to your local data directory
BASE_DIR = "./data"
pprint.pprint(rank_data)

df_perf = pd.DataFrame(performance_data)
df_rank = pd.DataFrame(rank_data)

# Min-Max Scaling 함수
def min_max_scale(series):
    min_val = series.min()
    max_val = series.max()
    if max_val == min_val:
        return pd.Series([0.5] * len(series))
    return (series - min_val) / (max_val - min_val)

# 1. 메트릭 순서 (performance_data의 key와 정확히 일치해야 함)
metrics = [
    'cell_PCC', 'Pseudobulk_AUPRC', 'AUPRC (cell-level)', 'logRMSE', # Accuracy
    'ARI', 'NMI', 'Purity', 'SCC_FC',                               # Structural
    'Specificity', 'Coexpression', 'DEG F1', 'DEG SCC', 'DEG AUPRC'  # Biological
]

# 2. 시각화용 레이블 (그래프 하단에 표시될 이름)
metric_labels = {
    'cell_PCC': 'cell PCC',
    'Pseudobulk_AUPRC': 'Pseudobulk\nAUPRC',
    'AUPRC (cell-level)': 'cell-level\nAUPRC',
    'logRMSE': 'logRMSE',
    'ARI': 'ARI',
    'NMI': 'NMI',
    'Purity': 'Purity',
    'SCC_FC': 'SCC-FC',
    'Specificity': 'Specificity',
    'Coexpression': 'Coexpression',
    'DEG F1': 'F1',
    'DEG SCC': 'SCC',
    'DEG AUPRC': 'AUPRC'
}

# 3. 카테고리 매핑 (performance_data의 key 기준)
category_map = {
    'cell_PCC': 'Accuracy',
    'Pseudobulk_AUPRC': 'Accuracy',
    'AUPRC (cell-level)': 'Accuracy',
    'logRMSE': 'Accuracy',
    'ARI': 'Structural',
    'NMI': 'Structural',
    'Purity': 'Structural',
    'SCC_FC': 'Structural',
    'Specificity': 'Biological',
    'Coexpression': 'Biological',
    'DEG F1': 'Biological',
    'DEG SCC': 'Biological',
    'DEG AUPRC': 'Biological'
}

# 4. 대분류별 기본 색상 (기존 유지)
colors = {
    'Accuracy': '#5e3c99',   # 보라색
    'Structural': '#e66101', # 주황색
    'Biological': '#2b83ba'  # 파란색
}

# 각 메트릭별로 Min-Max Scaling
df_normalized = df_perf.copy()
for metric in metrics:
    df_normalized[metric] = min_max_scale(df_perf[metric])

# Figure 생성
fig, axes = plt.subplots(len(df_perf), len(metrics), figsize=(16, 8), 
                         gridspec_kw={'hspace': 0.025, 'wspace': 0.15})

fig.patch.set_alpha(0)

# 각 셀에 바 그래프 그리기
for j, metric in enumerate(metrics):
    for i, tool in enumerate(df_perf['tool']):
        ax = axes[i, j]
        
        # 정규화된 성능 값 (바의 길이)
        normalized_value = df_normalized.loc[df_normalized['tool'] == tool, metric].values[0]
        
        # 순위 (색상)
        rank = df_rank.loc[df_rank['tool'] == tool, metric].values[0]
        
        # 순위를 0-1로 정규화 (1등=0.0 흰색, 11등=1.0 진한색)
        rank_normalized = (rank - 1) / 10
        
        # 대분류에 따른 기본 색상
        base_color = colors[category_map[metric]]
        rgb = to_rgb(base_color)
        
        # 흰색에서 base_color로 그라데이션
        final_color = tuple(1 - rank_normalized * (1 - rgb[k]) for k in range(3))
        
        # 바 그래프 그리기 (height 0.7로 증가)
        ax.barh(0, normalized_value, height=0.7, color=final_color, 
                edgecolor='gray', linewidth=0.5)
        
        # Top 3는 순위 표시
        if rank <= 3:
            ax.text(-0.08, 0, str(int(rank)), ha='center', va='center', 
                   fontsize=13, fontweight='normal', color='black')
        
        # 축 설정
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 0.5)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # 모든 테두리 제거
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        # 맨 왼쪽 컬럼에만 tool 이름 표시
        if j == 0:
            ax.text(-0.35, 0, tool, ha='right', va='center', fontsize=14)
        
        # 맨 아래 행에만 레이블 표시
        if i == len(df_perf) - 1:
            ax.text(0.5, -1.2, metric_labels[metric], ha='right', va='top', 
                   fontsize=14, rotation=45, transform=ax.transData)

# 여백 조정
plt.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.08)

# 저장
output_path = f'{BASE_DIR}/benchmark_bargraph_normalized.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=True, facecolor='none')
plt.close()

print(f"Bar graph saved to: {output_path}")
print("\nNormalized value range check:")
for metric in metrics:
    print(f"{metric}: {df_normalized[metric].min():.3f} - {df_normalized[metric].max():.3f}")