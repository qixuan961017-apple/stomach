# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

####################01 单细胞共定位分析作图 R1和R2求均值############
# outdir = '/public/home/s20223040710/rumen/11_stereo_seq/03_colocalization/'  # 根据需要设置路径
outdir = 'D:/Research/Rumen/11_stereo_seq/03_colocalization/'  # 根据需要设置路径
os.chdir(outdir)

species = "sheep" #horse  sheep
tissue = "bw"  #lw ww bw zw
view = "para" #intra 和  和 juxta 和 para
R1_data = pd.read_csv(species +'_R1_'+ tissue+'_misty_data_' + view + '.csv')
R2_data = pd.read_csv(species +'_R2_'+ tissue+'_misty_data_' + view + '.csv')

# 查看数据的前几行，确保数据被正确读取
print(R1_data)
print(R2_data)

# 合并两个数据框，确保按目标（target）和预测（predictor）进行合并
merged_data = pd.merge(R1_data, R2_data, on=['target', 'predictor', 'view'], suffixes=('_R1', '_R2'))
# 计算每对细胞类型的均值（平均共定位分数）
merged_data['mean_importance'] = merged_data[['importances_R1', 'importances_R2']].mean(axis=1)
#print(merged_data[['target', 'predictor', 'mean_importance']].head())
# 创建一个数据透视表，目标和预测变量作为行和列，均值作为数值
pivot_data = merged_data.pivot_table(index='target', columns='predictor', values='mean_importance')
# 对行索引（target）按升序排序，对列索引（predictor）按降序排序
pivot_data = pivot_data.sort_index(axis=0, ascending=True)  # 排序行（target）
pivot_data = pivot_data.sort_index(axis=1, ascending=False)  # 排序列（predictor）
pivot_data.columns = pivot_data.columns.str.replace('q05cell_abundance_w_sf_', '', regex=False)
pivot_data.index = pivot_data.index.str.replace('q05cell_abundance_w_sf_', '', regex=False)

#####绘制聚类的热图，对角线换成0######
pivot_filled = pivot_data.fillna(0)
plt.figure(figsize=(9, 12))
sns.clustermap(pivot_filled, cmap='Blues', metric='euclidean', method='average')
plt.title(f'Colocalization {species}_{tissue}')
plt.xlabel('Predictor')
plt.ylabel('Target')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'{species}_{tissue}_' + view + '_colocalization_heatmap_cluster.pdf', format='pdf')
plt.show()
