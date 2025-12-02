# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 15:47:54 2025
@author: Yahui Zhang
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pymn

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'

#outdir = '/public/home/s20223040710/rumen/11_stereo_seq/08_pymn'  # 根据需要设置路径
outdir = 'D:/Research/Rumen/11_stereo_seq/08_pymn/'  # 根据需要设置路径
os.chdir(outdir)

###################02 整合样本做相关性（绵羊前胃、绵羊皱胃、马胃）#######################
sample_names = ["sheep_R1_bw", "sheep_R1_lw", "sheep_R1_ww", "sheep_R2_lw", "sheep_R2_bw", "sheep_R2_ww"]
title = "sheep_qw"
sample_names = ["sheep_R1_zw","sheep_R2_zw"]
title = "sheep_zw"
sample_names = ["horse_R2_zw","horse_R1_zw"]
title = "horse_stomach"

sp_list = []
for sample in sample_names:    
    adata_sp = sc.read("D:/Research/Rumen/11_stereo_seq/02_cell2loc_anno/" + sample + "_cell2loc_anno.h5ad")
    print(f'sp基因的数量是：{len(adata_sp.var)}')
    # 读取细胞类型预测数据
    cell_type_predictions = pd.read_csv('D:/Research/Rumen/11_stereo_seq/02_cell2loc_anno/' + sample + '_predicted_cell_types.csv', index_col=0)
    predicted_cell_types_idx = np.argmax(cell_type_predictions.to_numpy(), axis=1)
    cell_type_names = cell_type_predictions.columns
    predicted_cell_types = cell_type_names[predicted_cell_types_idx]
    # 将预测的细胞类型添加到 adata_sp.obs 中
    adata_sp.obs['celltype'] = predicted_cell_types
    adata_sp.obs['celltype'] = adata_sp.obs['celltype'].str.replace('^q05cell_abundance_w_sf_', '', regex=True)
    adata_sp.obs['celltype'] = adata_sp.obs['celltype'].astype(str)
    sp_list.append(adata_sp)
    
adata_sp_all = sc.concat(sp_list, label='sample', join='inner', index_unique=None)

adata_sc = sc.read_h5ad('D:/Research/Rumen/10_scdata_add/zw.h5ad')
adata_sc.X = adata_sc.layers['counts']
adata_sc.obs['celltype'] = adata_sc.obs['major_celltype'].astype(str)
        
# 合并单细胞数据和空转数据
adata_combined = adata_sc.concatenate(adata_sp_all, batch_key="source", batch_categories=["sc", "sp"], index_unique=None, join='outer')
sc.pp.normalize_total(adata_combined, inplace=True)
sc.pp.log1p(adata_combined)
    
#sc.pp.highly_variable_genes(adata_combined, flavor='seurat', n_top_genes=3000)  # 选择最具变异的3000个基因
adata_combined.obs["source"] = adata_combined.obs["source"].astype(str)
adata_combined.obs["celltype"] = adata_combined.obs["celltype"].astype(str)

# 进行MetaNeighborUS分析
pymn.variableGenes(adata_combined, study_col='source')
len(adata_combined.var[adata_combined.var['highly_variable']])
pymn.MetaNeighborUS(adata_combined, study_col='source', ct_col='celltype', fast_version=True)
    
# 提取AUROC矩阵
auroc_combined = adata_combined.uns['MetaNeighborUS']

# 定义行和列的顺序
row_order = ["sp|Basal cell", "sp|Neck cell", "sp|Parietal cell", "sp|Pit cell", "sp|Spinous cell", "sp|Tuft cell", 
             "sp|Chief cell", "sp|Enteroendocrine cell", "sp|Capillary endothelial", "sp|Lymphatic endothelial", 
             "sp|Artery endothelial", "sp|Vein endothelial", "sp|Fibroblast", "sp|Smooth muscle cell", "sp|Pericytes", 
             "sp|Proliferative T cell", "sp|CD4_T cell", "sp|CD8_T cell", "sp|Activated T cell", "sp|Plasma_B", "sp|B cell", 
             "sp|Mast cell", "sp|Marcphages", "sp|Monocytes"]
col_order = ["sc|Basal cell", "sc|Neck cell", "sc|Parietal cell", "sc|Pit cell", "sc|Spinous cell", "sc|Tuft cell", 
             "sc|Chief cell", "sc|Enteroendocrine cell", "sc|Capillary endothelial", "sc|Lymphatic endothelial", 
             "sc|Artery endothelial", "sc|Vein endothelial", "sc|Fibroblast", "sc|Smooth muscle cell", "sc|Pericytes", 
             "sc|Proliferative T cell", "sc|CD4_T cell", "sc|CD8_T cell", "sc|Activated T cell", "sc|Plasma_B", "sc|B cell", 
             "sc|Mast cell", "sc|Marcphages", "sc|Monocytes"]


# 过滤有效的行和列
valid_rows = [r for r in row_order if r in auroc_combined.index]
valid_cols = [c for c in col_order if c in auroc_combined.columns]
subset = auroc_combined.loc[valid_rows, valid_cols]
    
# 可视化AUROC矩阵
plt.figure(figsize=(9, 8))
plt.title(title)
sns.heatmap(subset, cbar_kws={'label': 'AUROC'}, cmap="coolwarm", xticklabels=True, yticklabels=True)
plt.tight_layout()
plt.savefig(title + '_auroc_heatmap.pdf', format='pdf')
plt.show()
