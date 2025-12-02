# -*- coding: utf-8 -*-
import scanpy as sc
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import numpy as np
import matplotlib.image as mpimg
from matplotlib.patches import Patch


outdir = 'D:/Research/Rumen/11_stereo_seq/02_cell2loc_anno/'  # 根据需要设置路径
os.chdir(outdir)

sample_names = ["sheep_R1_zw", "horse_R2_zw", "sheep_R2_zw", "horse_R1_zw","sheep_R1_bw", "sheep_R1_lw", "sheep_R1_ww", "sheep_R2_lw", "sheep_R2_bw", "sheep_R2_ww"]  # 样本名称
color_dict = {
    "q05cell_abundance_w_sf_Basal cell": "#f8dfe2",
    "q05cell_abundance_w_sf_Neck cell": "#fca7ac",
    "q05cell_abundance_w_sf_Parietal cell": "#f3868c",
    "q05cell_abundance_w_sf_Pit cell": "#f12929",
    "q05cell_abundance_w_sf_Spinous cell": "#7f1114",
    "q05cell_abundance_w_sf_Tuft cell": "#651b1c",
    "q05cell_abundance_w_sf_Enteroendocrine cell": "#f77904",
    "q05cell_abundance_w_sf_Capillary endothelial": "#fdf8c9",
    "q05cell_abundance_w_sf_Lymphatic endothelial": "#f3cb45",	
    "q05cell_abundance_w_sf_Fibroblast": "#b8ce83",
    "q05cell_abundance_w_sf_Smooth muscle cell": "#91ae6c",
    "q05cell_abundance_w_sf_Pericytes": "#BFE1E5",
    "q05cell_abundance_w_sf_Proliferative T cell": "#95D0E2",
    "q05cell_abundance_w_sf_CD4_T cell": "#5AC3E1",
    "q05cell_abundance_w_sf_CD8_T cell": "#479DCB",
    "q05cell_abundance_w_sf_Activated T cell": "#3D8AA9",
    "q05cell_abundance_w_sf_Plasma_B": "#8AA2DD",
    "q05cell_abundance_w_sf_B cell": "#5F7DC7",
    "q05cell_abundance_w_sf_Mast cell": "#C7D1EE",
    "q05cell_abundance_w_sf_Chief cell": "#D3D3D3",
    "q05cell_abundance_w_sf_Marcphages": "#A89ADF",    
    
    "q05cell_abundance_w_sf_Schwann cell": "#000000",
    "q05cell_abundance_w_sf_Monocytes": "#907FD3",
    "q05cell_abundance_w_sf_Artery endothelial": "#f5e091",
    "q05cell_abundance_w_sf_Vein endothelial": "#f9ef9a",
}#颜色板

# 循环绘图
for sample in sample_names:
    print(f"Processing {sample}...")
    cell_type_predictions = pd.read_csv(f"{sample}_predicted_cell_types.csv", index_col=0)
    adata_sp = sc.read(f"{sample}_cell2loc_anno.h5ad")
    coords = adata_sp.obsm['spatial']
    x, y = coords[:, 0], coords[:, 1]

    # 计算图像比例
    x_span = x.max() - x.min()
    y_span = y.max() - y.min()
    aspect_ratio = y_span / x_span if x_span > 0 else 1.0

    fig, ax = plt.subplots(figsize=(12, 10 * aspect_ratio))

    for i, spot in enumerate(cell_type_predictions.index):
        proportions = cell_type_predictions.loc[spot]
        proportions = proportions[proportions > 0]

        # 👉 只保留你手动定义了颜色的类型
        proportions = proportions[proportions.index.isin(color_dict.keys())]
        if proportions.sum() == 0:
            continue

        # 归一化 + 最多保留 4 类
        proportions = proportions.sort_values(ascending=False).head(4)
        proportions = proportions / proportions.sum()
        spot_colors = [color_dict[ct] for ct in proportions.index]

        axins = inset_axes(ax, width=0.045, height=0.045, loc='center',
                           bbox_to_anchor=(x[i], y[i]),
                           bbox_transform=ax.transData,
                           borderpad=0)
        axins.pie(proportions, labels=None, colors=spot_colors)
        axins.set_aspect('equal')
        axins.axis('off')

    # 设置主图
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.max(), y.min())
    ax.set_aspect('equal')
    ax.axis('off')

    # 图例（只显示你手动设定的类型）
    legend_elements = [Patch(facecolor=color_dict[ct], label=ct) for ct in color_dict.keys()]
    plt.subplots_adjust(right=0.82)
    ax.legend(handles=legend_elements, loc='center left',
              bbox_to_anchor=(1.02, 0.5), title="Cell Types")

    plt.savefig(f"{sample}_cell2location_piecharts.pdf", format="pdf", bbox_inches="tight")
    plt.close()

print("✅ 所有样本绘图完成。")