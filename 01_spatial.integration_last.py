import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scanpy.external as sce

########## 02 逐个切片展示基因表达量#############
dir = 'D:/Research/Rumen/11_stereo_seq/01_integration/'  # 根据需要设置路径
os.chdir(dir)
outdir = "D:\\Research\\Rumen\\11_stereo_seq\\matrix\\"
sample_names = ["sheep_R1_bw", "sheep_R1_lw", "sheep_R1_ww", "sheep_R1_zw", "horse_R2_zw", "sheep_R2_zw", "sheep_R2_lw", "sheep_R2_bw", "sheep_R2_ww", "horse_R1_zw"]  # 样本名称

for sample in sample_names:
    adata_sp= sc.read("D:/Research/Rumen/11_stereo_seq/matrix/" + sample + ".h5ad")
    # 输出点数
    print(f"Sample: {sample}, Number of points: {adata_sp.shape[0]}")
    adata_sp.layers['counts'] = adata_sp.X.copy()
    sc.pp.normalize_total(adata_sp, target_sum=1e4)
    sc.pp.log1p(adata_sp)
    adata_sp.layers["logcounts"] = adata_sp.X.copy()
    sc.pl.spatial(
        adata_sp,
        color=["KRT6A","TSPYL4", "LUC7L"],  
        spot_size=50,
        size=1.3,
        ncols=1,
        alpha_img=0.5,
        cmap='viridis',  # 可选 colormap  viridis plasma inferno magma
        save=f"_markers_{sample}_granular.pdf"  # 自动保存图像
    )