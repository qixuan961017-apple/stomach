"""
Created on Fri Aug  8 16:58:43 2025
@author: Yahui Zhang
"""

# conda create -n SOAPy_st python=3.9
# pip install SOAPy_st
# conda activate SOAPy_st

import os
import cv2
import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import pandas as pd
import numpy as np
import SOAPy_st as sp
import matplotlib.pyplot as plt
from anndata import AnnData
from tqdm import tqdm
plt.rcParams['pdf.fonttype'] = 42  # 42 表示使用 TrueType 字体

sample = 'horse_R1_zw'
outdir = '/public/home/s20223040710/rumen/11_stereo_seq/06_velocity/'  # 根据需要设置路径
adata = sc.read_h5ad("/public/home/s20223040710/rumen/11_stereo_seq/02_cell2loc_anno/" + sample + "_cell2loc_anno.h5ad")
os.chdir(outdir)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace= True, subset = False)#  subset = True

sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver="arpack")
# sce.pp.harmony_integrate(adata_concat, 'sample')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30) #use_rep="X_pca_harmony",
sc.tl.leiden(adata, resolution=0.6, key_added="leiden")#horse_R2_zw horse_R1_zw 0.6

sc.pl.spatial(adata, img_key="hires", color='leiden', spot_size=50)
plt.savefig(f'{sample}.spatial.niche.png', bbox_inches='tight', dpi=500)

mask = sp.tl.get_mask_from_domain(adata, clusters='4', KSize=51, cluster_key='leiden')

# plt.imshow(mask, cmap='gray')
plt.imshow(mask, cmap='gray', vmin=0, vmax=1)  # plt只能显示0到1之间的值，将0-255进行压缩
plt.savefig(f'{sample}.spatial.mask.png', bbox_inches='tight', dpi=500)


###Regression
sp.tl.spatial_tendency(adata,mask,method= 'poly', radius=1000,location='out',frac=5)
# sp.tl.spatial_tendency(adata,mask,method= 'poly',gene_name = 'LAP3',radius=10,location='out')
sp.pl.show_tendency(adata, gene_name = ['KRT6A','TSPYL4','LUC7L'], show=True)
plt.savefig(f'{sample}.spatial.gene.trajectory.pdf', bbox_inches='tight', dpi=500)

####细胞轨迹
adata.obsm['q05_cell_abundance_w_sf'] = adata.obsm['q05_cell_abundance_w_sf'].rename( columns=lambda col: col.replace('q05cell_abundance_w_sf_', ''))
Spatial = AnnData(X = adata.obsm['q05_cell_abundance_w_sf'],obs = adata.obs,obsm = adata.obsm ,uns = adata.uns)
sp.tl.spatial_tendency(Spatial,mask,method='poly', radius=1000,location='out',frac=5)  #mask1

for cell in Spatial.var.index :
    sp.pl.show_tendency(adata, gene_name = cell, show=True)
    plt.savefig(f'{sample}.spatial.{cell}.trajectory.png', bbox_inches='tight', dpi=500)
    
    
sp.pl.show_tendency(adata, gene_name = Spatial.var.index, show=True)
plt.legend(title='celltype', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title('cell trajectory')
plt.ylabel('cell compositions')
plt.tight_layout()
plt.savefig(f'{sample}.spatial.cell.trajectory.pdf', bbox_inches='tight', dpi=500)


