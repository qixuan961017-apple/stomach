# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 15:24:15 2025
@author: Yahui Zhang
"""
####https://liana-py.readthedocs.io/en/latest/notebooks/misty.html\

# conda activate stereo
import os
import scanpy as sc
import decoupler as dc
import plotnine as p9
import liana as li
from liana.method import MistyData, genericMistyData, lrMistyData
from liana.method.sp import RandomForestModel, LinearModel, RobustLinearModel

outdir = 'D:/Research/Rumen/11_stereo_seq/03_colocalization/'  # 根据需要设置路径
os.chdir(outdir)

sample_names = ["horse_R1_zw","sheep_R1_ww", "sheep_R1_bw","sheep_R1_lw","horse_R2_zw", "sheep_R1_zw", "sheep_R2_zw", "sheep_R2_lw", "sheep_R2_bw", "sheep_R2_ww", ]  # 
for sample in sample_names:
    adata = sc.read("D:/Research/Rumen/11_stereo_seq/02_cell2loc_anno/" + sample + "_cell2loc_anno.h5ad")
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["logcounts"] = adata.X.copy()
    comps = li.ut.obsm_to_adata(adata, 'q05_cell_abundance_w_sf')
    #### 01 cell cell 不同细胞类型的共定位
    misty = genericMistyData(intra=comps, extra=comps, cutoff=0.05, bandwidth=200, n_neighs=6) #不同的平台是不一样的，分析临近几个点呢 200是通用的
    misty(model=RandomForestModel, n_jobs=10, verbose = True)
    #view= juxta 是点间，主图放了点间的 可以抽取自己想看的细胞类型，在q05_cell_abundance_w_sf里面就抽取，还有 intra 和  和 juxta 和 para
    # intrinsic (intraview)：内部视角，基于单个spot内部的标记表达，relate the expression of other markers to a specifc marker of interest within the same location；
    # local niche view (juxtaview)：邻近视角，通过构建图并计算每个空间单元周围邻居的标记表达总和来定义，relates the expression from the immediate neighborhood of a cell to the observed expression within that cell；
    # the broader, tissue view (paraview)：旁观视角，通过加权总和计算所有空间单元的标记表达，权重可以根据距离函数（如高斯、指数、线性或常数）进行调整，relates the expression of markers measured in cells within a radius around a given cell；
    view='juxta'
    fig = li.pl.interactions(misty, view=view, return_fig=True,figure_size=(10,10)) + p9.scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + p9.ggtitle(f"Cell Type Colocalization  {sample}") +p9.theme(plot_title=p9.element_text(size=12))
    fig.data.to_csv(outdir + '/' + sample + '_misty_data_'+ view +'.csv', index=True)  # index=False 可以避免保存行号
    fig.save(outdir + '/' + sample + '.spatial.misty.celltype.'+ view +'.colocalization.pdf',bbox_inches = 'tight',dpi = 500)
