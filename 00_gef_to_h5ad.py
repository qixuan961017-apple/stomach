# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 09:32:06 2024
@author: Yahui Zhang
"""
import stereo as st
import warnings
warnings.filterwarnings('ignore')
# %config InlineBackend.figure_format = 'retina'

##### 02 循环处理多个文件 #####
dir = 'D:/Research/Rumen/11_stereo_seq/matrix/'
# 遍历并读取每个文件的 meta 信息
with open('rename.txt', encoding='utf-8') as f:
    header = next(f)  # 跳过表头
    for line in f:
        line = line.strip()
        if not line:
            continue
        sid, species, treat, tissue, original = line.split('\t')
        
        data = st.io.read_gef(file_path= sid + '.tissue.gef', bin_size=50)  # 次数，分辨率发生变化，bin_type:bins, bin_size=50
        data.cells.cell_name, data.genes.gene_name
        #QC 计算每个细胞（严格来讲是 bin，为了方便，称为细胞）的 count 数量、基因数量以及线粒体比例
        #data.var['mt'] = data.var_names.str.upper().str.startswith('MT-') 
        data.tl.cal_qc()
        data.plt.violin()
        data.plt.spatial_scatter(cells_key=['total_counts', 'n_genes_by_counts','pct_counts_mt'], dot_size=None, palette='rainbow')
        data.plt.genes_count()
        data
        #过滤细胞和基因
        data = data.tl.filter_cells(
                min_gene=3,
                pct_counts_mt=20
                )
        data = data.tl.filter_genes(min_cell=1)  # 这里在实际分析中，推荐使用更高的参数;这里尊重教程
        data
        #Scanpy
        adata = st.io.stereo_to_anndata(data,flavor='scanpy')
        adata.obs['species'] = species
        adata.obs['treat'] = treat
        adata.obs['tissue'] = tissue
        sample = species+'_'+treat+'_'+tissue
        adata.obs['sample'] = sample
        adata.write(dir + sample + '.h5ad')
