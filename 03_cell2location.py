# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 23:15:56 2025
@author: Yahui Zhang
"""

###需要用高算平台的GPU进行计算
###https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html

defined_cols = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080']
import os
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import matplotlib as mpl
import pandas as pd
import cell2location
from cell2location.utils.filtering import filter_genes
from matplotlib import rcParams
from cell2location.models import RegressionModel
rcParams['pdf.fonttype'] = 42

outdir = '/public/home/s20223040710/rumen/11_stereo_seq/'  # 根据需要设置路径
os.chdir(outdir)
print("新的工作目录:", os.getcwd())

#####01. 数据预处理 #########
sample = "sheep_R1_ww" #horse_R1_zw
adata_sp = sc.read_h5ad("matrix/"+ sample +".h5ad")
#adata_sp.var['SYMBOL'] = adata_sp.var_names
print(adata_sp.var_names.duplicated().sum()) #用gene symbol之后有很多基因重复项

# find mitochondria-encoded (MT) genes
adata_sp.var['MT_gene'] = [gene.startswith('KEF53-') for gene in adata_sp.var_names]
# remove MT genes for spatial mapping (keeping their counts in the object)
adata_sp.obsm['MT'] = adata_sp[:, adata_sp.var['MT_gene'].values].X.toarray()
adata_sp = adata_sp[:, ~adata_sp.var['MT_gene'].values]

####scrna
adata_ref = sc.read("ww.h5ad")
#需要用原始矩阵
adata_ref.X = adata_ref.layers['counts']


# 调用示例 plot_and_save_QC(mod=sc_model, save_dir="./figures", prefix="mouse_brain")
def plot_and_save_QC(mod, summary_name="means", use_n_obs=1000, scale_average_detection=True,
                     save_dir="./", prefix="QC"):
    """
    绘制并保存 cell2location RegressionModel 的两个 QC 图：
    1. 重建准确性（Reconstruction Accuracy）
    2. 表达 signature vs 聚类平均表达（Signature vs Cluster Mean）
    参数：
        mod: 训练好的 RegressionModel 模型
        summary_name: 使用的 posterior 统计量名称，如 "means", "q05"
        use_n_obs: 子抽样数量用于重建准确性评估
        scale_average_detection: 是否乘以 detection_y_c 平均值
        save_dir: 图像保存路径
        prefix: 图像文件名前缀
    """
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    # ------------------ 图1：重建准确性 ------------------
    print("绘制图1：Reconstruction Accuracy")
    mod_super = super(type(mod), mod)  # 获取父类对象
    mod_super.plot_QC(summary_name=summary_name, use_n_obs=use_n_obs)
    plt.title("Reconstruction Accuracy")
    plt.tight_layout()
    path1 = save_dir / f"{prefix}_predict_accuracy_ref.QC.pdf"
    plt.savefig(path1, dpi=300)
    #plt.show()
    print(f"已保存图1至: {path1}")
    # ------------------ 图2：表达签名 vs 聚类均值 ------------------
    print("绘制图2：Signature vs Cluster Mean")
    inf_aver = mod.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T.copy()
    if scale_average_detection and ("detection_y_c" in mod.samples[f"post_sample_{summary_name}"]):
        inf_aver *= mod.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()
    aver = mod._compute_cluster_averages(key="labels")
    aver = aver[mod.factor_names_]
    plt.figure(figsize=(6,6))
    plt.hist2d(
        np.log10(np.asarray(aver).flatten() + 1),
        np.log10(np.asarray(inf_aver).flatten() + 1),
        bins=50,
        norm=mcolors.LogNorm()
    )
    plt.xlabel("Mean expression for every gene in every cluster")
    plt.ylabel("Estimated expression for every gene in every cluster")
    plt.title("Expression Signature vs Cluster Average")
    plt.colorbar(label='Density')
    plt.tight_layout()
    path2 = save_dir / f"{prefix}_predict_expression_ref.QC.pdf"
    plt.savefig(path2, dpi=300)
    #plt.show()
    print(f"已保存图2至: {path2}")

#####02. Fitting the reference model拟合参考模型#########
#注意：在估计参考细胞类型的signature之前，作者建议进行宽松一点的基因选择。作者更推荐于这种方式，而不是标准的高变异基因选择，因为作者的程序保留了罕见基因的标记，同时去除了大多数无信息的基因。默认参数cell_count_cutoff=5，cell_percentage_cutoff2=0.03，nonz_mean_cutoff=1.12比较合适，但是用户可以增加截断值以排除更多的基因。为了保留罕见细胞类型的标记基因，作者建议将cell_count_cutoff设置为较低的值，例如5，但是cell_percentage_cutoff2和nonz_mean_cutoff可以适当增加，可以控制在8,000至16,000个基因。
selected = filter_genes(adata_ref, cell_count_cutoff=10, cell_percentage_cutoff2=0.1, nonz_mean_cutoff=1.12)
# filter the object
adata_ref = adata_ref[:, selected].copy()
#adata_sp = adata_sp[:, selected].copy()

###Estimation of reference cell type signatures (NB regression) setup_anndata将创建实际用于训练模型的数据对象
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,batch_key='sample',labels_key='major_celltype',layer="counts",) #categorical_covariate_keys其他分类协变量（例如组织部位、处理条件等），模型会对其“校正”（视为非生物因素），不建议包含感兴趣的生物变量。
# create the regression model
mod = cell2location.models.RegressionModel(adata_ref)
####mod.train(max_epochs=250, use_gpu=True) 限速 演示用10
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002,)
#检查训练历史，可以确定模型是否需要更多的训练。该图应具有下降趋势并趋于平稳，否则需要增加max_epochs参数
# 创建一个图像和坐标轴对象
fig, ax = plt.subplots()
mod.plot_history(20, ax=ax)
plt.savefig(outdir + '02_cell2loc_anno/'+sample+'_mod_ref_plot_history.png') # 保存为PNG格式
plt.close()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
#下一步，cell2location通过将学习到的细胞丰度导出到anndata对象来总结后验分布。
# 具体来说，它计算参考数据集中每种细胞类型的表达式签名的5%、50%和95%分位数。
adata_ref = mod.export_posterior(adata_ref,sample_kwargs={"num_samples": 1000, "batch_size": 2500},)
plot_and_save_QC(mod, save_dir="./02_cell2loc_anno", prefix= sample )

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
	inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
	inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}' for i in adata_ref.uns['mod']['factor_names']]].copy()

inf_aver.columns = adata_ref.uns['mod']['factor_names']

inf_aver.head()
#如果要对同一参考数据集的多个载玻片运行反卷积，则保存细胞类型的签名可能是合理的。
inf_aver.to_csv("02_cell2loc_anno/"+ sample +"_inf_aver.csv")

####03 Cell2location: spatial mapping 细胞类型映射
# 找到单细胞和空转都有的基因，然后分别取子集。
intersect = np.intersect1d(adata_sp.var_names, inf_aver.index)
# 对adata_vis和inf_aver进行基因维度的对齐
adata_sp = adata_sp[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_sp, batch_key="sample",)
# create and train the model
mod = cell2location.models.Cell2location(adata_sp, cell_state_df=inf_aver,N_cells_per_location=20,detection_alpha=20) #作者最初选择高正则化（detection_alpha=200）作为默认值，因为作者在论文中使用的小鼠大脑和人类淋巴结数据集具有较低的技术效应，并且使用高正则化强度提高了每个位置的总估计细胞丰度与组织学定量核数之间的一致性（请参见cell2location论文的图S8F）。然而，在许多合作中，作者发现在人类组织的Visium实验中存在技术效应。这促使了将detection_alpha的新默认值设为20以及建议在用户的数据上测试这两个设置（detection_alpha=20和detection_alpha=200）。
mod.view_anndata_setup() #查看模型
#
####限速  演示小于50
# max_epochs=30000
mod.train(max_epochs=30000,batch_size=None,train_size=1)
# plot ELBO loss history during training, removing first 100 epochs from the plot
# 创建一个图像和坐标轴对象
fig, ax = plt.subplots()
mod.plot_history(1000, ax=ax) #iter_start=1000,从1000次开始展示
plt.legend(labels=['full data training']);
plt.savefig(outdir + '02_cell2loc_anno/'+sample+'_mod_sp_plot_history.png') # 保存为PNG格式
plt.close()
#再次，我们看到训练历史趋于平稳，模型收敛。我们现在可以导出每个点中估计的细胞类型丰度，并通过从后验分布中采样将其保存到空间数据中。

# adata_sp = mod.export_posterior(adata_sp)1
adata_sp = mod.export_posterior(adata_sp,sample_kwargs={'num_samples':1000,'batch_size':mod.adata.n_obs})
mod.plot_QC() #应该大致显示出噪声对角线周围的点
plt.savefig(outdir +"02_cell2loc_anno/"+ sample +'_predict_sp.QC.pdf')

adata_sp.obs[adata_sp.uns['mod']['factor_names']] = adata_sp.obsm['q05_cell_abundance_w_sf']

adata_sp.write(outdir + '02_cell2loc_anno/'+sample+'_cell2loc_anno.h5ad') # 保存为PNG格式