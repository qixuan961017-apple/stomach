import pandas as pd
import scanpy as sc
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs)
from samalg import SAM
import holoviews
import bokeh
import os
import matplotlib.pyplot as plt

os.chdir('/public/home/s20223040710/rumen/05_SAMap/')
fn2="/public/home/s20223040710/rumen/05_SAMap/data/human_esophagus.h5ad"
fn1="/public/home/s20223040710/rumen/05_SAMap/data/deer_lwqw.h5ad"
fn3="/public/home/s20223040710/rumen/05_SAMap/data/SN_lwqw.h5ad"
filenames = {'de':fn1,'hu':fn2,'SN':fn3}
sm = SAMAP(filenames,f_maps = '/public/home/s20223040710/rumen/05_SAMap/01_qw/maps/')

sm.run(pairwise=True)
samap = sm.samap
keys = {'hu':'subcell','de':'subcell','SN':'subcell'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
MappingTable.to_csv('/public/home/s20223040710/rumen/05_SAMap/qw_2.txt', sep='\t', index=True)

sm.scatter()
plt.savefig('/public/home/s20223040710/rumen/05_SAMap/scatter_plot.png', dpi=300)  # 以高分辨率保存为 PNG 文件
plt.show()

#基因对查找与可视化
#gpf = GenePairFinder(sm,k1=k1,k2=k2)
gpf = GenePairFinder(sm,keys)
gene_pairs = gpf.find_all(thr=0.1)
gene_pairs.head()
gene_pairs.to_csv("Result_SAMap_gene_pairs_qw_2.csv")#将跨物种的基因对保存为 CSV 文件