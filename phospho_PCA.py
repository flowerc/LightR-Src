'''The following script performs principal component analysis (PCA) on supplementary phosphoproteomic data as reported in:
Shaaya M, et al. Light-regulated allosteric switch enables temporal and subcellular control of enzyme activity. Elife. 2020;9:e60647.

Run this script in the same directory as the Excel (xlsx)-formatted phosphoproteomic data.
'''

from os import getcwd
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# import phosphoproteomic data
# (assumes abundances are already normalized to basal abundances in LightR-Src cells for each replicate)
path_data = getcwd() + '\\82058_0_supp_1617734_qd653y.xlsx'
df_phos = pd.read_excel(path_data).set_index('Annotated Sequence')

# average the three biological replicates, drop all other columns
conditions = ['10 s', '30 s', '1 min', '5 min', '60 min', '10 s ctrl', '30 s ctrl', '1 min ctrl', '5 min ctrl', '60 min ctrl']
for c in conditions: df_phos[c] = df_phos[['Sum of ' + c, 'Sum of ' + c + '2', 'Sum of ' + c + '3']].mean(axis=1)
df_phos_avg = df_phos[conditions]

# log2-transform and recenter data around 0 prior to PCA
df_phos_avg = df_phos_avg.apply(np.log2)
df_phos_avg = df_phos_avg.sub(df_phos_avg.mean(axis=1), axis=0)
df_phos_avg = df_phos_avg.T

# perform PCA
labels = pd.Series(df_phos_avg.index)
pca = PCA(n_components=2)
df_pca = pd.DataFrame(pca.fit_transform(df_phos_avg))
df_pca = pd.concat([df_pca, labels], axis=1)
var = pca.explained_variance_ratio_
df_pca.columns = ['PC1','PC2','Label']
PC1 = pca.components_[0]
PC2 = pca.components_[1]

# plot PCA-transformed data with labels
for i in range(df_pca.shape[0]):
    pc1 = df_pca.iloc[i]['PC1']
    pc2 = df_pca.iloc[i]['PC2']
    plt.scatter(pc1, pc2, c='k', s=40)
    label = df_pca.iloc[i]['Label']
    # arrange each label for readability
    if 'ctrl' in label:
        label = label.split(' ctrl')[0]
        if label == '10 s': plt.annotate(label, (-14,0.6), size=12)
        else: plt.annotate(' ' + label, (pc1, pc2), size=12)  
    else:
        if label == '10 s': plt.annotate(label, (-5,-3), size=12)
        elif label == '30 s': plt.annotate(label, (-1,-4), size=12)
        elif label == '1 min': plt.annotate(label, (2.7,-2.9), size=12)
        elif label == '5 min': plt.annotate(label, (15.5,1.2), size=12)
        elif label == '60 min': plt.annotate(label, (26,0.5), size=12)

# draw arrows between points to depict temporal trajectories
arrows = {
    '10 s': ('30 s', 0),
    '30 s': ('1 min', 0),
    '1 min': ('5 min', 0),
    '5 min': ('60 min', 0),
    '10 s ctrl': ('30 s ctrl', 1),
    '30 s ctrl': ('1 min ctrl', 1),
    '1 min ctrl': ('5 min ctrl', 1),
    '5 min ctrl': ('60 min ctrl', 1)
}
arrow_colors = ['red', 'blue']
arrow_labels = ['LightR-Src HeLa', 'HeLa']
legend_dict = {}
for sender in arrows.keys():
    receivers = arrows[sender]
    if type(receivers) is list:
        for receiver, cindex in receivers:
            x1 = df_pca[df_pca['Label'] == sender]['PC1'].values[0]
            y1 = df_pca[df_pca['Label'] == sender]['PC2'].values[0]
            x2 = df_pca[df_pca['Label'] == receiver]['PC1'].values[0]
            y2 = df_pca[df_pca['Label'] == receiver]['PC2'].values[0]
            a = plt.arrow(x1, y1, x2-x1, y2-y1, color=arrow_colors[cindex], width=0.001, head_width=0.5, length_includes_head=True, overhang=1)
            if len(arrow_labels) > 0: legend_dict[arrow_labels[cindex]] = a
    else:
        receiver, cindex = receivers
        x1 = df_pca[df_pca['Label'] == sender]['PC1'].values[0]
        y1 = df_pca[df_pca['Label'] == sender]['PC2'].values[0]
        x2 = df_pca[df_pca['Label'] == receiver]['PC1'].values[0]
        y2 = df_pca[df_pca['Label'] == receiver]['PC2'].values[0]
        a = plt.arrow(x1, y1, x2-x1, y2-y1, color=arrow_colors[cindex], width=0.001, head_width=0.5, length_includes_head=True, overhang=1)
        if len(arrow_labels) > 0: legend_dict[arrow_labels[cindex]] = a
plt.xlabel('PC1 (' + str(round(var[0]*100,1)) + '%)', fontsize=16)
plt.ylabel('PC2 (' + str(round(var[1]*100,1)) + '%)', fontsize=16)
if len(arrow_labels) > 0: plt.legend(legend_dict.values(), legend_dict.keys(), fontsize=12)
plt.show()

# print the gene IDs of all phosphopeptides with loading score greater than 0.05 in principal component 1
gene_dict = df_phos['Gene Name'].to_dict()
for PC1_loading, peptide in zip(PC1, df_phos_avg.columns):
    if PC1_loading > 0.05: print(gene_dict[peptide])
