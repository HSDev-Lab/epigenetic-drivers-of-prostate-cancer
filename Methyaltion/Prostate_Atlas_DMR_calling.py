from os import listdir
from os.path import isfile, join


mypath = "/Users/tobymilne-clark/Documents/Data/Atlas/Prostate_Atlas/Beta"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

os.chdir(mypath)

beta_files = load_methylation_frequencies(onlyfiles,'beta')

op.splitext(onlyfiles[-1])

betas = beta_files[:, :, 0] / beta_files[:, :, 1]


beta_cov = covariance_blocks(betas)
blks = extract_blocks(beta_cov, DMR_length_threshold=4, cov_threshold=0.05)


DMRS = []
for j in range(len(blks)):
    DMRS.append(np.mean(betas[range(blks[j,0],blks[j,1]+1),],axis=0))

DMRS_np = np.array(DMRS)

laplacian_score_DMRs = compute_laplacian_score(DMRS_np.transpose())

meta = np.column_stack((['endothelial_P8100.dedup.sorted.beta', 'basal_P8100.dedup.sorted.beta', 'Luminal_P8116.dedup.sorted.beta', 'SmoothMuscle_P8116.dedup.sorted.beta', 'Immune_P8113.dedup.sorted.beta', 'SmoothMuscle_P8113.dedup.sorted.beta', 'Immune_P8116.dedup.sorted.beta', 
 'Basal_P8107.dedup.sorted.beta', 'Fibroblast_P8113.dedup.sorted.beta', 'Immune_P8107.dedup.sorted.beta', 'luminal_P8100.dedup.sorted.beta', 'Basal_P8113.dedup.sorted.beta', 'immune_P8100.dedup.sorted.beta', 'smooth-muscle_P8100.dedup.sorted.beta', 
 'Endo_P8113.dedup.sorted.beta', 'Luminal_P8107.dedup.sorted.beta', 'Endo_P8116.dedup.sorted.beta', 'fibroblast_P8100.dedup.sorted.beta'],
[4,1,2,6,3,6,3,1,5,3,2,1,3,6,4,2,4,5],
["Endothelial", "Basal", "Luminal", "SmoothMuscle", "Immune", "SmoothMuscle", "Immune", "Basal", "Fibroblast", "Immune", "Luminal", "Basal", "Immune", "SmoothMuscle", "Endothelial", "Luminal", "Endothelial", "Fibroblast"],
["Stromal", "Epithelial", "Epithelial", "Stromal", "Immune", "Stromal", "Immune", "Epithelial", "Stromal", "Immune", "Epithelial", "Epithelial", "Immune", "Stromal", "Stromal", "Epithelial", "Stromal", "Stromal"],
["Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate","Prostate"],))

def fdr(p_vals):
    ranked_p_values = stats.rankdata(p_vals)
    fdr = p_vals * len(p_vals) / ranked_p_values
    fdr[fdr > 1] = 1
    return fdr

def within_mean(data,meta,col):
    within = np.empty((0))
    for i in np.unique(meta[:,col]):
#        min_idx = np.min(np.argwhere(meta[:,col]==i))
#        max_idx = np.max(np.argwhere(meta[:,col]==i))
        within = np.append(within,np.mean(data[:,np.argwhere(meta[:,col]==i)]))
    return within


def find_seperators(df,region,meta,group = False, groups_meta = None ,col = 2,thres=0.5):
    within = within_mean(df[region,],meta,1)
    if group == True:
        avg = within
        meta = groups_meta
    else:
        avg = np.mean(df[region,],axis=0)
    index = np.argsort(avg)
    values = np.atleast_1d(index[0])
    cpg_start = region[0]
    cpg_end = region[-1]
    for i in range(len(index)-1):
        diff = np.abs(avg[index[i]]-avg[index[i+1]])
        if diff < thres:
            values = np.append(values, index[i+1])
        if diff > thres:
            break
    if len(values) == len(index):
        values = "Uncertain"
        remaining = "Uncertain"
        ct = "Uncertain"
        cell = "Uncertain"
        origin = "Uncertain"
        direction = "Uncertain"
        tp = "Uncertain"
        diff = None
        pval = 1 # maybe do first group identified vs Rest rather than assigning 1
    else:
        remaining = np.delete(index,range(len(values)))
        origin = np.apply_along_axis(','.join, 0,np.unique(meta[values,4]))
        cell = np.apply_along_axis(','.join, 0,np.unique(meta[values,3]))
        if len(remaining) > len(values):
            direction = "Hypomethylation"
        else :
            tmp = values
            values = remaining
            remaining = tmp
            cell = np.apply_along_axis(','.join, 0,np.unique(meta[values,3]))
            origin = np.apply_along_axis(','.join, 0,np.unique(meta[values,4]))
            direction = "Hypermethylation"
        if len(np.unique(meta[values,col])) == 1 :
            ct = np.apply_along_axis(','.join, 0, np.unique(meta[values,col]))
            tp = 'type1'
        else:
            ct = "Mixed"
            if np.atleast_1d(len(np.unique(meta[values,3]))) == 1 and np.atleast_1d(len(np.unique(meta[values,4]))) == 1:
                tp = 'type2'
            elif np.atleast_1d(len(np.unique(meta[values,4]))) == 1:
                tp = 'type3'
            elif np.atleast_1d(len(np.unique(meta[values,3]))) == 1:
                tp = 'type4'
            else:
                tp = 'Uncertain'              
        pval = stats.ttest_ind(avg[values],avg[remaining], nan_policy='omit')[1] 
    return {'CpG_start' : cpg_start ,
            'CpG_end' : cpg_end,
            'sample_index' : values,
            'remaining_index' : remaining, 
            'celltype_1' : ct,
            'celltype_2' : cell,
            'tissue' : origin,
            'type' : tp,
            'direction' : direction, 
            'delta_minmax' : diff,
            'pvalue' : pval,
            'averages' : within}



DMRS_anno = []
for j in range(len(blks)):
    DMRS_anno.append(find_seperators(betas,range(blks[j,0],blks[j,1]+1),meta = meta,group = False ,col = 2 ,thres=0.4))

DMRS_anno_df = pd.DataFrame(DMRS_anno)
DMRS_anno_df.celltype_1 = DMRS_anno_df.celltype_1.astype(str)
DMRS_anno_df.celltype_2 = DMRS_anno_df.celltype_2.astype(str)
DMRS_anno_df.tissue = DMRS_anno_df.tissue.astype(str)
DMRS_anno_df.delta_minmax = DMRS_anno_df.delta_minmax.astype(float)
DMRS_anno_df['fdr'] = fdr(DMRS_anno_df['pvalue'].to_numpy(dtype='float'))
DMRS_anno_df['laplacian_score'] = laplacian_score_DMRs

"/Users/tobymilne-clark/Documents/Data/Atlas/hg38/CpG.bed.gz"

location = pd.read_csv('/Users/tobymilne-clark/Documents/Data/Atlas/hg38/CpG.bed.gz', compression='gzip', header=None, sep='\t', quotechar='"')

DMR_location = []
for i in range(len(DMRS_anno_df)):
   DMR_location.append((location.iloc[DMRS_anno_df.loc[i,'CpG_start'],0],
    location.iloc[DMRS_anno_df.loc[i,'CpG_start'],1],
    location.iloc[DMRS_anno_df.loc[i,'CpG_end'],1]))
   
DMR_location_df = pd.DataFrame(DMR_location)

new_dataframe = pd.concat([DMR_location_df, DMRS_anno_df], axis=1)

new_dataframe.to_csv('DMRs_with_coords_Prostate_Atlas.csv', index=False)



top25 = new_dataframe.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & type == "type1" & direction == "Hypomethylation"').groupby('celltype_1',group_keys=False).apply(pd.DataFrame.nlargest, n=25, columns='delta_minmax')
top25_vals = pd.DataFrame(np.concatenate(top25['averages'].to_numpy()).reshape(-1,6), columns = ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"])
top25_ordered = top25_vals.reindex(sorted(top25_vals.columns), axis=1)
sns.heatmap(top25_ordered,cmap=sns.color_palette("viridis", as_cmap=True))
plt.show()

top25 = new_dataframe.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & type == "type1" & direction == "Hypermethylation"').groupby('celltype_1',group_keys=False).apply(pd.DataFrame.nlargest, n=25, columns='delta_minmax')
top25_vals = pd.DataFrame(np.concatenate(top25['averages'].to_numpy()).reshape(-1,6), columns = ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"])
top25_ordered = top25_vals.reindex(sorted(top25_vals.columns), axis=1)
sns.heatmap(top25_ordered,cmap=sns.color_palette("viridis", as_cmap=True))
plt.show()



top25 = new_dataframe.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & type == "type2" & direction == "Hypomethylation"').groupby('celltype_1',group_keys=False).apply(pd.DataFrame.nlargest, n=25, columns='delta_minmax')
top25_vals = pd.DataFrame(np.concatenate(top25['averages'].to_numpy()).reshape(-1,6), columns = ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"])
top25_ordered = top25_vals.reindex(sorted(top25_vals.columns), axis=1)
sns.heatmap(top25_ordered,cmap=sns.color_palette("viridis", as_cmap=True))
plt.show()


groups_meta = np.column_stack(([1,2,3,4,5,6],[1,2,3,4,5,6],
                        ["Basal","Luminal","Immune","Endothelial","Fibroblast","SmoothMuscle"],
                        ["Epithelial","Epithelial","Immune","Stromal","Stromal","Stromal"],
                        ["Prostate","Prostate","Prostate","Prostate","Prostate","Prostate"]))


DMRS_anno_gr = []
for j in range(len(blks)):
    DMRS_anno_gr.append(find_seperators(betas,range(blks[j,0],blks[j,1]+1),meta = meta,group = True ,groups_meta = groups_meta,col = 2 ,thres=0.4))

DMRS_anno_gr_df = pd.DataFrame(DMRS_anno_gr)
DMRS_anno_gr_df.celltype_1 = DMRS_anno_gr_df.celltype_1.astype(str)
DMRS_anno_gr_df.celltype_2 = DMRS_anno_gr_df.celltype_2.astype(str)
DMRS_anno_gr_df.tissue = DMRS_anno_gr_df.tissue.astype(str)
DMRS_anno_gr_df.delta_minmax = DMRS_anno_gr_df.delta_minmax.astype(float)
DMRS_anno_gr_df['fdr'] = fdr(DMRS_anno_gr_df['pvalue'].to_numpy(dtype='float'))
DMRS_anno_gr_df['laplacian_score'] = laplacian_score_DMRs

"/Users/tobymilne-clark/Documents/Data/Atlas/hg38/CpG.bed.gz"

location = pd.read_csv('/Users/tobymilne-clark/Documents/Data/Atlas/hg38/CpG.bed.gz', compression='gzip', header=None, sep='\t', quotechar='"')

DMR_location = []
for i in range(len(DMRS_anno_df)):
   DMR_location.append((location.iloc[DMRS_anno_df.loc[i,'CpG_start'],0],
    location.iloc[DMRS_anno_df.loc[i,'CpG_start'],1],
    location.iloc[DMRS_anno_df.loc[i,'CpG_end'],1]))
   
DMR_location_df = pd.DataFrame(DMR_location)

new_dataframe_gr = pd.concat([DMR_location_df, DMRS_anno_gr_df], axis=1)

new_dataframe_gr.to_csv('DMRs_grouped_mean_with_coords_Prostate_Atlas.csv', index=False)

top25 = new_dataframe_gr.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & type == "type1" & direction == "Hypomethylation"').groupby('celltype_1',group_keys=False).apply(pd.DataFrame.nlargest, n=25, columns='delta_minmax')
top25_vals = pd.DataFrame(np.concatenate(top25['averages'].to_numpy()).reshape(-1,6), columns = ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"])
top25_ordered = top25_vals.reindex(sorted(top25_vals.columns), axis=1)
sns.heatmap(top25_ordered,cmap=sns.color_palette("viridis", as_cmap=True))
plt.show()

top25 = new_dataframe_gr.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & type == "type1" & direction == "Hypermethylation"').groupby('celltype_1',group_keys=False).apply(pd.DataFrame.nlargest, n=25, columns='delta_minmax')
top25_vals = pd.DataFrame(np.concatenate(top25['averages'].to_numpy()).reshape(-1,6), columns = ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"])
top25_ordered = top25_vals.reindex(sorted(top25_vals.columns), axis=1)
sns.heatmap(top25_ordered,cmap=sns.color_palette("viridis", as_cmap=True))
plt.show()

top25 = new_dataframe_gr.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & type == "type2"').groupby('celltype_1',group_keys=False).apply(pd.DataFrame.nlargest, n=25, columns='delta_minmax')
top25_vals = pd.DataFrame(np.concatenate(top25['averages'].to_numpy()).reshape(-1,6), columns = ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"])
top25_ordered = top25_vals.reindex(sorted(top25_vals.columns), axis=1)
sns.heatmap(top25_ordered,cmap=sns.color_palette("viridis", as_cmap=True))
plt.show()



from matplotlib.patches import Patch

data = pd.read_csv('/Users/tobymilne-clark/Documents/Data/Atlas/Prostate_Atlas/DMRs_grouped_mean_with_coords_Prostate_Atlas_matched.csv', delimiter=',')

sorted = data.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & direction == "Hypomethylation"').groupby('celltypes_expanded',group_keys=False).apply(pd.DataFrame.nlargest, n=200, columns='delta_minmax')
sorted_vals = sorted.iloc(axis=1)[49:55]
sorted_vals.columns= ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"]

lut = dict(zip(sorted['celltypes_expanded'].unique(), sns.hls_palette(len(sorted['celltypes_expanded'].unique()), l=0.5, s=0.8)))
row_colors = sorted['celltypes_expanded'].map(lut)
plot = sns.clustermap(sorted_vals,cmap=sns.color_palette("viridis", as_cmap=True),col_cluster=False,row_cluster=False,row_colors=row_colors)

handles = [Patch(facecolor=lut[name]) for name in lut]

plot = sns.clustermap(sorted_vals,cmap=sns.color_palette("viridis", as_cmap=True),col_cluster=False,row_cluster=False,row_colors=row_colors)
plot.fig.subplots_adjust(right=0.8)
plot.ax_cbar.set_position((0.9, .3, .03, .3))
plt.legend(handles, lut, title='celltype',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='best')
plt.show()



sorted = data.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & direction == "Hypermethylation"').groupby('celltypes_expanded',group_keys=False).apply(pd.DataFrame.nlargest, n=200, columns='delta_minmax')
sorted_vals = sorted.iloc(axis=1)[49:55]
sorted_vals.columns= ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"]

lut = dict(zip(sorted['celltypes_expanded'].unique(), sns.hls_palette(len(sorted['celltypes_expanded'].unique()), l=0.5, s=0.8)))
row_colors = sorted['celltypes_expanded'].map(lut)
plot = sns.clustermap(sorted_vals,cmap=sns.color_palette("viridis", as_cmap=True),col_cluster=False,row_cluster=False,row_colors=row_colors)

handles = [Patch(facecolor=lut[name]) for name in lut]

plot = sns.clustermap(sorted_vals,cmap=sns.color_palette("viridis", as_cmap=True),col_cluster=False,row_cluster=False,row_colors=row_colors)
plot.fig.subplots_adjust(right=0.8)
plot.ax_cbar.set_position((0.9, .3, .03, .3))
plt.legend(handles, lut, title='celltype',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='best')
plt.show()


sorted = data.query('fdr < 0.05 & abs(delta_minmax) > 0.4 & type == "type1"  & direction == "Hypomethylation"').groupby('celltypes_expanded',group_keys=False).apply(pd.DataFrame.nlargest, n=200, columns='delta_minmax')
sorted_vals = sorted.iloc(axis=1)[49:55]
sorted_vals.columns= ["Basal","Luminal","Immune","Endothelial","Fibroblast","Smooth Muscle"]

lut = dict(zip(sorted['celltypes_expanded'].unique(), sns.hls_palette(len(sorted['celltypes_expanded'].unique()), l=0.5, s=0.8)))
row_colors = sorted['celltypes_expanded'].map(lut)
plot = sns.clustermap(sorted_vals,cmap=sns.color_palette("viridis", as_cmap=True),col_cluster=False,row_cluster=False,row_colors=row_colors)

handles = [Patch(facecolor=lut[name]) for name in lut]

plot = sns.clustermap(sorted_vals,cmap=sns.color_palette("viridis", as_cmap=True),col_cluster=False,row_cluster=False,row_colors=row_colors)
plot.fig.subplots_adjust(right=0.8)
plot.ax_cbar.set_position((0.9, .3, .03, .3))
plt.legend(handles, lut, title='celltype',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='best')
plt.show()
