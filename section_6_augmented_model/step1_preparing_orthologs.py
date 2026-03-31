####################################################
### This new model is trained based on the candidate windows from mouse (n = 354,450 windows),
### adding orthologs from other mammals but with mouse access. as input;
### Compared to the aware model, the only difference is adding more orthologs sequences


###################################
### check the correlation between individual windows and their orthologs 
### (not just the median corr across 240 mammals, but rather all the correlations for each window)

source("~/work/scripts/utils.R")
x = readRDS('/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/window_list_uniq.rds')
y = read.table('/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/candidate_region_exclude_blacklist.bed')
colnames(y) = c("chr", "start", "end", "celltype")
y$window_id = paste0(y$chr, "_", y$start, "_", y$end)
x = x %>% left_join(y[,c("window_id", "celltype")], by = "window_id")
x = x[!is.na(x$celltype),]
x = x[,c("window_ID", "window_id", "celltype", "chr", "start", "end")]
options(scipen = 999)
write.table(x, "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/window_include.txt", row.names=F, col.names=F, sep="\t", quote=F)



### 1. Correlation between each of 354K mouse windows with their orthologs

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

res = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/result_corr.rds"))

candidate_window = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/candidate_region_exclude_blacklist.bed"))
colnames(candidate_window) = c("chr", "start", "end", "cell_class")
candidate_window$window_id = paste0(candidate_window$chr, "_", candidate_window$start, "_", candidate_window$end)

window_id = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/umap_Mus_musculus/window_list.bed"))
colnames(window_id) = c("chr", "start", "end", "window_ID")
window_id$window_id = paste0(window_id$chr, "_", window_id$start, "_", window_id$end)

candidate_window = candidate_window %>% left_join(window_id[,c("window_id", "window_ID")], by = "window_id")

res_sub = res[res$window_ID %in% as.vector(candidate_window$window_ID),]

options(scipen = 999)
write.table(candidate_window[,c("chr","start","end","window_ID")], 
    paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/candidate_window_354K.bed"), row.names=F, col.names=F, sep="\t", quote=F)


### 2. extracting the PhyloP

cd /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/conservation
wget https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Mus_musculus/phylop.bw
mv phylop.bw 241-mammalian-2020v2.bigWig

/net/shendure/vol10/projects/cxqiu/nobackup/install/bigWigAverageOverBed \
    241-mammalian-2020v2.bigWig \
    candidate_window_354K.bed \
    phyloP_scores.354K_windows.tab

pp_score = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/phyloP_scores.354K_windows.tab"))
pp_score = pp_score[,c(1,6)]
colnames(pp_score) = c("window_ID", "phyloP_scores")

df = res_sub %>% left_join(pp_score, by = "window_ID")

p = ggplot(data = df, aes(x=corr, y=phyloP_scores)) +
    geom_hex(bins = 70) +
    labs(x="Correlation coefficient", y="PhyloP score") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave("~/share/corr_with_orthologs_for_individual_354K_windows.pdf", p, width = 5, height = 5)



### 3. instead of using all the orthologs, we only retain those orthologs with corr > 0.6 from other species
### removing noisy sequences

res_save = res_sub[res_sub$corr >= 0.6,] 
# focusing on n = 354,450 in mouse and their orthologs (n = 73,800,896)
# n = 60,223,282 highly functional conserved orthologs sequences from other 240 mammals; 81.6% of all orthologs including correlation analysis
# 60,223,282 / 354,450 = 169.9, we expanded 170 times more training sequences
# Of note, the real number could be even bigger, since a mouse sequence might have multiple orthologs in one other species

write.table(res_save[,c(1,3)], paste0(work_path, "/14_crested/mouse_fake_track_15/orthologs_seq_include.all.txt"), row.names=F, col.names=F, sep="\t", quote=F)

candidate_window_x = candidate_window[candidate_window$chr %in% paste0("chr", c(1:19,'X')),]
res_save = res_save[res_save$window_ID %in% as.vector(candidate_window_x$window_ID),]
# 57,995,373 real orthologs really included, after filtering by chrs, corresponding to 341164 windows in mouse

write.table(res_save[,c(1,3)], paste0(work_path, "/14_crested/mouse_fake_track_15/orthologs_seq_include.txt"), row.names=F, col.names=F, sep="\t", quote=F)





##############################################
### creating subset of genome for each species

import sys, os
import gzip
import pysam
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"

model_id = "mouse_fake_track_15"

with open(os.path.join(work_path, "genome", "mamm_list.txt")) as f:
    mamm_list = [line.rstrip() for line in f if line.strip()]

mamm = mamm_list[int(sys.argv[1]) - 1]
print(mamm)

fasta = pysam.FastaFile(os.path.join(work_path, "genome", mamm, f"{mamm}.fa"))

window_include = set()
with open(f"{work_path}/{model_id}/orthologs_seq_include.txt") as f:
    for line in f:
        l = line.rstrip().split('\t')
        if l[1] == mamm:
            window_include.add(l[0])

with gzip.open(f"{work_path}/genome/{mamm}/candidate_windows_sub.txt.gz", "rt") as f, \
     open(f"{work_path}/{model_id}/genome/{mamm}.fa", "w") as o, \
     open(f"{work_path}/{model_id}/genome/{mamm}.chrom.sizes", "w") as o2:
    for line in f:
        chrom, start, end, window_id = line.rstrip().split("\t")
        if window_id not in window_include:
            continue
        try:
            seq = fasta.fetch(chrom, int(start) - 10, int(end) + 10)
        except Exception as e:
            print(f"Skip {chrom}:{start}-{end} ({window_id}) -> {e}")
            continue
        header = f"{window_id}_{mamm}_{chrom}_{start}_{end}"
        o.write(f">{header}\n{seq}\n")
        o2.write(f"{header}\t{len(seq)}\n")



#############################
### mouse is a little special

mamm = "Mus_musculus"
fasta = pysam.FastaFile(os.path.join(work_path, "genome", mamm, f"{mamm}.fa"))

window_include = set()
with open(f"{work_path}/{model_id}/orthologs_seq_include.txt") as f:
    for line in f:
        l = line.rstrip().split('\t')
        window_include.add(l[0])

with open(f"{work_path}/mouse_fake_track_14/cross_species/window_include.txt") as f, \
     open(f"{work_path}/{model_id}/genome/{mamm}.fa", 'w') as o, \
     open(f"{work_path}/{model_id}/genome/{mamm}.chrom.sizes", 'w') as o2:
    for line in f:
        l = line.rstrip().split('\t')
        if l[0] in window_include:
            window_id = l[0]
            chrom = l[3]
            mid = int((int(l[4]) + int(l[5]))/2)
            start = mid - 1057; end = mid + 1057
            try:
                seq = fasta.fetch(chrom, int(start) - 10, int(end) + 10)
            except Exception as e:
                print(f"Skip {chrom}:{start}-{end} ({window_id}) -> {e}")
                continue
            header = f"{window_id}_{mamm}_{chrom}_{start}_{end}"
            o.write(f">{header}\n{seq}\n")
            o2.write(f"{header}\t{len(seq)}\n")


#############################################
### preparing the pseudo-genome and adata.var

import sys, os
import random

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
model_id = "mouse_fake_track_15"

with open(os.path.join(work_path, "genome", "mamm_list.txt")) as f:
    mamm_list = [line.rstrip() for line in f if line.strip()]

with open(f"{work_path}/{model_id}/create_genome.sh", "w") as o:
    o.write("cat \\\n")
    for mamm in mamm_list:
        o.write(f"{work_path}/{model_id}/genome/{mamm}.fa \\\n")
    o.write(f"> {work_path}/{model_id}/manu_genome.fa\n")
    o.write("cat \\\n")
    for mamm in mamm_list:
        o.write(f"{work_path}/{model_id}/genome/{mamm}.chrom.sizes \\\n")
    o.write(f"> {work_path}/{model_id}/manu_genome.chrom.sizes\n")

random.seed(0)
mamm_include = ['Mus_musculus', 'Homo_sapiens']
mamm_include_x = random.sample([i for i in mamm_list if i not in mamm_include], len(mamm_list) - len(mamm_include))
mamm_include.extend(mamm_include_x)

for mamm_num in [1,2,4,8,16,32,64,128,241]:
    print(mamm_num)
    with open(f"{work_path}/{model_id}/mamm_{mamm_num}/manu_var_{mamm_num}.txt", "w") as o:
        for mamm in mamm_include[0:mamm_num]:
            with open(f"{work_path}/{model_id}/genome/{mamm}.chrom.sizes") as f:
                for line in f:
                    l = line.rstrip().split('\t')
                    window_id = l[0].split('_')[0]
                    o.write(f"{l[0]}:10-2124\t{l[0]}\t10\t2124\t{window_id}\n")






##################################################################
### Training CREsted model using cut-site based bigwig file, 
### and peaks called on profiles of individual level-2 cell types

import anndata as ad
import crested
import numpy as np
import pandas as pd
import os, sys
import matplotlib

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
model_id = "mouse_fake_track_15"

bigwigs_folder = "/net/shendure/vol2/projects/cxqiu/JAX_atac/Novaseq/call_peaks_celltype_L2/BigWig_cut_site_norm"
regions_file = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/candidate_region_exclude_blacklist.bed"

# Set the genome (this only includes regular chrs)
genome = crested.Genome(
    "/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10/mm10.fa",
    "/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10/chromosome_sizes.txt"
)
crested.register_genome(
    genome
)  # Register the genome so that it can be used by the package

print(genome.fetch("chr1", 10000000, 10000010))

adata = crested.import_bigwigs(
    bigwigs_folder=bigwigs_folder,
    regions_file=regions_file,
    target_region_width=1000,  # optionally, use a different width than the consensus regions file (500bp) for the .X values calculation
    target="count",  # or "max", "count", "logcount" --> what we will be predicting
)
adata
#AnnData object with n_obs × n_vars = 36 × 341164 (some windows were excluded because they are not on the regular chrs)

### merging cell classes to create window clusters
celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]
celltype_dict = {}
with open(f"{work_path}/mouse_fake_track_12/celltype_list_convert.txt") as f:
    for line in f:
        l = line.rstrip().split('\t')
        celltype_dict[l[0]] = celltype_dict.get(l[0], [])
        celltype_dict[l[0]].append(l[2])

X_new = []
for celltype in celltype_list:
    print(celltype)
    mask = adata.obs.index.isin(celltype_dict[celltype])
    X_sub = adata.X[mask]
    if X_sub.shape[0] > 0:
        X_new.append(np.asarray(X_sub).max(axis=0))

X_new = np.array(X_new)

chr_list = [f"chr{i}" for i in range(1, 20)]
chr_list.append("chrX")
window_include = {}
with open(f"{work_path}/mouse_fake_track_14/window_include.txt") as f:
    for line in f:
        l = line.rstrip().split('\t')
        if l[3] in chr_list:
            window_include[f"{l[3]}:{l[4]}-{l[5]}"] = l[0]

adata.var['window_id'] = [window_include[i] for i in adata.var.index]

adata_new = ad.AnnData(X=X_new, obs=pd.DataFrame(index=celltype_list, data={"file_path": celltype_list}), var=adata.var)
### 32 × 341164

# Choose the chromosomes for the validation and test sets
crested.pp.train_val_test_split(
    adata_new, strategy="chr", val_chroms=["chr8", "chr10"], test_chroms=["chr9", "chr18"]
)

### without any normalization
adata_new.write_h5ad(os.path.join(work_path, model_id, "data_window_cluster_mouse.h5ad"))


adata = ad.read_h5ad(os.path.join(work_path, model_id, "data_window_cluster_mouse.h5ad"))
crested.pp.normalize_peaks(
    adata, top_k_percent=0.03
)  # The top_k_percent parameters can be tuned based on potential bias towards cell types. If some weights are overcompensating too much, consider increasing the top_k_percent. Default is 0.01

### subset 10K top peaks in each cell type
def compute_softmax_stable(x):
    x = np.array(x, dtype=np.float64)
    exp_x = np.exp(x - np.max(x))
    return exp_x / np.sum(exp_x)

X_softmax = np.apply_along_axis(compute_softmax_stable, axis=0, arr=adata.X)

topk = 10000
selected_columns = set()
for row in X_softmax:
    top_idx = np.argsort(row)[-topk:]  # Top k indices
    selected_columns.update(top_idx)

selected_columns = sorted(list(selected_columns))
adata_sub = adata[:, selected_columns].copy()
adata_sub.write_h5ad(os.path.join(work_path, model_id, "data_window_cluster_mouse.top10K.h5ad"))


topk = 3000
selected_columns = set()
for row in X_softmax:
    top_idx = np.argsort(row)[-topk:]  # Top k indices
    selected_columns.update(top_idx)

selected_columns = sorted(list(selected_columns))
adata_sub = adata[:, selected_columns].copy()
adata_sub.write_h5ad(os.path.join(work_path, model_id, "data_window_cluster_mouse.top3K.h5ad"))




##################################################
### Adding orthologs sequences to the training set

import anndata as ad
import crested
import numpy as np
import pandas as pd
import os, sys
import matplotlib
import tensorflow as tf
print(tf.config.list_physical_devices('GPU'))

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
model_id = "mouse_fake_track_15"

# Set the genome (this only includes regular chrs)
genome = crested.Genome(
    f"{work_path}/{model_id}/manu_genome.fa",
    f"{work_path}/{model_id}/manu_genome.chrom.sizes"
)
crested.register_genome(
    genome
)  # Register the genome so that it can be used by the package
print(genome.fetch("I841531_Acinonyx_jubatus_LLWD01000002.1_233225_235339", 10, 2124))



################
### PREPARE DATA

for mamm_num in [1,2,4,8,16,32,64,128,241]:
    print(mamm_num)
    adata_orig = ad.read_h5ad(os.path.join(work_path, model_id, "data_window_cluster_mouse.h5ad"))
    adata_orig_top10K = ad.read_h5ad(os.path.join(work_path, model_id, "data_window_cluster_mouse.top10K.h5ad"))
    manu_var = pd.read_csv(f"{work_path}/{model_id}/mamm_{mamm_num}/manu_var_{mamm_num}.txt", sep="\t", header=None, index_col=0, names=["chr", "start", "end", "window_id"])
    manu_var_subset = manu_var[manu_var["window_id"].isin(adata_orig_top10K.var["window_id"])]
    adata_orig.var_names = adata_orig.var["window_id"]
    adata = adata_orig[:, manu_var_subset["window_id"]]
    (adata.var['window_id'].values == manu_var_subset['window_id'].values).all()
    adata = adata.copy()
    adata.var['chr'] = manu_var_subset['chr'].values
    adata.var['start'] = manu_var_subset['start'].values
    adata.var['end'] = manu_var_subset['end'].values
    adata.var["region"] = (
        adata.var["chr"].astype(str) + ":" +
        adata.var["start"].astype(str) + "-" +
        adata.var["end"].astype(str))
    adata.var_names = adata.var["region"]
    crested.pp.normalize_peaks(
        adata, top_k_percent=0.03)
    adata.write_h5ad(f"{work_path}/{model_id}/mamm_{mamm_num}/data_window_cluster_mouse.mamm_{mamm_num}.h5ad", compression="gzip")
    print(f"mamm_num = {mamm_num}: N = {adata.shape[1]} peaks")


# mamm_num = 1,   N = 238,538 peaks
# mamm_num = 2,   N = 431,420 peaks
# mamm_num = 4,   N = 771,468 peaks
# mamm_num = 8,   N = 1,477,022 peaks
# mamm_num = 16,  N = 2,922,771 peaks
# mamm_num = 32,  N = 5,525,475 peaks
# mamm_num = 64,  N = 10,987,299 peaks
# mamm_num = 128, N = 21,964,576 peaks
# mamm_num = 241, N = 40,772,560 peaks



mamm_num = 32
print(mamm_num)
adata_orig = ad.read_h5ad(os.path.join(work_path, model_id, "data_window_cluster_mouse.h5ad"))
adata_orig_top3K = ad.read_h5ad(os.path.join(work_path, model_id, "data_window_cluster_mouse.top3K.h5ad"))
manu_var = pd.read_csv(f"{work_path}/{model_id}/mamm_{mamm_num}/manu_var_{mamm_num}.txt", sep="\t", header=None, index_col=0, names=["chr", "start", "end", "window_id"])
manu_var_subset = manu_var[manu_var["window_id"].isin(adata_orig_top3K.var["window_id"])]
adata_orig.var_names = adata_orig.var["window_id"]
adata = adata_orig[:, manu_var_subset["window_id"]]
(adata.var['window_id'].values == manu_var_subset['window_id'].values).all()
adata = adata.copy()
adata.var['chr'] = manu_var_subset['chr'].values
adata.var['start'] = manu_var_subset['start'].values
adata.var['end'] = manu_var_subset['end'].values
adata.var["region"] = (
    adata.var["chr"].astype(str) + ":" +
    adata.var["start"].astype(str) + "-" +
    adata.var["end"].astype(str))
adata.var_names = adata.var["region"]
crested.pp.normalize_peaks(
    adata, top_k_percent=0.03)
adata.write_h5ad(f"{work_path}/{model_id}/mamm_{mamm_num}/data_window_cluster_mouse.mamm_{mamm_num}.top3K.h5ad", compression="gzip")
print(f"mamm_num = {mamm_num}: N = {adata.shape[1]} peaks")

# mamm_num = 32,  N = 2,212,495 peaks


