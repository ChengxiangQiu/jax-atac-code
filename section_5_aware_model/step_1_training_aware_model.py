
####################################################
### Preparing the candidate window list for training

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

window_list_uniq = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/window_list_uniq.rds"))
pd = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/umap_Mus_musculus/pd.rds"))
print(sum(window_list_uniq$window_ID == pd$window_ID))
window_list_uniq$celltype = as.vector(pd$anno_L2)

window_list_uniq_no_promoter = window_list_uniq[window_list_uniq$celltype != "Promoters",]

options(scipen = 999)
write.table(window_list_uniq_no_promoter[,c("chr", "start", "end", "celltype")],
            paste0(work_path, "/14_crested/mouse_fake_track_14/candidate_region.bed"), row.names=F, col.names=F, sep="\t", quote=F)
### n = 357,282 windows

### excluding windows which are overlapped with ENCODE blacklist
### n = 354,450 windows

code_path=/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10
work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested

bedtools intersect -v \
-a "$work_path"/mouse_fake_track_14/candidate_region.bed \
-b "$code_path"/mm10-blacklist-extend-1057.bed \
> "$work_path"/mouse_fake_track_14/candidate_region_exclude_blacklist.bed

dat = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/candidate_region_exclude_blacklist.bed"))
celltype_list = names(table(dat$V4))
for(celltype in celltype_list){
    print(celltype)
    dat_sub = dat[dat$V4 == celltype,]
    write.table(dat_sub, paste0(work_path, "/14_crested/mouse_fake_track_14/tmp/", celltype, ".bed"), row.names=F, col.names=F, sep="\t", quote=F)
}

celltype_list=(Adipocyte_cells Adipocyte_cells_Cyp2e1 B_cells Brain_capillary_endothelial_cells Cardiomyocytes CNS_neurons Corticofugal_neurons Endocardial_cells Endothelium Epithelial_cells Erythroid_cells Eye Glia Glomerular_endothelial_cells Gut_epithelial_cells Hepatocytes Intermediate_neuronal_progenitors Kidney Lateral_plate_and_intermediate_mesoderm Liver_sinusoidal_endothelial_cells Lung_and_airway Lymphatic_vessel_endothelial_cells Melanocyte_cells Mesoderm Neural_crest_PNS_neurons Neuroectoderm_and_glia Olfactory_ensheathing_cells Olfactory_neurons Oligodendrocytes Skeletal_muscle_cells T_cells White_blood_cells)
for celltype in "${celltype_list[@]}"; do
    echo "${celltype}"
    bedtools sort -i "${celltype}.bed" | bedtools merge -i - > "${celltype}.merge.bed"
done

dat = NULL
for(celltype in celltype_list){
    print(celltype)
    dat_sub = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/tmp/", celltype, ".merge.bed"))
    dat_sub$celltype = celltype
    dat = rbind(dat, dat_sub)
}
### 74,621 regions after merging windows within each cell class; 35,445,000

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

adata_new = ad.AnnData(X=X_new, obs=pd.DataFrame(index=celltype_list, data={"file_path": celltype_list}), var=adata.var)

adata = adata_new
### 32 × 341164

# Choose the chromosomes for the validation and test sets
crested.pp.train_val_test_split(
    adata, strategy="chr", val_chroms=["chr8", "chr10"], test_chroms=["chr9", "chr18"]
)

print(adata.var["split"].value_counts())
adata.var.head(3)

#train    280060
#val       30595
#test      30509

crested.pp.change_regions_width(
    adata, 2114
)  # change the adata width of the regions to 2114bp

crested.pp.normalize_peaks(
    adata, top_k_percent=0.03
)  # The top_k_percent parameters can be tuned based on potential bias towards cell types. If some weights are overcompensating too much, consider increasing the top_k_percent. Default is 0.01

# Save the final preprocessing results
adata.write_h5ad(os.path.join(work_path, "mouse_fake_track_14", "data_window_cluster.h5ad"))

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
adata_sub.write_h5ad(os.path.join(work_path, "mouse_fake_track_14", "data_window_cluster_top10K.h5ad"))

adata = ad.read_h5ad(os.path.join(work_path, "mouse_fake_track_14", "data_window_cluster_top10K.h5ad"))
adata

# AnnData object with n_obs × n_vars = 32 × 238538

print(adata.var["split"].value_counts())
adata.var.head(3)

#train    195207
#test      22014
#val       2131

datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=256,  # lower this if you encounter OOM errors
    max_stochastic_shift=3,  # optional data augmentation to slightly reduce overfitting
    always_reverse_complement=True,  # default True. Will double the effective size of the training dataset.
)

# Load chrombpnet-like architecture for a dataset with 2114bp regions and 19 cell types
model_architecture = crested.tl.zoo.dilated_cnn(
    seq_len=2114, num_classes=len(list(adata.obs_names))
)

# Load the default configuration for training a peak regression model
config = crested.tl.default_configs(
    "peak_regression"
)  # or "topic_classification" for topic classification
print(config)

import keras

# Create your own configuration
# I recommend trying this for peak regression with a weighted cosine mse log loss function
optimizer = keras.optimizers.Adam(learning_rate=1e-3)
loss = crested.tl.losses.CosineMSELogLoss(max_weight=100, multiplier=1)
metrics = [
    keras.metrics.MeanAbsoluteError(),
    keras.metrics.MeanSquaredError(),
    keras.metrics.CosineSimilarity(axis=1),
    crested.tl.metrics.PearsonCorrelation(),
    crested.tl.metrics.ConcordanceCorrelationCoefficient(),
    crested.tl.metrics.PearsonCorrelationLog(),
    crested.tl.metrics.ZeroPenaltyMetric(),
]

alternative_config = crested.tl.TaskConfig(optimizer, loss, metrics)
print(alternative_config)

# setup the trainer
trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=alternative_config,
    project_name="window_cluster",  # change to your liking
    run_name="basemodel",  # change to your liking
    logger="wandb",  # or None, 'dvc', 'tensorboard'
    seed=7,  # For reproducibility
)

# train the model
trainer.fit(
    epochs=60,
    learning_rate_reduce_patience=3,
    early_stopping_patience=6,
    model_checkpointing_best_only=False
)

wandb: Run summary:
wandb:                              batch/batch_step 24470
wandb:     batch/concordance_correlation_coefficient 0.86476
wandb:                       batch/cosine_similarity 0.92731
wandb:                           batch/learning_rate 6e-05
wandb:                                    batch/loss -0.80743
wandb:                     batch/mean_absolute_error 0.84097
wandb:                      batch/mean_squared_error 2.82984
wandb:                     batch/pearson_correlation 0.89474
wandb:                 batch/pearson_correlation_log 0.80038
wandb:                     batch/zero_penalty_metric 1497.17725
wandb:     epoch/concordance_correlation_coefficient 0.86477
wandb:                       epoch/cosine_similarity 0.92729
wandb:                                   epoch/epoch 15
wandb:                           epoch/learning_rate 2e-05
wandb:                                    epoch/loss -0.80741
wandb:                     epoch/mean_absolute_error 0.84093
wandb:                      epoch/mean_squared_error 2.82958
wandb:                     epoch/pearson_correlation 0.89473
wandb:                 epoch/pearson_correlation_log 0.80035
wandb: epoch/val_concordance_correlation_coefficient 0.3563
wandb:                   epoch/val_cosine_similarity 0.83408
wandb:                                epoch/val_loss -0.41267
wandb:                 epoch/val_mean_absolute_error 1.7789
wandb:                  epoch/val_mean_squared_error 12.7852
wandb:                 epoch/val_pearson_correlation 0.5587
wandb:             epoch/val_pearson_correlation_log 0.42541
wandb:                 epoch/val_zero_penalty_metric 1769.44568
wandb:                     epoch/zero_penalty_metric 1497.28674


############################################
### Finetuning on cell type-specific regions

X_softmax = np.apply_along_axis(compute_softmax_stable, axis=0, arr=adata.X)

topk = 3000
selected_columns = set()
for row in X_softmax:
    top_idx = np.argsort(row)[-topk:]  # Top k indices
    selected_columns.update(top_idx)

selected_columns = sorted(list(selected_columns))
adata_sub = adata[:, selected_columns].copy()
adata_sub.write_h5ad(os.path.join(work_path, "mouse_fake_track_14", "data_window_cluster_top3K.h5ad"))

adata = ad.read_h5ad(os.path.join(work_path, "mouse_fake_track_14", "data_window_cluster_top3K.h5ad"))
adata
# AnnData object with n_obs × n_vars = 32 × 95139

datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=64,  # Recommended to go for a smaller batch size than in the basemodel
    max_stochastic_shift=3,
    always_reverse_complement=True,
)

# First load the pretrained model on all peaks
model_architecture = keras.models.load_model(
    f"{work_path}/mouse_fake_track_14/window_cluster/basemodel/checkpoints/15.keras",
    compile=False,  # Choose the basemodel with best validation loss/performance metrics
)

# Use the same config you used for the pretrained model. EXCEPT THE LEARNING RATE, make sure that is lower than it was on the epoch you select the model from.
optimizer = keras.optimizers.Adam(learning_rate=1e-4)  # Lower LR!
loss = crested.tl.losses.CosineMSELogLoss(max_weight=100, multiplier=1)
metrics = [
    keras.metrics.MeanAbsoluteError(),
    keras.metrics.MeanSquaredError(),
    keras.metrics.CosineSimilarity(axis=1),
    crested.tl.metrics.PearsonCorrelation(),
    crested.tl.metrics.ConcordanceCorrelationCoefficient(),
    crested.tl.metrics.PearsonCorrelationLog(),
    crested.tl.metrics.ZeroPenaltyMetric(),
]

alternative_config = crested.tl.TaskConfig(optimizer, loss, metrics)
print(alternative_config)

# setup the trainer
trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=alternative_config,
    project_name="window_cluster",  # change to your liking
    run_name="finetuned_model",  # change to your liking
    logger="wandb",  # or 'wandb', 'tensorboard'
)

trainer.fit(
    epochs=60,
    learning_rate_reduce_patience=3,
    early_stopping_patience=6,
    model_checkpointing_best_only=False
)

wandb: Run summary:
wandb:                              batch/batch_step 17210
wandb:     batch/concordance_correlation_coefficient 0.88564
wandb:                       batch/cosine_similarity 0.96575
wandb:                           batch/learning_rate 2e-05
wandb:                                    batch/loss -0.85889
wandb:                     batch/mean_absolute_error 0.93462
wandb:                      batch/mean_squared_error 3.71274
wandb:                     batch/pearson_correlation 0.91373
wandb:                 batch/pearson_correlation_log 0.81949
wandb:                     batch/zero_penalty_metric 93.38518
wandb:     epoch/concordance_correlation_coefficient 0.88562
wandb:                       epoch/cosine_similarity 0.96576
wandb:                                   epoch/epoch 6
wandb:                           epoch/learning_rate 1e-05
wandb:                                    epoch/loss -0.8589
wandb:                     epoch/mean_absolute_error 0.93462
wandb:                      epoch/mean_squared_error 3.71397
wandb:                     epoch/pearson_correlation 0.91371
wandb:                 epoch/pearson_correlation_log 0.81948
wandb: epoch/val_concordance_correlation_coefficient 0.40204
wandb:                   epoch/val_cosine_similarity 0.84401
wandb:                                epoch/val_loss -0.44335
wandb:                 epoch/val_mean_absolute_error 2.12861
wandb:                  epoch/val_mean_squared_error 18.677
wandb:                 epoch/val_pearson_correlation 0.62174
wandb:             epoch/val_pearson_correlation_log 0.6374
wandb:                 epoch/val_zero_penalty_metric 86.98798
wandb:                     epoch/zero_penalty_metric 93.38563

######################
### Evaluate the model

adata = ad.read_h5ad(os.path.join(work_path, "mouse_fake_track_14", "data_window_cluster_top3K.h5ad"))

datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=256,  # lower this if you encounter OOM errors
)

# load an existing model
evaluator = crested.tl.Crested(data=datamodule)
model_path = f"{work_path}/mouse_fake_track_14/window_cluster/finetuned_model/checkpoints/06.keras"

evaluator.load_model(
    model_path,
    compile=True,
)

# evaluate the model on the test set
evaluator.test()

2026-01-22T18:25:08.107853-0800 INFO Test concordance_correlation_coefficient: 0.3854
2026-01-22T18:25:08.108033-0800 INFO Test cosine_similarity: 0.8430
2026-01-22T18:25:08.108084-0800 INFO Test loss: 0.1457
2026-01-22T18:25:08.108124-0800 INFO Test mean_absolute_error: 2.0451
2026-01-22T18:25:08.108162-0800 INFO Test mean_squared_error: 16.9867
2026-01-22T18:25:08.108212-0800 INFO Test pearson_correlation: 0.6103
2026-01-22T18:25:08.108282-0800 INFO Test pearson_correlation_log: 0.6171
2026-01-22T18:25:08.108344-0800 INFO Test zero_penalty_metric: 325.1969


############################
### Output the access./pred.

import anndata as ad
import crested
import numpy as np
import pandas as pd
import os, sys
import matplotlib
import keras
import pysam
from pathlib import Path
from scipy.stats import pearsonr

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
experiment_id = "mouse_fake_track_14"

fasta = pysam.FastaFile(f"{work_path}/genome/Mus_musculus/Mus_musculus.fa")

adata = ad.read_h5ad(f"{work_path}/{experiment_id}/data_window_cluster_top3K.h5ad")

adata_sub = adata[:, adata.var["split"] == "test"].copy()
X_sub = adata_sub.X.toarray() if hasattr(adata_sub.X, "toarray") else adata_sub.X
obs_flat = np.log1p(X_sub.T.flatten())
var_sub = adata_sub.var.copy()

candidate_sequences = []
for idx, row in var_sub.iterrows():
    candidate_sequences.append(fasta.fetch(row["chr"], int(row["start"]), int(row["end"])))

result = []
model_dir = Path(f"{work_path}/{experiment_id}/window_cluster/finetuned_model/checkpoints/")
model_list = list(model_dir.glob("*.keras"))
model_list_sorted = sorted(model_list, key=lambda p: int(p.stem))
for model_path in model_list_sorted:
    model = keras.models.load_model(model_path, compile=False)
    predictions = crested.tl.predict(input=candidate_sequences, model=model)
    pre_flat = np.log1p(predictions.flatten())
    r, p_value = pearsonr(obs_flat, pre_flat)
    model_id = int(model_path.stem)
    print(f"{model_id} / {round(r, 2)}")
    result.append((model_id, round(r, 2)))

df_result = pd.DataFrame(result, columns=["model_id", "pearson_r"])
df_result.to_csv(f"{work_path}/{experiment_id}/model_correlation_results.txt", sep="\t", index=False)

### corr = 0.68 for model 06

######### BACKUP ####################################

#########################################################
### presenting prediction performance on test set, compared to downsamplied training data

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

dat_orig = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/model_correlation_results.txt"), header=T)
dat_orig = dat_orig %>%
  slice_max(model_id, n = 2, with_ties = FALSE) %>%
  slice(2) %>% select(-model_id)

downsample_list = c(10, 20, 30, 40, 50, 60, 70, 80, 90)

dat = NULL
for(downsample_frac in downsample_list){
    dat_i = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/downsample_", downsample_frac, "/model_correlation_results.txt"), header=T)
    dat = rbind(dat, dat_i)
}

dat_downsample = dat %>% group_by(downsample) %>%
  slice_max(model_id, n = 2, with_ties = FALSE) %>%
  slice(2) %>% select(-model_id)
dat_downsample = rbind(dat_downsample, data.frame(downsample = 100, pearson_r = 0.7))
dat_downsample$downsample = factor(dat_downsample$downsample, levels = rev(c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)))

p = ggplot(data=dat_downsample, aes(x=downsample, y=pearson_r, color=downsample)) +
    geom_point() + geom_line(aes(group = 1)) +
    labs(x="% of subset", y="pearson_r", title="") +
    scale_color_viridis(discrete=TRUE) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave("~/share/pearson_r_downsample.pdf", p, height = 4, width = 3)


exclude_chrs_num_list = c(2, 4, 6, 8, 10, 12, 14)

dat = NULL
for(exclude_chrs_num in exclude_chrs_num_list){
    dat_i = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/exclude_", exclude_chrs_num, "_chrs/model_correlation_results.txt"), header=T)
    dat = rbind(dat, dat_i)
}

dat_downsample = dat %>% group_by(exclude_chrs_num) %>%
  slice_max(model_id, n = 2, with_ties = FALSE) %>%
  slice(2) %>% select(-model_id)
dat_downsample = rbind(data.frame(exclude_chrs_num = 0, pearson_r = 0.7), dat_downsample)
dat_downsample$exclude_chrs_num = factor(dat_downsample$exclude_chrs_num, levels = c(0, 2, 4, 6, 8, 10, 12, 14))

p = ggplot(data=dat_downsample, aes(x=exclude_chrs_num, y=pearson_r, color=exclude_chrs_num)) +
    geom_point() + geom_line(aes(group = 1)) +
    labs(x="# of chrs excluded", y="pearson_r", title="") +
    scale_color_viridis(discrete=TRUE) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave("~/share/pearson_r_exclude_chrs.pdf", p, height = 4, width = 3)





