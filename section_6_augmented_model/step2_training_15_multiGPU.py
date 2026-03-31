
####################################################
### This new model is trained based on the candidate windows from mouse (n = 354,450 windows),
### adding orthologs from other mammals but with mouse access. as input;
### Compared to the aware model, the only difference is adding more orthologs sequences

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

for gpu in tf.config.list_physical_devices('GPU'):
    try:
        tf.config.experimental.set_memory_growth(gpu, True)
    except:
        pass

from tensorflow import keras

work_path = "/gpfs/projects/shendurelabcre/cxqiu/atac_seq/14_crested"
model_id = "mouse_fake_track_15"

mamm_num = int(sys.argv[1])

# Set the genome (this only includes regular chrs)
genome = crested.Genome(
    f"{work_path}/{model_id}/manu_genome.fa",
    f"{work_path}/{model_id}/manu_genome.chrom.sizes"
)
crested.register_genome(
    genome
)  # Register the genome so that it can be used by the package
print(genome.fetch("I841531_Acinonyx_jubatus_LLWD01000002.1_233225_235339", 10, 2124))


##################
### START TRAINING

adata = ad.read_h5ad(f"{work_path}/{model_id}/mamm_{mamm_num}/data_window_cluster_mouse.mamm_{mamm_num}.h5ad")

os.chdir(f"{work_path}/{model_id}/mamm_{mamm_num}")
print(os.getcwd())

strategy = tf.distribute.MirroredStrategy()

datamodule = crested.tl.data.AnnDataModule(
        adata,
        batch_size=256,
        max_stochastic_shift=3,
    )

with strategy.scope():
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
    model_architecture = crested.tl.zoo.dilated_cnn(
        seq_len=2114,
        num_classes=len(adata.obs_names)
    )
    trainer = crested.tl.Crested(
        data=datamodule,
        model=model_architecture,
        config=alternative_config,
        project_name="window_cluster",
        run_name="basemodel",
        logger="wandb",
        seed=7,
    )
    trainer.fit(
        epochs=20,
        learning_rate_reduce_patience=3,
        early_stopping_patience=5,
        model_checkpointing_best_only=False,
    )




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

work_path = "/gpfs/projects/shendurelabcre/cxqiu/atac_seq/14_crested"
experiment_id = "mouse_fake_track_15"

fasta = pysam.FastaFile(f"{work_path}/mm10.fa")

### this has been normalized, but only normalized within mouse
adata = ad.read_h5ad(f"{work_path}/{experiment_id}/data_window_cluster_mouse.top10K.h5ad")

adata_sub = adata[:, adata.var["split"] == "test"].copy()
X_sub = adata_sub.X.toarray() if hasattr(adata_sub.X, "toarray") else adata_sub.X
obs_flat = np.log1p(X_sub.T.flatten())
var_sub = adata_sub.var.copy()

### n = 22014 test sequences
candidate_sequences = []
for idx, row in var_sub.iterrows():
    mid = (int(row["start"]) + int(row["end"]))/2
    candidate_sequences.append(fasta.fetch(row["chr"], mid - 1057, mid + 1057))

for mamm_num in [1,2,4,8,16,32,64,128,241]:
    print(mamm_num)
    model_dir = Path(f"{work_path}/{experiment_id}/mamm_{mamm_num}/window_cluster/basemodel/checkpoints/")
    model_list = list(model_dir.glob("*.keras"))
    model_list_sorted = sorted(model_list, key=lambda p: int(p.stem))
    result = []
    for model_path in model_list_sorted:
        model = keras.models.load_model(model_path, compile=False)
        predictions = crested.tl.predict(input=candidate_sequences, model=model)
        pre_flat = np.log1p(predictions.flatten())
        r, p_value = pearsonr(obs_flat, pre_flat)
        model_id = int(model_path.stem)
        print(f"{model_id} / {round(r, 2)}")
        result.append((model_id, round(r, 2)))
    df_result = pd.DataFrame(result, columns=["model_id", "pearson_r"])
    df_result.to_csv(f"{work_path}/{experiment_id}/mamm_{mamm_num}/model_correlation_results.txt", sep="\t", index=False)



###################################
### Plot prediction on held-out set

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

model_id = "mouse_fake_track_15"

mamm_num_list = c(1,2,4,8,16,32,64,128,241)

dat = NULL

for(mamm_num in mamm_num_list){
    dat_i = read.table(paste0(work_path, '/14_crested/', model_id, '/mamm_', mamm_num, '/model_correlation_results.txt'), header=T)
    colnames(dat_i) = c("epoch", "pearson_r")
    dat_i$mamm_num = mamm_num
    dat = rbind(dat, dat_i)
}

dat$num_species = factor(dat$mamm_num)

p = ggplot(dat, aes(x = epoch, y = pearson_r, color = num_species, group = num_species)) +
    geom_point() +
    geom_smooth(method = 'loess', formula = 'y ~ x', se = FALSE) +
    theme_classic(base_size = 10)
ggsave("~/share/mamm_num_pearson_r.pdf", p, height = 4, width = 5)






#######################################
### Output the access./pred. on val set

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

work_path = "/gpfs/projects/shendurelabcre/cxqiu/atac_seq/14_crested"
experiment_id = "mouse_fake_track_15"

fasta = pysam.FastaFile(f"{work_path}/{experiment_id}/manu_genome.fa")

for mamm_num in [1,2,4,8,16,32,64,128,241]:
    print(mamm_num)
    adata = ad.read_h5ad(f"{work_path}/{experiment_id}/mamm_{mamm_num}/data_window_cluster_mouse.mamm_{mamm_num}.h5ad")
    adata_sub = adata[:, adata.var["split"] == "val"].copy()
    X_sub = adata_sub.X.toarray() if hasattr(adata_sub.X, "toarray") else adata_sub.X
    obs_flat = np.log1p(X_sub.T.flatten())
    var_sub = adata_sub.var.copy()
    candidate_sequences = []
    for idx, row in var_sub.iterrows():
        mid = (int(row["start"]) + int(row["end"]))/2
        candidate_sequences.append(fasta.fetch(row["chr"], mid - 1057, mid + 1057))
    model_dir = Path(f"{work_path}/{experiment_id}/mamm_{mamm_num}/window_cluster/basemodel/checkpoints/")
    model_list = list(model_dir.glob("*.keras"))
    model_list_sorted = sorted(model_list, key=lambda p: int(p.stem))
    result = []
    for model_path in model_list_sorted:
        model = keras.models.load_model(model_path, compile=False)
        predictions = crested.tl.predict(input=candidate_sequences, model=model)
        pre_flat = np.log1p(predictions.flatten())
        r, p_value = pearsonr(obs_flat, pre_flat)
        model_id = int(model_path.stem)
        print(f"{model_id} / {round(r, 2)}")
        result.append((model_id, round(r, 2)))
    df_result = pd.DataFrame(result, columns=["model_id", "pearson_r"])
    df_result.to_csv(f"{work_path}/{experiment_id}/mamm_{mamm_num}/model_correlation_results.val.txt", sep="\t", index=False)



###################################
### Plot prediction on val set

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

model_id = "mouse_fake_track_15"

mamm_num_list = c(1,2,4,8,16,32,64,128,241)

dat = NULL

for(mamm_num in mamm_num_list){
    dat_i = read.table(paste0(work_path, '/14_crested/', model_id, '/mamm_', mamm_num, '/model_correlation_results.val.txt'), header=T)
    colnames(dat_i) = c("epoch", "pearson_r")
    dat_i$mamm_num = mamm_num
    dat = rbind(dat, dat_i)
}

dat$num_species = factor(dat$mamm_num)

p = ggplot(dat, aes(x = epoch, y = pearson_r, color = num_species, group = num_species)) +
    geom_point() +
    geom_smooth(method = 'loess', formula = 'y ~ x', se = FALSE) +
    theme_classic(base_size = 10)
ggsave("~/share/mamm_num_pearson_r_val.pdf", p, height = 4, width = 5)





