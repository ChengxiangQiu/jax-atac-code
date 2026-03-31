
###########################################################
### Fine-tuning the model using top 3K peaks per cell class

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

mamm_num = 32

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

adata = ad.read_h5ad(f"{work_path}/{model_id}/mamm_{mamm_num}/data_window_cluster_mouse.mamm_{mamm_num}.top3K.h5ad")

# mamm_num = 32,  N = 2,212,495 peaks

os.chdir(f"{work_path}/{model_id}/mamm_{mamm_num}")
print(os.getcwd())

strategy = tf.distribute.MirroredStrategy()

datamodule = crested.tl.data.AnnDataModule(
        adata,
        batch_size=128,
        max_stochastic_shift=3,
        always_reverse_complement=True,
    )

with strategy.scope():
    model_architecture = keras.models.load_model(
        f"{work_path}/{model_id}/mamm_{mamm_num}/window_cluster/basemodel/checkpoints/04.keras",
        compile=False,  # Choose the basemodel with best validation loss/performance metrics
    )
    optimizer = keras.optimizers.Adam(learning_rate=1e-4)
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
    trainer = crested.tl.Crested(
        data=datamodule,
        model=model_architecture,
        config=alternative_config,
        project_name="window_cluster",
        run_name="finetuned_model",
        logger="wandb",
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

adata = ad.read_h5ad(f"{work_path}/{experiment_id}/data_window_cluster_mouse.top3K.h5ad")

adata_sub = adata[:, adata.var["split"] == "test"].copy()
X_sub = adata_sub.X.toarray() if hasattr(adata_sub.X, "toarray") else adata_sub.X
obs_flat = np.log1p(X_sub.T.flatten())
var_sub = adata_sub.var.copy()

### n = 22014 test sequences (top10K)
### n = 8223 test sequences (top3K)

candidate_sequences = []
for idx, row in var_sub.iterrows():
    mid = (int(row["start"]) + int(row["end"]))/2
    candidate_sequences.append(fasta.fetch(row["chr"], mid - 1057, mid + 1057))

mamm_num = 32
print(mamm_num)
model_dir = Path(f"{work_path}/{experiment_id}/mamm_{mamm_num}/window_cluster/finetuned_model/checkpoints/")
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
df_result.to_csv(f"{work_path}/{experiment_id}/mamm_{mamm_num}/model_correlation_results.finetune_top3K.txt", sep="\t", index=False)



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

for mamm_num in [32]:
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
    model_dir = Path(f"{work_path}/{experiment_id}/mamm_{mamm_num}/window_cluster/finetuned_model/checkpoints/")
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
    df_result.to_csv(f"{work_path}/{experiment_id}/mamm_{mamm_num}/model_correlation_results.finetune_val.txt", sep="\t", index=False)


