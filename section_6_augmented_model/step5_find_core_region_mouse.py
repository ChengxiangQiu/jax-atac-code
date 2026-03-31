

###############################################
### convert predicted score to Phred like value

import gzip
import os, sys
import numpy as np

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
model_id = "mouse_fake_track_15"
mamm = "Mus_musculus"

dat = np.loadtxt(f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/dat.txt.gz")

def phred_from_column(values):
    sorted_vals = np.sort(values)
    N = len(sorted_vals)
    ranks = np.searchsorted(sorted_vals, values, side="right")
    q = (ranks - 0.5) / N
    q = np.minimum(q, 1 - 1e-12)
    phred = -10 * np.log10(1 - q)
    return phred

phred_dat = np.zeros(dat.shape, dtype=np.float64)

for j in range(dat.shape[1]):
    print(j)
    phred_dat[:, j] = phred_from_column(dat[:, j])

file_name = f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/dat.phred_score.txt.gz"
with gzip.open(file_name, 'wt') as f:
    np.savetxt(f, phred_dat, fmt="%.3f", delimiter="\t")


>>> (in R)
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
mamm = "Mus_musculus"
celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")
dat_phred = read.table(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.phred_score.txt.gz"))
colnames(dat_phred) = celltype_list
saveRDS(dat_phred, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.phred_score.txt.rds"))


dat = read.table(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.txt.gz"))
dat_loc = read.table(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat_loc.txt.gz"))

colnames(dat) = celltype_list
colnames(dat_loc) = c("chr", "start")
dat_loc$end = dat_loc$start + 100

saveRDS(dat, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.txt.rds"))
saveRDS(dat_loc, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat_loc.txt.rds"))



##################################
### step-1: call peaks on model_14

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")
library(GenomicRanges)

model_id = "mouse_fake_track_15"

mamm = "Mus_musculus"

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

dat = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.txt.rds"))
dat_loc = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat_loc.txt.rds"))

dat_phred = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.phred_score.txt.rds"))

dat = as.matrix(dat)
dat_phred = as.matrix(dat_phred)
a = 10000000/100

window_list = NULL
for(i in 1:ncol(dat)){
    celltype = celltype_list[i]
    cutoff = quantile(dat[,i], 1-a/nrow(dat))
    print(paste0(celltype, ": ", cutoff))
    keep = dat[,i] > cutoff
    dat_loc_sub = dat_loc[keep,]
    dat_loc_sub$score = dat[keep,i]
    dat_loc_sub$phred_score = dat_phred[keep,i]
    dat_loc_sub$celltype = celltype
    window_list = rbind(window_list, dat_loc_sub)
}
saveRDS(window_list, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/window_list_10Mb.rds"))

[1] "Adipocyte_cells: 5.83"
[1] "Adipocyte_cells_Cyp2e1: 5.398"
[1] "B_cells: 2.844"
[1] "Brain_capillary_endothelial_cells: 4.379"
[1] "CNS_neurons: 9.18300354827941"
[1] "Cardiomyocytes: 5.058"
[1] "Corticofugal_neurons: 6.774"
[1] "Endocardial_cells: 4.977"
[1] "Endothelium: 4.821"
[1] "Epithelial_cells: 7.1950035482794"
[1] "Erythroid_cells: 10.905"
[1] "Eye: 4.392"
[1] "Glia: 5.824"
[1] "Glomerular_endothelial_cells: 3.988"
[1] "Gut_epithelial_cells: 5.571"
[1] "Hepatocytes: 6.13"
[1] "Intermediate_neuronal_progenitors: 8.359"
[1] "Kidney: 5.962"
[1] "Lateral_plate_and_intermediate_mesoderm: 5.723"
[1] "Liver_sinusoidal_endothelial_cells: 4.678"
[1] "Lung_and_airway: 5.185"
[1] "Lymphatic_vessel_endothelial_cells: 4.061"
[1] "Melanocyte_cells: 3.498"
[1] "Mesoderm: 6.425"
[1] "Neural_crest_PNS_neurons: 8.784"
[1] "Neuroectoderm_and_glia: 5.364"
[1] "Olfactory_ensheathing_cells: 4.435"
[1] "Olfactory_neurons: 8.29800354827941"
[1] "Oligodendrocytes: 5.872"
[1] "Skeletal_muscle_cells: 6.035"
[1] "T_cells: 5.019"
[1] "White_blood_cells: 8.423"

x = window_list %>% group_by(celltype) %>% slice_min(order_by = phred_score, n = 1, with_ties = F)
### Q = 24.5
### q = 0.99645, top 0.355%

dat_all = NULL
for(i in 1:length(celltype_list)){
    celltype = celltype_list[i]
    print(paste0(i, "/", celltype))
    dat_loc_sub = window_list[window_list$celltype == celltype, 1:4]

    gr = GRanges(
      seqnames = dat_loc_sub$chr,
      ranges   = IRanges(start = dat_loc_sub$start + 1, end = dat_loc_sub$end),
      score    = dat_loc_sub$score,
    )

    gr_merged <- reduce(sort(gr), with.revmap = TRUE)

    mcols(gr_merged)$mean_score <- sapply(
      mcols(gr_merged)$revmap,
      function(i) mean(mcols(gr)$score[i], na.rm = TRUE)
    )

    df_merged <- data.frame(
      chr        = as.character(seqnames(gr_merged)),
      start      = start(gr_merged) - 1,
      end        = end(gr_merged),
      mean_score = round(mcols(gr_merged)$mean_score, 3),
      celltype   = celltype
    )
    dat_all = rbind(dat_all, df_merged)

    write.table(df_merged, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/", celltype, ".bed"), row.names=F, col.names=F, sep="\t", quote=F)
}

write.table(dat_all, paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/window_list_10Mb.merge.bed"), row.names=F, col.names=F, sep="\t", quote=F)

### n = 341,676

### this is after merging 100 bp windows, but before detecting core regions by left-right trimming







##############################################################################
### step-1: the first round of identifying core region

import random
import sys, os
import numpy as np
import gzip
import pysam
from collections import Counter

import anndata as ad
import crested
import keras

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"

celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]

celltype = celltype_list[int(sys.argv[1]) - 1]
print(celltype)

########################################
# Constants
########################################

BASES = np.array(['A','C','G','T'])
BASE_TO_INT = {b:i for i,b in enumerate(BASES)}

mamm = "Mus_musculus"
model_id = "mouse_fake_track_15"
step_size = 20
seq_length = 2114
n_reps = 10

########################################
# FASTA + model
########################################

fasta = pysam.FastaFile(f"{work_path}/genome/{mamm}/{mamm}.fa")

model_path = f"{work_path}/{model_id}/mamm_32/window_cluster/finetuned_model_1e5/checkpoints/02.keras"
model = keras.models.load_model(model_path, compile=False)

########################################
# Chrom sizes
########################################

chr_size = {}
with open(f"{work_path}/genome/{mamm}/{mamm}.chrom.sizes.update") as f:
    for line in f:
        c, s = line.rstrip().split('\t')
        chr_size[c] = int(s)

########################################
# Dinucleotide model
########################################

def prepare_dinuc_probabilities(seq):
    seq = seq.upper()
    counts = Counter()
    total = 0
    for a, b in zip(seq[:-1], seq[1:]):
        if a in BASE_TO_INT and b in BASE_TO_INT:
            counts[a+b] += 1
            total += 1
    dinuc_freq = {
        a+b: counts.get(a+b, 0) / total
        for a in BASES for b in BASES
    }
    mono = np.zeros(4)
    for i, a in enumerate(BASES):
        mono[i] = sum(dinuc_freq[a+b] for b in BASES)
    mono /= mono.sum()
    cond = np.zeros((4,4))
    for i,a in enumerate(BASES):
        row = np.array([dinuc_freq[a+b] for b in BASES])
        cond[i] = row / row.sum() if row.sum() > 0 else 0.25
    return mono, cond

########################################
# FAST dinucleotide sequence generator
########################################

def generate_dinuc_sequence(length, first_probs, cond_mat):
    if length == 0:
        return ""
    seq = np.empty(length, dtype=np.int8)
    seq[0] = np.random.choice(4, p=first_probs)
    r = np.random.rand(length - 1)
    cdf = np.cumsum(cond_mat, axis=1)
    for i in range(1, length):
        seq[i] = np.searchsorted(cdf[seq[i-1]], r[i-1])
    return ''.join(BASES[seq])

########################################
# Precompute replacement schedule
########################################

replacement_schedule = []
left = right = 0
while left + right < seq_length:
    if left == right:
        left += step_size
    else:
        right += step_size
    if left + right > seq_length:
        break
    replacement_schedule.append((left, right))

########################################
# Load regions
########################################

all_regions = []
bed_path = f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/{celltype}.bed"

with open(bed_path) as f:
    for line in f:
        l = line.rstrip().split('\t')
        loc = (int(l[1]) + int(l[2])) // 2
        if loc-1057 > 0 and loc+1057 < chr_size[l[0]]:
            all_regions.append((l[0], int(l[1]), int(l[2])))

with gzip.open(f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region/{celltype}_loc.gz", "wt") as out_loc:
    for chrom, start, end in all_regions:
        out_loc.write(f"{chrom}\t{start}\t{end}\n")

########################################
# Main loop
########################################

batch_size = 100
batches = [(i, min(i+batch_size, len(all_regions)))
           for i in range(0, len(all_regions), batch_size)]

for batch_idx, (start, end) in enumerate(batches, 1):
    print(f"Batch {batch_idx}/{len(batches)}")
    candidate_sequences = []
    for chrom, s, e in all_regions[start:end]:
        loc = (s + e) // 2
        local_seq = fasta.fetch(chrom, max(1, loc - 2500), min(chr_size[chrom], loc + 2500))
        first_probs, cond_mat = prepare_dinuc_probabilities(local_seq)
        orig_seq = fasta.fetch(chrom, loc-1057, loc+1057)
        for left, right in replacement_schedule:
            middle = orig_seq[left:seq_length-right]
            for _ in range(n_reps):
                lp = generate_dinuc_sequence(left, first_probs, cond_mat)
                rp = generate_dinuc_sequence(right, first_probs, cond_mat)
                candidate_sequences.append(lp + middle + rp)
    predictions = crested.tl.predict(input=candidate_sequences, model=model)
    predictions_reshape = np.median(predictions.reshape(-1, 10, predictions.shape[1]), axis=1)
    out_path = f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region/{celltype}_batch_{batch_idx}.gz"
    with gzip.open(out_path, 'wt') as out:
        np.savetxt(out, predictions_reshape, fmt="%.3f")







################################################################################
### step-2: the second round of identifying core region

import random
import sys, os
import numpy as np
import gzip
import pysam
from collections import Counter

import anndata as ad
import crested
import keras

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"

celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]

celltype = celltype_list[int(sys.argv[1]) - 1]
print(celltype)

########################################
# Constants
########################################

BASES = np.array(['A','C','G','T'])
BASE_TO_INT = {b:i for i,b in enumerate(BASES)}

mamm = "Mus_musculus"
model_id = "mouse_fake_track_15"
step_size = 20
seq_length = 2114
n_reps = 10

adata = ad.read_h5ad(f"{work_path}/{model_id}/data_window_cluster_mouse.h5ad")
celltype_index = np.where(np.isin(adata.obs.index, celltype))[0]

########################################
# FASTA + model
########################################

fasta = pysam.FastaFile(f"{work_path}/genome/{mamm}/{mamm}.fa")

model_path = f"{work_path}/{model_id}/mamm_32/window_cluster/finetuned_model_1e5/checkpoints/02.keras"
model = keras.models.load_model(model_path, compile=False)

########################################
# Chrom sizes
########################################

chr_size = {}
with open(f"{work_path}/genome/{mamm}/{mamm}.chrom.sizes.update") as f:
    for line in f:
        c, s = line.rstrip().split('\t')
        chr_size[c] = int(s)

########################################
# Dinucleotide model
########################################

def prepare_dinuc_probabilities(seq):
    seq = seq.upper()
    counts = Counter()
    total = 0
    for a, b in zip(seq[:-1], seq[1:]):
        if a in BASE_TO_INT and b in BASE_TO_INT:
            counts[a+b] += 1
            total += 1
    dinuc_freq = {
        a+b: counts.get(a+b, 0) / total
        for a in BASES for b in BASES
    }
    mono = np.zeros(4)
    for i, a in enumerate(BASES):
        mono[i] = sum(dinuc_freq[a+b] for b in BASES)
    mono /= mono.sum()
    cond = np.zeros((4,4))
    for i,a in enumerate(BASES):
        row = np.array([dinuc_freq[a+b] for b in BASES])
        cond[i] = row / row.sum() if row.sum() > 0 else 0.25
    return mono, cond

########################################
# FAST dinucleotide sequence generator
########################################

def generate_dinuc_sequence(length, first_probs, cond_mat):
    if length == 0:
        return ""
    seq = np.empty(length, dtype=np.int8)
    seq[0] = np.random.choice(4, p=first_probs)
    r = np.random.rand(length - 1)
    cdf = np.cumsum(cond_mat, axis=1)
    for i in range(1, length):
        seq[i] = np.searchsorted(cdf[seq[i-1]], r[i-1])
    return ''.join(BASES[seq])

########################################
# Precompute replacement schedule
########################################

replacement_schedule = []
left = right = 0
while left + right < seq_length:
    if left == right:
        left += step_size
    else:
        right += step_size
    if left + right > seq_length:
        break
    replacement_schedule.append((left, right))

########################################
# Load regions
########################################

all_regions = []
with gzip.open(f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region/{celltype}_loc.gz", "rt") as file:
    for line in file:
        l = line.rstrip().split('\t')
        all_regions.append((l[0], int(l[1]), int(l[2])))

########################################
# Main loop
########################################

batch_size = 100
batches = [(i, min(i+batch_size, len(all_regions)))
           for i in range(0, len(all_regions), batch_size)]

for batch_idx, (start, end) in enumerate(batches, 1):
    print(f"Batch {batch_idx}/{len(batches)}")
    orig_sequences = []
    for chrom, s, e in all_regions[start:end]:
        loc = (s + e) // 2
        orig_sequences.append(fasta.fetch(chrom, loc-1057, loc+1057))
    orig_sequences_pred = crested.tl.predict(input=orig_sequences, model=model)[:,celltype_index]
    result = np.loadtxt(f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region/{celltype}_batch_{batch_idx}.gz")[:,celltype_index]
    candidate_loc = []
    candidate_sequences = []
    for cnt, (chrom, s, e) in enumerate(all_regions[start:end], 1):
        loc = (s + e) // 2
        local_seq = fasta.fetch(chrom, max(1, loc - 2500), min(chr_size[chrom], loc + 2500))
        first_probs, cond_mat = prepare_dinuc_probabilities(local_seq)
        orig_seq = fasta.fetch(chrom, loc-1057, loc+1057)
        orig_seq_score = orig_sequences_pred[cnt-1]
        value = result[(cnt-1)*len(replacement_schedule) : cnt*len(replacement_schedule)]
        for replacement in range(len(replacement_schedule)):
            if value[replacement] < orig_seq_score*0.9:
                left, right = replacement_schedule[replacement]
                if left == right:
                    right = max(0, right - step_size)
                    rp = [generate_dinuc_sequence(right, first_probs, cond_mat) for _ in range(n_reps)]
                    for left_ in range(left + step_size, seq_length - right, step_size):
                        candidate_loc.append((chrom, str(s), str(e), str(round(orig_seq_score[0],3)), str(left_), str(right)))
                        middle = orig_seq[left_:seq_length-right]
                        for x in range(n_reps):
                            lp = generate_dinuc_sequence(left_, first_probs, cond_mat)
                            candidate_sequences.append(lp + middle + rp[x])
                else:
                    left = max(0, left - step_size)
                    lp = [generate_dinuc_sequence(left, first_probs, cond_mat) for _ in range(n_reps)]
                    for right_ in range(right + step_size, seq_length - left, step_size):
                        candidate_loc.append((chrom, str(s), str(e), str(round(orig_seq_score[0],3)), str(left), str(right_)))
                        middle = orig_seq[left:seq_length-right_]
                        for x in range(n_reps):
                            rp = generate_dinuc_sequence(right_, first_probs, cond_mat)
                            candidate_sequences.append(lp[x] + middle + rp)
                break
    predictions = crested.tl.predict(input=candidate_sequences, model=model)
    predictions_reshape = np.median(predictions.reshape(-1, n_reps, predictions.shape[1]), axis=1)
    out_path = f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region_2/{celltype}_batch_{batch_idx}.gz"
    with gzip.open(out_path, 'wt') as out:
        np.savetxt(out, predictions_reshape, fmt="%.3f")
    out_path = f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region_2/{celltype}_batch_{batch_idx}_loc.gz"
    with gzip.open(out_path, 'wt') as out:
        np.savetxt(out, candidate_loc, fmt='%s', delimiter='\t')








######################################
### step-3: identifying core regions

import random
import sys, os
import numpy as np
import pandas as pd
import gzip
import pysam
from collections import Counter

import anndata as ad

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"

celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]

celltype = celltype_list[int(sys.argv[1]) - 1]
print(celltype)

mamm = "Mus_musculus"
model_id = "mouse_fake_track_15"
step_size = 20
seq_length = 2114

adata = ad.read_h5ad(f"{work_path}/{model_id}/data_window_cluster_mouse.h5ad")
celltype_index = np.where(np.isin(adata.obs.index, celltype))[0]

replacement_schedule = []
left = right = 0
while left + right < seq_length:
    if left == right:
        left += step_size
    else:
        right += step_size
    if left + right > seq_length:
        break
    replacement_schedule.append((left, right))

all_regions = []
with gzip.open(f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region/{celltype}_loc.gz", "rt") as file:
    for line in file:
        l = line.rstrip().split('\t')
        all_regions.append((l[0], int(l[1]), int(l[2])))

batch_size = 100
batches = [(i, min(i+batch_size, len(all_regions)))
           for i in range(0, len(all_regions), batch_size)]
batch_num = len(batches)

dat_loc_old = []
for batch_idx, (start, end) in enumerate(batches, 1):
    print(f"Batch {batch_idx}/{len(batches)}")
    for chrom, s, e in all_regions[start:end]:
        for replacement in range(len(replacement_schedule)):
            left, right = replacement_schedule[replacement]
            dat_loc_old.append((chrom, s, e, left, right))

dat_old = []
for batch_idx in range(batch_num):
    print(f"{batch_idx+1}/{batch_num}")
    arr = np.loadtxt(f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region/{celltype}_batch_{batch_idx+1}.gz")[:, celltype_index]
    dat_old.append(arr)

dat_loc_old = pd.DataFrame(dat_loc_old, columns=["chrom", "start", "end", "left", "right"])
dat_loc_old['score'] = np.concatenate(dat_old)

dat_loc = []
dat = []
for batch_idx in range(batch_num):
    print(f"{batch_idx+1}/{batch_num}")
    df = pd.read_csv(
        f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region_2/{celltype}_batch_{batch_idx+1}_loc.gz",
        sep="\t",
        header=None)
    dat_loc.append(df)
    arr = np.loadtxt(f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region_2/{celltype}_batch_{batch_idx+1}.gz")[:, celltype_index]
    dat.append(arr)

dat_loc = pd.concat(dat_loc, ignore_index=True)
dat_loc['score'] = np.concatenate(dat)
dat_loc.columns = ['chrom', 'start', 'end', 'orig_score', 'left', 'right', 'score']

result = []
result_err = []
for cnt, (chrom, s, e) in enumerate(all_regions, 1):
    print(f"{cnt}/{len(all_regions)}")
    loc = (s + e) // 2
    dat = dat_loc[(dat_loc['chrom'] == chrom) & 
                  (dat_loc['start'] == s) & 
                  (dat_loc['end'] == e)]
    dat_old = dat_loc_old[(dat_loc_old['chrom'] == chrom) & 
                          (dat_loc_old['start'] == s) & 
                          (dat_loc_old['end'] == e)]
    dat_all = pd.concat([dat.drop(columns=['orig_score']), dat_old], axis=0, ignore_index=True)
    dat_filtered = dat[dat['score'] < dat['orig_score'] * 0.9]
    if len(dat) == 0:
        result_err.append(dat_old)
    else:
        if len(dat_filtered) == 0:
            row = dat.iloc[0]
        else:
            row = dat_filtered.iloc[0]
        if row['left'] > row['right']:
            left = row['left']-step_size
            right = row['right']
        else:
            left = row['left']
            right = row['right']-step_size
        if left == 0 and right == 0:
            result_err.append(dat_old)
        else:
            score = dat_all[(dat_all['left'] == left) & (dat_all['right'] == right)]['score'].iloc[0]
            result.append((chrom, s, e, left, right, loc-1057+left, loc+1057-right, score))


out_path = f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region_3/{celltype}_result.gz"
with gzip.open(out_path, 'wt') as out:
    np.savetxt(out, result, fmt='%s', delimiter='\t')

if len(result_err) != 0:
    out_path = f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/core_region_3/{celltype}_result.err.gz"
    with gzip.open(out_path, 'wt') as out:
        np.savetxt(out, pd.concat(result_err, axis=0, ignore_index=True), fmt='%s', delimiter='\t')







##############################################
### merging cell types to create a single file

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")
options(scipen = 999)

model_id = "mouse_fake_track_15"

mamm = "Mus_musculus"

celltype_list=c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

dat = list()
for(celltype in celltype_list){
    print(celltype)
    dat_i = read.table(paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/core_region_3/', celltype, '_result.gz'))
    colnames(dat_i) = c('chrom', 'start', 'end', 'left', 'right', 'core_start', 'core_end', 'score')
    dat_i = dat_i[,c('chrom', 'core_start', 'core_end', 'score')]
    dat_i$celltype = celltype
    dat[[celltype]] = dat_i
    dat_i$peak_id = paste0("peak_", 1:nrow(dat_i))
    dat_i$strand = 0
    write.table(dat_i[,c("peak_id", "chrom", "core_start", "core_end", "strand")], 
        paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/core_region_3/', celltype, '.txt'), row.names=F, col.names=F, sep='\t', quote=F)
}
dat = do.call(rbind, dat)
rownames(dat) = NULL

write.table(dat, paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/window_list_10Mb.core_region.bed'), row.names=F, col.names=F, sep='\t', quote=F)
### n = 341552
### avg size = 470, median size = 414; sd = 288

cut -f1-3 window_list_10Mb.core_region.bed \
  | bedtools sort -i - \
  | bedtools merge -i - > window_list_10Mb.core_region.merge_across_celltype.bed
### 84,870,322 bps footprints
### 84870322/2818974548 = 3.0%



########################################################################
### as jay suggested, converting the predicted score to Phred-like value

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

model_id = "mouse_fake_track_15"

mamm = "Mus_musculus"

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
celltype = celltype_list[kk]

dat = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "_trf/dat.txt.rds"))

dat_i = read.table(paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/core_region_3/', celltype, '_result.gz'))
colnames(dat_i) = c('chrom', 'start', 'end', 'left', 'right', 'core_start', 'core_end', 'score')
dat_i = dat_i[,c('chrom', 'core_start', 'core_end', 'score')]

values = dat[, kk]
sorted = sort(values)
N = length(sorted)

Q_from_x <- function(x) {
    r <- findInterval(x, sorted)
    q <- (r - 0.5) / N
    q <- min(q, 1 - 1e-12)
    -10 * log10(1 - q)
}

dat_i$Phred_score = round(vapply(dat_i$score, Q_from_x, numeric(1)), 3)
dat_i$score = NULL
dat_i$celltype = celltype

options(scipen = 999)
write.table(dat_i, paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/core_region_3/', celltype, '_result.Phred.txt'), row.names=F ,col.names=F, sep="\t", quote=F)


### merge the result
celltype_list="Adipocyte_cells Adipocyte_cells_Cyp2e1 B_cells Brain_capillary_endothelial_cells CNS_neurons Cardiomyocytes Corticofugal_neurons Endocardial_cells Endothelium Epithelial_cells Erythroid_cells Eye Glia Glomerular_endothelial_cells Gut_epithelial_cells Hepatocytes Intermediate_neuronal_progenitors Kidney Lateral_plate_and_intermediate_mesoderm Liver_sinusoidal_endothelial_cells Lung_and_airway Lymphatic_vessel_endothelial_cells Melanocyte_cells Mesoderm Neural_crest_PNS_neurons Neuroectoderm_and_glia Olfactory_ensheathing_cells Olfactory_neurons Oligodendrocytes Skeletal_muscle_cells T_cells White_blood_cells"
for celltype in ${celltype_list}; do
    cat "${celltype}_result.Phred.txt" >> ../window_list_10Mb.core_region.phred_score.bed
done




### combining predicted score and Phred-like score, to create a bed file for downloading
### Supplementary_File_2_Peaks_Specific_to_32_Clusters.bed.gz
mamm = "Mus_musculus"
dat = read.table(paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/window_list_10Mb.core_region.bed'))
dat_x = read.table(paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/window_list_10Mb.core_region.phred_score.bed'))

dat$V6 = as.vector(dat_x$V4)
colnames(dat) = c("chr", "start", "end", "orig_pred_score", "cell_class", "phred_like_score")
dat = dat[,c("chr", "start", "end", "orig_pred_score", "phred_like_score", "cell_class")]

options(scipen = 999)
write.table(dat, "/net/shendure/vol10/www/content/members/cxqiu/public/backup/jax_atac/download/Supplementary_File_3_Evolution_Augmented_Model_Predict_On_Mouse_Genome.bed", row.names=F, col.names=T, sep="\t", quote=F)




