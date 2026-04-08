
############################################################
### Identifying enhancers predicted by STEAM-v1-model (hg38)

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


#######################################################
### Step-1: convert predicted score to Phred like value

import gzip
import os, sys
import numpy as np

model_id = "mouse_fake_track_15"
mamm = "Homo_sapiens"

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
model_id = "mouse_fake_track_15"
mamm = "Homo_sapiens"
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



#########################################
### step-2: Identifying candidate windows

source("~/work/scripts/utils.R")
library(GenomicRanges)

model_id = "mouse_fake_track_15"

mamm = "Homo_sapiens"

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

[1] "Adipocyte_cells: 6.736"
[1] "Adipocyte_cells_Cyp2e1: 6.35"
[1] "B_cells: 2.928"
[1] "Brain_capillary_endothelial_cells: 4.208"
[1] "CNS_neurons: 9.741"
[1] "Cardiomyocytes: 5.061"
[1] "Corticofugal_neurons: 6.778"
[1] "Endocardial_cells: 4.929"
[1] "Endothelium: 4.89300307156145"
[1] "Epithelial_cells: 7.838"
[1] "Erythroid_cells: 8.605"
[1] "Eye: 5.249"
[1] "Glia: 6.108"
[1] "Glomerular_endothelial_cells: 3.687"
[1] "Gut_epithelial_cells: 5.732"
[1] "Hepatocytes: 5.937"
[1] "Intermediate_neuronal_progenitors: 8.841"
[1] "Kidney: 6.581"
[1] "Lateral_plate_and_intermediate_mesoderm: 6.174"
[1] "Liver_sinusoidal_endothelial_cells: 4.19000307156146"
[1] "Lung_and_airway: 5.587"
[1] "Lymphatic_vessel_endothelial_cells: 4.251"
[1] "Melanocyte_cells: 3.627"
[1] "Mesoderm: 7.075"
[1] "Neural_crest_PNS_neurons: 9.471"
[1] "Neuroectoderm_and_glia: 6.514"
[1] "Olfactory_ensheathing_cells: 4.495"
[1] "Olfactory_neurons: 8.388"
[1] "Oligodendrocytes: 5.821"
[1] "Skeletal_muscle_cells: 6.759"
[1] "T_cells: 4.905"
[1] "White_blood_cells: 8.69100307156146"

x = window_list %>% group_by(celltype) %>% slice_min(order_by = phred_score, n = 1, with_ties = F)
### Q = 25.1
### q = 0.99691, top 0.309%

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

### n = 336,649

### this is after merging 100 bp windows, but before detecting core regions by left-right trimming







##############################################################################
### step-3: the first round of identifying core region

import random
import sys, os
import numpy as np
import gzip
import pysam
from collections import Counter

import anndata as ad
import crested
import keras


celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]

celltype = celltype_list[int(sys.argv[1]) - 1]
print(celltype)

########################################
# Constants
########################################

BASES = np.array(['A','C','G','T'])
BASE_TO_INT = {b:i for i,b in enumerate(BASES)}

mamm = "Homo_sapiens"
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
### step-4: the second round of identifying core region

import random
import sys, os
import numpy as np
import gzip
import pysam
from collections import Counter

import anndata as ad
import crested
import keras


celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]

celltype = celltype_list[int(sys.argv[1]) - 1]
print(celltype)

########################################
# Constants
########################################

BASES = np.array(['A','C','G','T'])
BASE_TO_INT = {b:i for i,b in enumerate(BASES)}

mamm = "Homo_sapiens"
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
### step-5: identifying core regions

import random
import sys, os
import numpy as np
import pandas as pd
import gzip
import pysam
from collections import Counter

import anndata as ad

celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]

celltype = celltype_list[int(sys.argv[1]) - 1]
print(celltype)

mamm = "Homo_sapiens"
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



######################################################
### Step-6: merging cell types to create a single file


options(scipen = 999)

model_id = "mouse_fake_track_15"

mamm = "Homo_sapiens"

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
### n = 336378
### avg size = 542, median size = 494; sd = 310

cut -f1-3 window_list_10Mb.core_region.bed \
  | bedtools sort -i - \
  | bedtools merge -i - > window_list_10Mb.core_region.merge_across_celltype.bed
### 90,044,608 bps footprints
### 90044608/3257347282 = 2.8%


bedtools intersect -a window_list_10Mb.core_region.merge_across_celltype.bed \
-b window_list_10Mb.core_region.bed \
-wa -wb > window_list_10Mb.core_region.merge_across_celltype.overlap_back.bed

bedtools intersect -a window_list_10Mb.core_region.merge_across_celltype.bed \
-b /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/TSS_2500.human_v25.bed \
-wa | uniq > window_list_10Mb.core_region.merge_across_celltype.overlap_promoter.bed



x = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/prediction_mammals/prediction_Homo_sapiens_trf/window_list_10Mb.core_region.merge_across_celltype.overlap_back.bed"))
y = unique(x[,c(1,2,3,8)]) %>% group_by(V1, V2, V3) %>% tally()
z = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/prediction_mammals/prediction_Homo_sapiens_trf/window_list_10Mb.core_region.merge_across_celltype.overlap_promoter.bed"))
z$category = "overlap_promoter"
print(nrow(y))
print(round(100*sum(y$n == 1)/nrow(y)))
print(round(100*sum(y$n == 2)/nrow(y)))
print(round(100*sum(y$n > 2)/nrow(y)))

[1] 158047
[1] 62
[1] 21
[1] 17

y_num_human = table(y$n)
y_num_mouse = readRDS("/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/prediction_mammals/prediction_Mus_musculus_trf/mouse_enhancer_species_num.rds")

df = data.frame(num_celltypes = c(1:32, 1:32),
    num_enhancers = c(as.vector(y_num_mouse), as.vector(y_num_human)),
    species = c(rep("mouse", 32), rep("human", 32)))
df$num_celltypes[df$num_celltypes > 10] = 10

df$num_celltypes = factor(df$num_celltypes, levels =1:10)
df$species = factor(df$species, levels = c("mouse", "human"))

p = ggplot(data=df, aes(x=num_celltypes, y=num_enhancers, fill = species)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic(base_size = 10) +
  scale_fill_manual(values=c("human" = "#658bca", "mouse" = "#d79c36"))
ggsave("~/share/hist_num_celltypes.pdf", p, width = 5, height = 5)



y = y %>% left_join(z, by = c("V1", "V2", "V3"))
y$category[is.na(y$category)] = "no_overlap"
y_x = y
colnames(y_x) = c("chr","start","end","num_celltype", "category")
y_x$num_celltype[y_x$num_celltype > 10] = 10
y_x = y_x %>% filter(category == "overlap_promoter") %>% group_by(num_celltype) %>% tally() %>%
left_join(y_x %>% group_by(num_celltype) %>% tally() %>% rename(total_n = n), by = "num_celltype") %>%
mutate(frac = n/total_n, species = "human")


y_x_mouse = readRDS("/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/prediction_mammals/prediction_Mus_musculus_trf/mouse_enhancer_overlap_promoter.rds")

y_x = rbind(y_x, y_x_mouse)
y_x$frac = 100 * y_x$frac

y_x$num_celltype = factor(y_x$num_celltype, levels =1:10)
y_x$species = factor(y_x$species, levels = c("mouse", "human"))


python /net/gs/vol1/home/cxqiu/bin/python_script/generate_background_region.py \
    --n 10000 \
    --mean 569.7 \
    --sd 355.9 \
    --chrom_sizes /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/genome/Homo_sapiens/Homo_sapiens.chrom.sizes.update \
    --out /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/background_region.human.bed

bedtools intersect -a background_region.human.bed -b TSS_2500.human_v25.bed -wa | uniq | wc -l
### 0.1009

p = ggplot(data=y_x, aes(x=num_celltype, y=frac, fill = species)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_hline(yintercept = 9) +
  geom_hline(yintercept = 10) +
  theme_classic(base_size = 10) +
  scale_fill_manual(values=c("human" = "#658bca", "mouse" = "#d79c36"))
ggsave("~/share/hist_frac_promoter.pdf", p, width = 5, height = 5)






celltype_list=c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

x_sub = unique(x[,c(1,2,3,8)])

y_sub = y[y$n <= 6,]

x_sub = x_sub %>% left_join(y_sub, by = c("V1", "V2", "V3")) %>% filter(!is.na(n))


df = NULL
for(celltype_i in celltype_list){
    print(celltype_i)
    for(celltype_j in celltype_list){
        x_i = x_sub %>% filter(V8 %in% c(celltype_i, celltype_j)) %>%
        group_by(V1, V2, V3) %>% tally()
        df = rbind(df, data.frame(celltype_i, celltype_j, sum(x_i$n > 1), nrow(x_i)))
    }
}

colnames(df) = c("celltype_i", "celltype_j", "overlap_n", "union_n")
df$jacc = df$overlap_n/df$union_n

mat = df[,c("celltype_i","celltype_j","jacc")] %>% dcast(celltype_i~celltype_j)
rownames(mat) = mat[,1]
mat = as.matrix(mat[,-1])

diag(mat) = max(mat)

celltype_list_order = c("Gut_epithelial_cells","Hepatocytes","Kidney","Lung_and_airway","Epithelial_cells","Lateral_plate_and_intermediate_mesoderm","Mesoderm","Skeletal_muscle_cells","Cardiomyocytes","Erythroid_cells","Neuroectoderm_and_glia","Eye","CNS_neurons","Neural_crest_PNS_neurons","Intermediate_neuronal_progenitors","B_cells","T_cells","White_blood_cells","Adipocyte_cells_Cyp2e1","Adipocyte_cells","Glia","Melanocyte_cells","Olfactory_ensheathing_cells","Oligodendrocytes","Olfactory_neurons","Corticofugal_neurons","Endothelium","Glomerular_endothelial_cells","Brain_capillary_endothelial_cells","Liver_sinusoidal_endothelial_cells","Endocardial_cells","Lymphatic_vessel_endothelial_cells")

library(gplots)
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)

bk <- seq(0, 0.66, length.out = 101)

pdf("~/share/human_enhancers_overlap.pdf", 8, 8)
p = heatmap.2(as.matrix(mat[celltype_list_order, celltype_list_order]), 
          #col=Colors, 
          col=viridis(100), 
          breaks = bk,
          scale="none", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(5,5))
dev.off()



##############################################################
### Step-7: converting the predicted score to Phred-like value


model_id = "mouse_fake_track_15"

mamm = "Homo_sapiens"

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

mamm = "Homo_sapiens"
dat = read.table(paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/window_list_10Mb.core_region.bed'))
dat_x = read.table(paste0(work_path, '/14_crested/', model_id, '/prediction_mammals/prediction_', mamm, '_trf/window_list_10Mb.core_region.phred_score.bed'))

dat$V6 = as.vector(dat_x$V4)
colnames(dat) = c("chr", "start", "end", "orig_pred_score", "cell_class", "phred_like_score")
dat = dat[,c("chr", "start", "end", "orig_pred_score", "phred_like_score", "cell_class")]

options(scipen = 999)
write.table(dat, "/net/shendure/vol10/www/content/members/cxqiu/public/backup/jax_atac/download/Supplementary_File_4_Evolution_Augmented_Model_Predict_On_Human_Genome.bed", row.names=F, col.names=T, sep="\t", quote=F)



######################################
### Step-8: using HOMER to call motifs



### identify enriched TF motifs by HOMER

celltype_list=(Adipocyte_cells Adipocyte_cells_Cyp2e1 B_cells Brain_capillary_endothelial_cells CNS_neurons Cardiomyocytes Corticofugal_neurons Endocardial_cells Endothelium Epithelial_cells Erythroid_cells Eye Glia Glomerular_endothelial_cells Gut_epithelial_cells Hepatocytes Intermediate_neuronal_progenitors Kidney Lateral_plate_and_intermediate_mesoderm Liver_sinusoidal_endothelial_cells Lung_and_airway Lymphatic_vessel_endothelial_cells Melanocyte_cells Mesoderm Neural_crest_PNS_neurons Neuroectoderm_and_glia Olfactory_ensheathing_cells Olfactory_neurons Oligodendrocytes Skeletal_muscle_cells T_cells White_blood_cells)
celltype="${celltype_list[$SGE_TASK_ID - 1]}"
echo $celltype
model_id=mouse_fake_track_15
work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/"$model_id"/prediction_mammals/prediction_Homo_sapiens_trf/core_region_3
/net/shendure/vol10/projects/cxqiu/nobackup/install/HOMER/bin/findMotifsGenome.pl \
"$work_path"/$celltype.txt \
hg38 \
"$work_path"/"$celltype"_denovo \
-noknown -cpg



### mv results to the web server

celltype_list=(Adipocyte_cells Adipocyte_cells_Cyp2e1 B_cells Brain_capillary_endothelial_cells CNS_neurons Cardiomyocytes Corticofugal_neurons Endocardial_cells Endothelium Epithelial_cells Erythroid_cells Eye Glia Glomerular_endothelial_cells Gut_epithelial_cells Hepatocytes Intermediate_neuronal_progenitors Kidney Lateral_plate_and_intermediate_mesoderm Liver_sinusoidal_endothelial_cells Lung_and_airway Lymphatic_vessel_endothelial_cells Melanocyte_cells Mesoderm Neural_crest_PNS_neurons Neuroectoderm_and_glia Olfactory_ensheathing_cells Olfactory_neurons Oligodendrocytes Skeletal_muscle_cells T_cells White_blood_cells)
model_id=mouse_fake_track_15
work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/"$model_id"/prediction_mammals/prediction_Homo_sapiens_trf/core_region_3
web_path=/net/shendure/vol10/www/content/members/cxqiu/public/nobackup/HOMER_motif/human_enhancers
for celltype in "${celltype_list[@]}"; do
cp $work_path/"$celltype"_denovo/homerResults.html $web_path/$celltype.html
done



