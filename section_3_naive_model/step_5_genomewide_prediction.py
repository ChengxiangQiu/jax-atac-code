
#######################################
### Preparing the genomewide prediction

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu

### The evolution naive model is available:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/evolution_naive.keras

###################################################
### Step-1: Genomewide prediction using naive model

import sys, os
import anndata as ad
import crested
import numpy as np
import matplotlib
import gzip
import keras
import pysam

work_path = "XXX"

with open("mamm_list.txt") as f:
    mamm_list = [line.rstrip() for line in f if line.strip()]

mamm = mamm_list[int(sys.argv[1]) - 1]
print(mamm)

# load a trained model
model_path = "evolution_naive.keras"
model = keras.models.load_model(model_path, compile=False)  # change to your model path

output_path = os.path.join(work_path, "celltype_L2_cut_norm", "prediction_mammals", f"prediction_{mamm}")
if not os.path.exists(output_path):
    os.makedirs(output_path)

fasta = pysam.FastaFile(os.path.join(work_path, "genome", mamm, f"{mamm}.fa"))

window_size = 2114
step_size = 100
central_size = 1000
offset = (window_size - central_size) // 2

all_regions = []
central_start = []

with gzip.open(os.path.join(work_path, "genome", mamm, "candidate_windows.txt.gz"), 'rt') as f:
    for line in f:
        chrom, pos_str = line.rstrip().split('\t')
        pos = int(pos_str)
        all_regions.append(f"{chrom}:{pos}-{pos + window_size}")
        central_start.append((chrom, pos + offset))

batch_size = 100000
batch_list = [(start, min(start + batch_size, len(all_regions))) for start in range(0, len(all_regions), batch_size)]
print(f"Number of batches: {len(batch_list)}")

batch_num = 1
for batch in batch_list:
    print(f"Processing batch {batch_num}/{len(batch_list)}: {batch[0]} - {batch[1]}")

    n_counts = []
    candidate_sequences = []
    for region in all_regions[batch[0]:batch[1]]:
        chrom, coords = region.split(":", 1)
        start, end = map(int, coords.split("-"))
        seq = fasta.fetch(chrom, start, end)
        candidate_sequences.append(seq)
        n_count = sum(seq.upper().count(b) for b in "AGCT")
        n_counts.append(n_count)

    nonzero_counts = np.array(n_counts).reshape(-1, 1)
    file_name = os.path.join(output_path, f"batch_{batch_num}.nonzero.txt.gz")
    with gzip.open(file_name, 'wt') as f:
        np.savetxt(f, nonzero_counts, fmt="%d")

    predictions = crested.tl.predict(input=candidate_sequences, model=model)
    file_name = os.path.join(output_path, f"batch_{batch_num}.txt.gz")
    with gzip.open(file_name, 'wt') as f:
        np.savetxt(f, predictions, fmt="%.3f")

    arr = np.array(central_start[batch[0]:batch[1]], dtype=object)
    file_name = os.path.join(output_path, f"batch_{batch_num}.loc.txt.gz")
    with gzip.open(file_name, 'wt') as f:
        np.savetxt(f, arr, fmt='%s\t%d')

    batch_num += 1

print("Completing analysis!")



###########################################################################
### Step-2: Merging batches, calculating 100-bp resolution predicted scores

import gzip
import os, sys
import numpy as np
import pandas as pd

with open("mamm_list.txt") as f:
    mamm_list = [line.rstrip() for line in f if line.strip()]

mamm = mamm_list[int(sys.argv[1]) - 1]
print(mamm)

work_path = ""

data_path = os.path.join(work_path, "celltype_L2_cut_norm", "prediction_mammals", f"prediction_{mamm}")

file_number = len([i for i in os.listdir(data_path) if "loc" in i and "batch" in i])
file_list = [f"batch_{i}.txt.gz" for i in range(1, file_number + 1)]
file_loc_list = [f"batch_{i}.loc.txt.gz" for i in range(1, file_number + 1)]
file_nonzero_list = [f"batch_{i}.nonzero.txt.gz" for i in range(1, file_number + 1)]

dat = []
dat_loc = []
dat_nonzero = []

for batch in range(file_number):
    print(f"Reading batch: {batch + 1}/{file_number}")
    with gzip.open(os.path.join(data_path, file_list[batch]), 'rt') as f:
        dat.append(np.loadtxt(f))
    with gzip.open(os.path.join(data_path, file_loc_list[batch]), 'rt') as f:
        dat_loc.extend(tuple(line.rstrip().split('\t')[:2]) for line in f)
    with gzip.open(os.path.join(data_path, file_nonzero_list[batch]), 'rt') as f:
        dat_nonzero.append(np.loadtxt(f).reshape(-1, 1))

dat = np.vstack(dat)
dat_nonzero = np.vstack(dat_nonzero)

chrom = dat_loc[0][0]
start_index = 0

with gzip.open(os.path.join(data_path, 'dat.txt.gz'), 'wt') as f, gzip.open(os.path.join(data_path, 'dat_loc.txt.gz'), 'wt') as f_loc:
    for cnt in range(dat.shape[0]):
        if dat_loc[cnt][0] != chrom:
            dat_sub = dat[start_index:cnt, :]
            dat_loc_sub = dat_loc[start_index:cnt]
            dat_nonzero_sub = dat_nonzero[start_index:cnt, :]
            mask = dat_nonzero_sub[:, 0] < 2000
            dat_sub[mask, :] = 0
            print(f"{chrom}:{dat_sub.shape[0]}")
            for i in range(9, dat_sub.shape[0]):
                window_sum = np.sum(dat_sub[(i-9):(i+1), :], axis=0)/10
                f.write('\t'.join(f'{x:.3f}' for x in window_sum) + '\n')
                f_loc.write('\t'.join(dat_loc_sub[i]) + '\n')
            chrom = dat_loc[cnt][0]
            start_index = cnt

    cnt += 1
    dat_sub = dat[start_index:cnt, :]
    dat_loc_sub = dat_loc[start_index:cnt]
    dat_nonzero_sub = dat_nonzero[start_index:cnt, :]
    mask = dat_nonzero_sub[:, 0] < 2000
    dat_sub[mask, :] = 0
    print(f"{chrom}:{dat_sub.shape[0]}")
    for i in range(9, dat_sub.shape[0]):
        window_sum = np.sum(dat_sub[(i-9):(i+1), :], axis=0)/10
        f.write('\t'.join(f'{x:.3f}' for x in window_sum) + '\n')
        f_loc.write('\t'.join(dat_loc_sub[i]) + '\n')




