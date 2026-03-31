


############################################################
### prediction on mouse/human genome using 100 bp resolution

import sys, os
import anndata as ad
import crested
import numpy as np
import matplotlib
import gzip
import keras
import pysam

work_path = "/gpfs/projects/shendurelabcre/cxqiu/atac_seq/14_crested"
work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
model_id = "mouse_fake_track_15"

mamm = "Homo_sapiens"

iter = int(sys.argv[1])

# load a trained model
model_path = f"{work_path}/{model_id}/mamm_32/window_cluster/finetuned_model_1e5/checkpoints/02.keras"
model = keras.models.load_model(model_path, compile=False)

output_path = os.path.join(f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf/trf_{iter}")
if not os.path.exists(output_path):
    os.makedirs(output_path)

fasta = pysam.FastaFile(f"{work_path}/genome/{mamm}/fasta_trf/{mamm}_trf_{iter}.fa")

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



#############
### merge ###

import gzip
import os, sys
import numpy as np
import pandas as pd

work_path = "/gpfs/projects/shendurelabcre/cxqiu/atac_seq/14_crested"
model_id = "mouse_fake_track_15"

mamm = sys.argv[1]

data_path = f"{work_path}/{model_id}/prediction_mammals/prediction_{mamm}_trf"

file_number = len([i for i in os.listdir(f"{data_path}/trf_1") if "loc" in i and "batch" in i])
file_list = [f"batch_{i}.txt.gz" for i in range(1, file_number + 1)]
file_loc_list = [f"batch_{i}.loc.txt.gz" for i in range(1, file_number + 1)]
file_nonzero_list = [f"batch_{i}.nonzero.txt.gz" for i in range(1, file_number + 1)]

dat_loc = []
dat_nonzero = []

for batch in range(file_number):
    print(f"Reading batch: {batch + 1}/{file_number}")
    with gzip.open(os.path.join(f"{data_path}/trf_1", file_loc_list[batch]), 'rt') as f:
        dat_loc.extend(tuple(line.rstrip().split('\t')[:2]) for line in f)
    with gzip.open(os.path.join(f"{data_path}/trf_1", file_nonzero_list[batch]), 'rt') as f:
        dat_nonzero.append(np.loadtxt(f).reshape(-1, 1))

dat_nonzero = np.vstack(dat_nonzero)

dat = []
for batch in range(file_number):
    print(f"Reading batch: {batch + 1}/{file_number}")
    dat_batch = []
    for iter in range(1, 11):
        with gzip.open(os.path.join(f"{data_path}/trf_{iter}", file_list[batch]), 'rt') as f:
            dat_batch.append(np.loadtxt(f))
    dat_batch_arr = np.stack(dat_batch, axis = 0)
    dat_batch_mean = np.mean(dat_batch_arr, axis = 0)
    dat.append(dat_batch_mean)

dat = np.vstack(dat)

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


