

######################################################################
### Genome-wide prediction with the tandem repeat mitigation heuristic

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


######################################################################################################
### Step-1: creating new mouse genome by replacing the TRF regions with locally dinu-matched sequences

### mm10_SimpleRepeats.bed was downloaded from UCSC table brower (mm10)
### bedtools sort -i mm10_SimpleRepeats.bed | bedtools merge -i - > mm10_SimpleRepeats.merged.bed

import pysam
import numpy as np
import os, sys
from collections import Counter

iter = int(sys.argv[1])

work_path = ""

mamm = "Mus_musculus"

if mamm == "Mus_musculus":
    mamm_short = "mm10"
else:
    mamm_short = "hg38"

chr_list = set()
with open(f"{work_path}/genome/{mamm}/{mamm}.chrom.sizes.update") as file:
    for line in file:
        l = line.rstrip().split('\t')
        chr_list.add(l[0])

def dinuc_freq_from_sequence(seq):
    seq = seq.upper()
    counts = Counter()
    total = 0
    for a, b in zip(seq[:-1], seq[1:]):
        if a in "ACGT" and b in "ACGT":
            counts[a + b] += 1
            total += 1
    all_dinucs = [a + b for a in "ACGT" for b in "ACGT"]
    freqs = {d: counts.get(d, 0) / total for d in all_dinucs}
    return freqs

def prepare_dinuc_probabilities(seq):
    dinuc_freq = dinuc_freq_from_sequence(seq)
    bases = np.array(['A','C','G','T'])
    mono_counts = {b:0 for b in bases}
    for d, f in dinuc_freq.items():
        mono_counts[d[0]] += f
    first_base_probs = np.array([mono_counts[b] for b in bases])
    first_base_probs /= first_base_probs.sum()
    cond_probs = {}
    for b1 in bases:
        probs = np.array([dinuc_freq[b1+b2] for b2 in bases])
        total = probs.sum()
        if total == 0:
            probs = np.ones(4)/4
        else:
            probs /= total
        cond_probs[b1] = probs
    return first_base_probs, cond_probs

def generate_dinuc_sequence(length, first_base_probs, cond_probs, seed=None):
    if seed is not None:
        np.random.seed(seed)
    bases = np.array(['A','C','G','T'])
    seq = []
    first_base = np.random.choice(bases, p=first_base_probs)
    seq.append(first_base)
    for _ in range(length-1):
        prev_base = seq[-1]
        next_base = np.random.choice(bases, p=cond_probs[prev_base])
        seq.append(next_base)
    return "".join(seq)

fa = pysam.FastaFile(f"{work_path}/genome/{mamm}/{mamm}.fa")

seqs = {name: list(fa.fetch(name)) for name in fa.references if name in chr_list}

seqs_length = {name: len(seqs[name]) for name in fa.references if name in chr_list}

with open(f"{work_path}/genome/{mamm}/{mamm_short}_SimpleRepeats.merged.bed") as bed:
    for line in bed:
        chrom, start, end = line.strip().split()[:3]
        if chrom in chr_list:
            start, end = int(start), int(end)
            length = end - start
            middle = (start + end)//2
            seq_local = fa.fetch(chrom, max(1, middle - 2500), min(seqs_length[chrom], middle + 2500))
            first_base_probs, cond_probs = prepare_dinuc_probabilities(seq_local)
            rand_bases = generate_dinuc_sequence(length, first_base_probs, cond_probs)
            seqs[chrom][start:end] = rand_bases

# Write masked fasta
with open(f"{work_path}/genome/{mamm}/fasta_trf/{mamm}_trf_{iter}.fa", "w") as out:
    for chrom in chr_list:
        out.write(f">{chrom}\n")
        seq = "".join(seqs[chrom])
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")





####################################################################
### Step-2: prediction on mouse/human genome using 100 bp resolution

import sys, os
import anndata as ad
import crested
import numpy as np
import matplotlib
import gzip
import keras
import pysam

work_path = ""
mamm = "Mus_musculus"
model_id = ""

iter = int(sys.argv[1])

# load a trained model
model_path = f"{web_path}/CREsted_model/evolution_aware_model.keras"
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


######################################
### Step-3: merging prediction results

import gzip
import os, sys
import numpy as np
import pandas as pd

mamm = "Mus_musculus"
model_id = "XXX"

work_path = ""
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


