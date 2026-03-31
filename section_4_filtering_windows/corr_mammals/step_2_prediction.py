
####################################################
### Second, predicting on multiple different mammals

import sys, os
import anndata as ad
import crested
import numpy as np
import matplotlib
import gzip
import keras
import pysam

work_path = "XXX"

chr_list = ["chr3", "chr5", "chr8", "chr10", "chr9", "chr18"]

# Get the mammal name based on index from argv
chr = chr_list[int(sys.argv[1]) - 1]
print(chr)

mamm = "Mus_musculus"

# load a trained model
model_path = os.path.join(work_path, "mouse_fake_track_3", "window_cluster", "basemodel", "checkpoints", "10.keras")
model = keras.models.load_model(model_path, compile=False)

output_path = f"{work_path}/mouse_fake_track_3/prediction_mammals/prediction_{mamm}/{chr}"
if not os.path.exists(output_path):
    os.makedirs(output_path)

fasta = pysam.FastaFile(os.path.join(work_path, "genome", mamm, f"{mamm}.fa"))

window_size = 2114
step_size = 10
central_size = 1000
offset = (window_size - central_size) // 2

all_regions = []
central_start = []

with gzip.open(os.path.join(work_path, "genome", mamm, f"candidate_windows_10bp_{chr}.txt.gz"), 'rt') as f:
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


