
#######################################
### Preparing the genomewide prediction

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


#########################
### Step-1: Download data

### Of note, here I downloaded genome references of 241 mammals from Zoonomia
### https://zoonomiaproject.org/

https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/mamm_list.txt


##############################################################################
### Step-2: Excluding chrs which are shorter than 2114 (at least one sequence)

import os, sys

origin_path = "XXX"
work_path = "XXX"

# Read all lines into a list
with open(os.path.join(work_path, "genome", "mamm_list.txt")) as f:
    mamm_list = [line.rstrip() for line in f if line.strip()]

# Get the mammal name based on index from argv
mamm = mamm_list[int(sys.argv[1]) - 1]
print(mamm)

input_file = os.path.join(origin_path, "genome_sizes", f"{mamm}.chrom.sizes")
output_file = os.path.join(work_path, "genome", mamm, f"{mamm}.chrom.sizes")

with open(input_file) as f_in, open(output_file, "w") as f_out:
    for line in f_in:
        chrom, size = line.rstrip().split("\t")
        if int(size) > 2200:
            f_out.write(line)


#################################################################
### Step-3: Creating the candidate sequence list for each species

<Python Script>

import sys, os
import gzip

work_path = "XXX"

# Read all lines into a list
with open(os.path.join(work_path, "genome", "mamm_list.txt")) as f:
    mamm_list = [line.rstrip() for line in f if line.strip()]

# Get the mammal name based on index from argv
mamm = mamm_list[int(sys.argv[1]) - 1]
print(mamm)

chrom_sizes = {}
with open(os.path.join(work_path, "genome", mamm, f"{mamm}.chrom.sizes")) as f:
    for line in f:
        chrom, size = line.rstrip().split("\t")
        chrom_sizes[chrom] = int(size)

window_size = 2114
step_size = 100
central_size = 1000
offset = (window_size - central_size) // 2

with gzip.open(os.path.join(work_path, "genome", mamm, "candidate_windows.txt.gz"), 'wt') as out:
    for chrom in chrom_sizes:
        print(f"Processing {chrom}")
        end_position = chrom_sizes[chrom]
        for pos in range(1, end_position - window_size + 1, step_size):
            out.write(f"{chrom}\t{pos}\n")

