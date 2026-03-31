###############################################################
### First, counting length or each chr

import sys, os

work_path = "XXX"

with open(os.path.join(work_path, "genome", "mamm_list.txt")) as f:
    mamm_list = [line.rstrip() for line in f if line.strip()]

# n = 240 mammals
mamm = mamm_list[int(sys.argv[1]) - 1]
print(mamm)

input_fasta = f"{work_path}/genome/{mamm}/{mamm}.fa"
output_file = f"{work_path}/genome/{mamm}/{mamm}.chrom.sizes.update"

lengths = {}
with open(input_fasta) as f:
    name = None
    seq_len = 0
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                lengths[name] = seq_len
            name = line[1:].split()[0]  # take first word after ">"
            seq_len = 0
        else:
            seq_len += len(line)
    if name is not None:
        lengths[name] = seq_len

with open(output_file, "w") as out:
    for name, seq_len in lengths.items():
        out.write(f"{name}\t{seq_len}\n")




###############################################################
### First, saving all the candidate windows (step size = 100bp)

import sys, os
import gzip

work_path = "XXX"

# Read all lines into a list
with open(os.path.join(work_path, "genome", "mamm_list.txt")) as f:
    mamm_list = [line.rstrip() for line in f if line.strip()]

mamm_list.remove('Mus_musculus')

# n = 240 mammals
mamm = mamm_list[int(sys.argv[1]) - 1]
print(mamm)

chrom_sizes = {}
with open(os.path.join(work_path, "genome", mamm, f"{mamm}.chrom.sizes.update")) as f:
    for line in f:
        chrom, size = line.rstrip().split("\t")
        chrom_sizes[chrom] = int(size)

peak_list = set()
with open(f"{work_path}/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/peak_list_filter.txt") as f:
    for line in f:
        peak_list.add(line.rstrip().split('\t')[3])

with open(f"{work_path}/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/liftover/stitchHalFrags_{mamm}.bed") as f, gzip.open(os.path.join(work_path, "genome", mamm, "candidate_windows_sub.txt.gz"), 'wt') as out:
    for line in f:
        l = line.rstrip().split("\t")
        if l[3] in peak_list and l[0] in chrom_sizes:
            center = round((int(l[1]) + int(l[2]))/2)
            if (center - 1057) > 0 and (center + 1057) < chrom_sizes[l[0]]:
                out.write(l[0] + '\t' + str(center-1057) + '\t' + str(center+1057) + '\t' + l[3] + '\n')

