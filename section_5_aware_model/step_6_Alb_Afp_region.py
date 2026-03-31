### contribution score of candidate regions

### In Hepatocytes six candidate regions:
chr5 90414331 90414805
chr5 90431071 90431405
chr5 90439301 90439475
chr5 90447621 90447895
chr5 90485691 90485925
chr5 90490571 90490725


import anndata as ad
import crested
import numpy as np
import pandas as pd
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import keras
import pysam
import logomaker
from pathlib import Path
from scipy.stats import pearsonr

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
mamm = "Mus_musculus"
model_id = "mouse_fake_track_14"

model_path = f"{work_path}/mouse_fake_track_14/window_cluster/finetuned_model/checkpoints/06.keras"
model = keras.models.load_model(model_path, compile=False)
adata = ad.read_h5ad(os.path.join(work_path, "mouse_fake_track_14", "data_window_cluster_top3K.h5ad"))

genome = crested.Genome(
    "/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10/mm10.fa",
    "/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10/chromosome_sizes.txt"
)
crested.register_genome(genome)

def recenter_region(region):
    chrom, coords = region.split(":", 1)
    start, end = map(int, coords.split("-"))
    mid = (start + end) // 2
    return f"{chrom}:{int(mid-1057)}-{int(mid+1057)}"

regions_of_interest = [
    "chr5:90414331-90414805",
    "chr5:90431071-90431405",
    "chr5:90439301-90439475",
    "chr5:90447621-90447895",
    "chr5:90485691-90485925",
    "chr5:90490571-90490725",
]
classes_of_interest = ["Hepatocytes"]
class_idx = list(adata.obs_names.get_indexer(classes_of_interest))

scores, one_hot_encoded_sequences = crested.tl.contribution_scores(
    [recenter_region(i) for i in regions_of_interest],
    target_idx=class_idx,
    model=model,
)

# ── plot ──────────────────────────────────────────────────────────────
fig, axes = plt.subplots(len(regions_of_interest), 1, figsize=(20, 4 * len(regions_of_interest)))

for i, region in enumerate(regions_of_interest):
    chrom, coords = region.split(":", 1)
    start, end = map(int, coords.split("-"))
    mid = (start + end) // 2
    rec_start = mid - 1057
    offset_start = start - rec_start
    offset_end   = end   - rec_start
    contrib = scores[i, 0, offset_start:offset_end, :] * one_hot_encoded_sequences[i, offset_start:offset_end, :]
    contrib_df = pd.DataFrame(contrib, columns=["A", "C", "G", "T"])
    logomaker.Logo(contrib_df, ax=axes[i], color_scheme="classic")
    axes[i].set_title(region, fontsize=9, fontweight="bold")
    axes[i].set_ylabel("Contribution score")
    axes[i].axhline(0, color="black", linewidth=0.5, linestyle="--")

plt.tight_layout()
plt.savefig("/net/gs/vol1/home/cxqiu/share/contribution_scores_hepatocytes.pdf", dpi=300, bbox_inches="tight")
plt.show()







#### Within these regions, model-trimmed elements constitute only XX% (Alb & Afp) of bases, 
#### yet account for XX%, respectively, of nucleotides with phyloP > 2.

#########################################################################
### creating the candidate region for calculating the contribution scores

cd /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/conservation
script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install/

$script_path/bigWigToBedGraph mm10.60way.phyloP60wayPlacental.bw -chrom=chr5 mm10.60way.phyloP60wayPlacental.chr5.bedgraph

region_name=Afp_Alb
cd /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/conservation
bedtools intersect -a ../prediction_mammals/prediction_Mus_musculus_trf/window_list_10Mb.core_region.bed \
-b ./Afp_Alb/"$region_name"_region.bed \
-wa > ./Afp_Alb/"$region_name"_region.aware_core_region.bed



#################################
### calculate contribution scores

import anndata as ad
import crested
import numpy as np
import pandas as pd
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import keras
import pysam
import logomaker
from pathlib import Path
from scipy.stats import pearsonr

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
mamm = "Mus_musculus"
model_id = "mouse_fake_track_14"

model_path = f"{work_path}/mouse_fake_track_14/window_cluster/finetuned_model/checkpoints/06.keras"
model = keras.models.load_model(model_path, compile=False)
adata = ad.read_h5ad(os.path.join(work_path, "mouse_fake_track_14", "data_window_cluster_top3K.h5ad"))

genome = crested.Genome(
    "/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10/mm10.fa",
    "/net/shendure/vol10/projects/cxqiu/nobackup/genome/atac_data/mm10/chromosome_sizes.txt"
)
crested.register_genome(genome)

for region_name in ['Afp_Alb']:
    print(region_name)
    dfs = []
    cnt = 1
    with open(f"{work_path}/mouse_fake_track_14/conservation/Afp_Alb/{region_name}_region.aware_core_region.bed") as file:
        for line in file:
            print(cnt); cnt += 1
            l = line.rstrip().split('\t')
            chrom = l[0]; start = int(l[1]); end = int(l[2])
            mid = (start + end) // 2
            width = (end - start) // 2
            scores, one_hot_encoded_sequences = crested.tl.contribution_scores(
                f"{chrom}:{int(mid-1057)}-{int(mid+1057)}",
                target_idx=list(adata.obs_names.get_indexer([l[4]])),
                model=model,
            )
            result = np.multiply(scores.reshape(2114, 4), one_hot_encoded_sequences.reshape(2114, 4))[(1057 - width):(1057 + width)].sum(axis=1)
            df = pd.DataFrame({
                'chr': l[0],
                'start': np.arange(int(l[1]), int(l[2])),
                'end': np.arange(int(l[1]) + 1, int(l[2]) + 1),
                'score': result,
                'celltype': l[4],
                'enhancer': f"{l[0]}:{l[1]}-{l[2]}"
            })
            dfs.append(df)
    big_df = pd.concat(dfs, ignore_index=True)
    big_df.to_csv(f"{work_path}/mouse_fake_track_14/conservation/Afp_Alb/{region_name}_region.aware_core_region.contribution_score.bed", sep='\t', index=False, header=True)


##############################################################################################
### filtering bases with contribution scores > 0.1 (whatever cell type it has been associated)

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")
options(scipen = 999)

for(i in c('Afp_Alb')){
    dat = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Afp_Alb/", i, "_region.aware_core_region.contribution_score.bed"), header=T)
    dat_x = dat %>% group_by(chr, start, end) %>% slice_max(order_by = score, n = 1)
    dat_y = dat_x[dat_x$score > 0.1, c(1:3)]
    write.table(dat_y, paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Afp_Alb/", i, "_region.aware_core_region.filter.bed"), row.names=F, col.names=F, sep="\t", quote=F)
}

i=Afp_Alb
bedtools sort -i "$i"_region.aware_core_region.filter.bed | \
    bedtools merge > "$i"_region.aware_core_region.filter_x.bed
mv "$i"_region.aware_core_region.filter_x.bed "$i"_region.aware_core_region.filter.bed



##############################################################################
### overlapping with high contribution score regions and excluding CDS regions

cd /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/conservation
i=Afp_Alb

### how many base pairs in the whole region, and in the high contribution score regions, after excluding CDS
bedtools subtract \
    -a ./Afp_Alb/"${i}_region.bed" \
    -b mouse.v12.CDS.merged.bed \
    > ./Afp_Alb/"${i}_region.exclude_CDS.bed"

bedtools subtract \
    -a ./Afp_Alb/"$i"_region.aware_core_region.filter.bed \
    -b mouse.v12.CDS.merged.bed \
    > ./Afp_Alb/"$i"_region.aware_core_region.filter.exclude_CDS.bed

### how many base pairs having phylo score > 2 in the whole region, and in the high contribution score regions, after excluding CDS

bedtools intersect \
    -a mm10.60way.phyloP60wayPlacental.chr5.bedgraph \
    -b ./Afp_Alb/"${i}_region.exclude_CDS.bed" \
    > ./Afp_Alb/"${i}_region.exclude_CDS.bedgraph"

bedtools intersect \
    -a mm10.60way.phyloP60wayPlacental.chr5.bedgraph \
    -b ./Afp_Alb/"$i"_region.aware_core_region.filter.exclude_CDS.bed \
    > ./Afp_Alb/"$i"_region.aware_core_region.filter.exclude_CDS.bedgraph

i = "Afp_Alb"
a = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Afp_Alb/", i, "_region.exclude_CDS.bed"))
b = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Afp_Alb/", i, "_region.aware_core_region.filter.exclude_CDS.bed"))
print(sum(b$V3 - b$V2)/sum(a$V3 - a$V2))
### 1%

c = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Afp_Alb/", i, "_region.exclude_CDS.bedgraph"))
d = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Afp_Alb/", i, "_region.aware_core_region.filter.exclude_CDS.bedgraph"))
cutoff = 2
print(sum(d$V3[d$V4 > cutoff] - d$V2[d$V4 > cutoff])/sum(c$V3[c$V4 > cutoff] - c$V2[c$V4 > cutoff]))
### 8%










