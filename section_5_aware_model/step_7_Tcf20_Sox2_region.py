
#######################################
### Focusing on Tcf20 and Sox2 regions 

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


#########################################################################
### Step-1: creating the candidate region for calculating the contribution scores

#### Within these regions, model-trimmed elements constitute only XX% (Tcf20) and XX% (Sox2) of bases, 
#### yet account for XX% and XX%, respectively, of nucleotides with phyloP > 2.

$script_path/bigWigToBedGraph mm10.60way.phyloP60wayPlacental.bw -chrom=chr15 mm10.60way.phyloP60wayPlacental.chr15.bedgraph
$script_path/bigWigToBedGraph mm10.60way.phyloP60wayPlacental.bw -chrom=chr3 mm10.60way.phyloP60wayPlacental.chr3.bedgraph

region_name=Tcf20
cd /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/conservation
bedtools intersect -a ../prediction_mammals/prediction_Mus_musculus_trf/window_list_10Mb.core_region.bed \
-b ./Tcf20_Sox2/"$region_name"_region.bed \
-wa > ./Tcf20_Sox2/"$region_name"_region.aware_core_region.bed

region_name=Sox2
cd /net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_14/conservation
bedtools intersect -a ../prediction_mammals/prediction_Mus_musculus_trf/window_list_10Mb.core_region.bed \
-b ./Tcf20_Sox2/"$region_name"_region.bed \
-wa > ./Tcf20_Sox2/"$region_name"_region.aware_core_region.bed


#########################################
### Step-2: calculate contribution scores

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

mamm = "Mus_musculus"

model_path = f"{web_path}/CREsted_model/evolution_aware_model.keras"
model = keras.models.load_model(model_path, compile=False)
adata = ad.read_h5ad(os.path.join(work_path, "data_window_cluster_top3K.h5ad"))

genome = crested.Genome(
    "mm10.fa",
    "chromosome_sizes.txt"
)
crested.register_genome(genome)

for region_name in ['Tcf20', 'Sox2']:
    print(region_name)
    dfs = []
    cnt = 1
    with open(f"{work_path}/mouse_fake_track_14/conservation/Tcf20_Sox2/{region_name}_region.aware_core_region.bed") as file:
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
    big_df.to_csv(f"{work_path}/mouse_fake_track_14/conservation/Tcf20_Sox2/{region_name}_region.aware_core_region.contribution_score.bed", sep='\t', index=False, header=True)


######################################################################################################
### Step-2: filtering bases with contribution scores > 0.1 (whatever cell type it has been associated)

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

for(i in c('Tcf20', 'Sox2')){
    dat = read.table(paste0(work_path, "/Tcf20_Sox2/", i, "_region.aware_core_region.contribution_score.bed"), header=T)
    dat_x = dat %>% group_by(chr, start, end) %>% slice_max(order_by = score, n = 1)
    dat_y = dat_x[dat_x$score > 0.1, c(1:3)]
    write.table(dat_y, paste0(work_path, "/Tcf20_Sox2/", i, "_region.aware_core_region.filter.bed"), row.names=F, col.names=F, sep="\t", quote=F)
}

i=Tcf20
bedtools sort -i "$i"_region.aware_core_region.filter.bed | \
    bedtools merge > "$i"_region.aware_core_region.filter_x.bed
mv "$i"_region.aware_core_region.filter_x.bed "$i"_region.aware_core_region.filter.bed

i=Sox2
bedtools sort -i "$i"_region.aware_core_region.filter.bed | \
    bedtools merge > "$i"_region.aware_core_region.filter_x.bed
mv "$i"_region.aware_core_region.filter_x.bed "$i"_region.aware_core_region.filter.bed


######################################################################################
### Step-3: overlapping with high contribution score regions and excluding CDS regions

i=Tcf20

### how many base pairs in the whole region, and in the high contribution score regions, after excluding CDS
bedtools subtract \
    -a ./Tcf20_Sox2/"${i}_region.bed" \
    -b mouse.v12.CDS.merged.bed \
    > ./Tcf20_Sox2/"${i}_region.exclude_CDS.bed"

bedtools subtract \
    -a ./Tcf20_Sox2/"$i"_region.aware_core_region.filter.bed \
    -b mouse.v12.CDS.merged.bed \
    > ./Tcf20_Sox2/"$i"_region.aware_core_region.filter.exclude_CDS.bed

### how many base pairs having phylo score > 2 in the whole region, and in the high contribution score regions, after excluding CDS

bedtools intersect \
    -a mm10.60way.phyloP60wayPlacental.chr15.bedgraph \
    -b ./Tcf20_Sox2/"${i}_region.exclude_CDS.bed" \
    > ./Tcf20_Sox2/"${i}_region.exclude_CDS.bedgraph"

bedtools intersect \
    -a mm10.60way.phyloP60wayPlacental.chr15.bedgraph \
    -b ./Tcf20_Sox2/"$i"_region.aware_core_region.filter.exclude_CDS.bed \
    > ./Tcf20_Sox2/"$i"_region.aware_core_region.filter.exclude_CDS.bedgraph

i = "Tcf20"
a = read.table(paste0(work_path, "/Tcf20_Sox2/", i, "_region.exclude_CDS.bed"))
b = read.table(paste0(work_path, "/Tcf20_Sox2/", i, "_region.aware_core_region.filter.exclude_CDS.bed"))
print(sum(b$V3 - b$V2)/sum(a$V3 - a$V2))
### 3%

c = read.table(paste0(work_path, "/Tcf20_Sox2/", i, "_region.exclude_CDS.bedgraph"))
d = read.table(paste0(work_path, "/Tcf20_Sox2/", i, "_region.aware_core_region.filter.exclude_CDS.bedgraph"))
cutoff = 2
print(sum(d$V3[d$V4 > cutoff] - d$V2[d$V4 > cutoff])/sum(c$V3[c$V4 > cutoff] - c$V2[c$V4 > cutoff]))
### 20%




i=Sox2

### how many base pairs in the whole region, and in the high contribution score regions, after excluding CDS
bedtools subtract \
    -a ./Tcf20_Sox2/"${i}_region.bed" \
    -b mouse.v12.CDS.merged.bed \
    > ./Tcf20_Sox2/"${i}_region.exclude_CDS.bed"

bedtools subtract \
    -a ./Tcf20_Sox2/"$i"_region.aware_core_region.filter.bed \
    -b mouse.v12.CDS.merged.bed \
    > ./Tcf20_Sox2/"$i"_region.aware_core_region.filter.exclude_CDS.bed

### how many base pairs having phylo score > 2 in the whole region, and in the high contribution score regions, after excluding CDS

bedtools intersect \
    -a mm10.60way.phyloP60wayPlacental.chr3.bedgraph \
    -b ./Tcf20_Sox2/"${i}_region.exclude_CDS.bed" \
    > ./Tcf20_Sox2/"${i}_region.exclude_CDS.bedgraph"

bedtools intersect \
    -a mm10.60way.phyloP60wayPlacental.chr3.bedgraph \
    -b ./Tcf20_Sox2/"$i"_region.aware_core_region.filter.exclude_CDS.bed \
    > ./Tcf20_Sox2/"$i"_region.aware_core_region.filter.exclude_CDS.bedgraph

i = "Sox2"
a = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Tcf20_Sox2/", i, "_region.exclude_CDS.bed"))
b = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Tcf20_Sox2/", i, "_region.aware_core_region.filter.exclude_CDS.bed"))
print(sum(b$V3 - b$V2)/sum(a$V3 - a$V2))
### 2.6%

c = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Tcf20_Sox2/", i, "_region.exclude_CDS.bedgraph"))
d = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/Tcf20_Sox2/", i, "_region.aware_core_region.filter.exclude_CDS.bedgraph"))
cutoff = 2
print(sum(d$V3[d$V4 > cutoff] - d$V2[d$V4 > cutoff])/sum(c$V3[c$V4 > cutoff] - c$V2[c$V4 > cutoff]))
### 6.6%









#####################################
### Step-4: A more general statement
### of the XX Mb of the genome with a PhyloP score of > 2, 
### XX Mb are accounted for by being protein-coding, while 
### XX Mb may be accounted for by a cell class-specific enhancer function nominated here

$script_path/bigWigToBedGraph mm10.60way.phyloP60wayPlacental.bw -chrom=chr3 mm10.60way.phyloP60wayPlacental.chr3.bedgraph

bedtools sort -i mouse.v12.CDS.merged.bed \
    | bedtools merge -i - \
    | bedtools intersect \
        -a mm10.60way.phyloP60wayPlacental.chr3.bedgraph \
        -b - \
        > mm10.60way.phyloP60wayPlacental.chr3.CDS.bedgraph

bedtools sort -i aware_core_region.bed \
    | bedtools merge -i - \
    | bedtools intersect \
        -a mm10.60way.phyloP60wayPlacental.chr3.bedgraph \
        -b - \
        > mm10.60way.phyloP60wayPlacental.chr3.aware_core_region.bedgraph

bedtools subtract \
    -a mm10.60way.phyloP60wayPlacental.chr3.aware_core_region.bedgraph \
    -b <(bedtools sort -i mouse.v12.CDS.merged.bed | bedtools merge -i -) \
    > mm10.60way.phyloP60wayPlacental.chr3.aware_core_region.excluded_CDS.bedgraph

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
source("~/work/scripts/utils.R")

a = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/mm10.60way.phyloP60wayPlacental.chr3.bedgraph"))
b = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/mm10.60way.phyloP60wayPlacental.chr3.CDS.bedgraph"))
c = read.table(paste0(work_path, "/14_crested/mouse_fake_track_14/conservation/mm10.60way.phyloP60wayPlacental.chr3.aware_core_region.excluded_CDS.bedgraph"))

for(cutoff in seq(0,1,0.1)){
    x = sum(a$V3[a$V4 > cutoff] - a$V2[a$V4 > cutoff])
    y = sum(b$V3[b$V4 > cutoff] - b$V2[b$V4 > cutoff])
    z = sum(c$V3[c$V4 > cutoff] - c$V2[c$V4 > cutoff])
    print(paste0("PhyloP score > ", cutoff, ": ", x, " / ", y, " / ", z))
}


bedtools sort -i mouse.v12.CDS.merged.bed \
    | bedtools merge -i - > tmp1

bedtools sort -i aware_core_region.bed \
    | bedtools merge -i - > tmp2

bedtools subtract \
    -a tmp1 \
    -b tmp2 \
    > tmp3




