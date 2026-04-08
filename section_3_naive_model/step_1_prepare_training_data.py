
######################################################
### Preparing the data used for CREsted model training

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


#################################################################################
### Step-1: Generating the candidate peak list (focusing on cell class [Level-2])

bedtools intersect -a merged_peaks.bed -b TSS_2500.bed -v > tmp.bed
bedtools intersect -a tmp.bed -b mm10-blacklist.bed -v > merged_peaks_2L_non_overlap_with_promoter_2500_black_regions.bed
rm "$DTATPATH"/tmp.bed

### Level-2: n = 692058 peaks

### This file can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/merged_peaks_2L_non_overlap_with_promoter_2500_black_regions.bed


#######################################################################
### Step-2: Generating the bigwig profile for each cell class (level-2)

data_path=XXX
genome_path=XXX
WORKPATH=XXX

python \
./help_code/aggregate_cut_site_count.py \
"$data_path"/get_unique_fragments/celltype_"$SGE_TASK_ID".bed.gz \
"$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".bedgraph.tmp

cat "$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".bedgraph.tmp | \
sort -k1,1 -k2,2n > "$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".bedgraph

bedtools intersect \
-a "$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".bedgraph \
-b "$WORKPATH"/14_crested/TSS_2500.bed \
-v > "$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".no_promoter.bedgraph

python \
./help_code/normalize_bedgraph_by_CPM.py \
"$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".no_promoter.bedgraph \
"$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".no_promoter.norm_CPM.bedgraph

### bedGraphToBigWig is UCSC tools:
### https://hgdownload.gi.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig

bedGraphToBigWig \
"$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".no_promoter.norm_CPM.bedgraph \
"$genome_path"/chromosome_sizes.txt \
"$data_path"/BigWig_cut_site_norm/celltype_"$SGE_TASK_ID".bw












