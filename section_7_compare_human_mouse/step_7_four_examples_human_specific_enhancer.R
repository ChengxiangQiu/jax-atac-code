
##########################################################
#### find examples for visulization in UCSC genome browser

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu




gene │ hg38 │ hg19 | display_region
FECH │ chr18:57583221-57584055 │ chr18:55250453-55251287 | chr18:57,579,132-57,591,656
TFRC │ chr3:196107971-196108885 │ chr3:195834842-195835756 | chr3:196,074,873-196,129,772
CYP2C19 │ chr10:94747931-94748145 │ chr10:96507688-96507902 | chr10:94,731,245-94,769,944
APOB │ chr2:21092501-21092875 │ chr2:21315373-21315747 | chr2:21,042,571-21,095,222


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"

dat = data.frame(chr = c("chr18", "chr3", "chr10", "chr2"),
start = c(55250453, 195834842, 96507688, 21315373),
end = c(55251287, 195835756, 96507902, 21315747),
gene = c("FECH", "TFRC", "CYP2C19", "APOB"))

write.table(dat, paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/example/candidate.hg19.bed"), row.names=F, col.names=F, sep="\t", quote=F)


###########################################################
### Step-1: extract value from bigwig from silvia's dataset

work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/human_mouse_enhancer
script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install

"$script_path"/bigWigToBedGraph -chrom=chr18 -start=55200453 -end=55301287 \
"$work_path"/human_fatal_data/liver_hepatoblasts.bw \
"$work_path"/example/Hepatocytes.FECH.hg19.bedgraph

"$script_path"/bigWigToBedGraph -chrom=chr3 -start=195784842 -end=195885756 \
"$work_path"/human_fatal_data/liver_hepatoblasts.bw \
"$work_path"/example/Hepatocytes.TFRC.hg19.bedgraph

"$script_path"/bigWigToBedGraph -chrom=chr10 -start=96457688 -end=96557902 \
"$work_path"/human_fatal_data/liver_hepatoblasts.bw \
"$work_path"/example/Hepatocytes.CYP2C19.hg19.bedgraph

"$script_path"/bigWigToBedGraph -chrom=chr2 -start=21265373 -end=21365747 \
"$work_path"/human_fatal_data/liver_hepatoblasts.bw \
"$work_path"/example/Hepatocytes.APOB.hg19.bedgraph


"$script_path"/bigWigToBedGraph -chrom=chr18 -start=55200453 -end=55301287 \
"$work_path"/human_fatal_data/liver_erythroblasts.bw \
"$work_path"/example/Erythroid_cells.FECH.hg19.bedgraph

"$script_path"/bigWigToBedGraph -chrom=chr3 -start=195784842 -end=195885756 \
"$work_path"/human_fatal_data/liver_erythroblasts.bw \
"$work_path"/example/Erythroid_cells.TFRC.hg19.bedgraph

"$script_path"/bigWigToBedGraph -chrom=chr10 -start=96457688 -end=96557902 \
"$work_path"/human_fatal_data/liver_erythroblasts.bw \
"$work_path"/example/Erythroid_cells.CYP2C19.hg19.bedgraph

"$script_path"/bigWigToBedGraph -chrom=chr2 -start=21265373 -end=21365747 \
"$work_path"/human_fatal_data/liver_erythroblasts.bw \
"$work_path"/example/Erythroid_cells.APOB.hg19.bedgraph


######################################
### Step-2: liftover from hg19 to hg38

work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/human_mouse_enhancer
script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install

celltype_list=(Hepatocytes Erythroid_cells)
gene_list=(FECH TFRC CYP2C19 APOB)

for gene in "${gene_list[@]}"; do
    for celltype in "${celltype_list[@]}"; do
        "$script_path"/liftOver -minMatch=0.1 \
        "$work_path"/example/"${celltype}.${gene}.hg19.bedgraph" \
        "$script_path"/liftOver_chain/hg19ToHg38.over.chain.gz \
        "$work_path"/example/"${celltype}.${gene}.hg38.bedgraph" \
        "$work_path"/example/"${celltype}.${gene}.hg38.unmapped.bed"
    done
done

rm *unmapped.bed
rm *hg19.bedgraph

script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install
genome_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/genome/Homo_sapiens

celltype=Hepatocytes

cat "$celltype".APOB.hg38.bedgraph \
"$celltype".CYP2C19.hg38.bedgraph \
"$celltype".FECH.hg38.bedgraph \
"$celltype".TFRC.hg38.bedgraph \
> "$celltype".hg38.bedgraph

cat "$celltype".hg38.bedgraph | LC_COLLATE=C sort -k1,1 -k2,2n > "$celltype".hg38.sorted.bedgraph
"$script_path"/bedGraphToBigWig "$celltype".hg38.sorted.bedgraph "$genome_path"/Homo_sapiens.chrom.sizes.update "$celltype".bw
rm "$celltype".hg38.sorted.bedgraph


celltype=Erythroid_cells

cat "$celltype".APOB.hg38.bedgraph \
"$celltype".CYP2C19.hg38.bedgraph \
"$celltype".FECH.hg38.bedgraph \
"$celltype".TFRC.hg38.bedgraph \
> "$celltype".hg38.bedgraph

cat "$celltype".hg38.bedgraph | LC_COLLATE=C sort -k1,1 -k2,2n > "$celltype".hg38.sorted.bedgraph
"$script_path"/bedGraphToBigWig "$celltype".hg38.sorted.bedgraph "$genome_path"/Homo_sapiens.chrom.sizes.update "$celltype".bw
rm "$celltype".hg38.sorted.bedgraph

mv *bw /net/shendure/vol10/www/content/members/cxqiu/public/nobackup/human_fetal_atlas/hg38


https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/human_fetal_atlas/hub.txt


https://shendure-web.gs.washington.edu/content/members/cxqiu/public/nobackup/jax_atac_augmented_human_prediction/hub.txt






