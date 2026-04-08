
#######################################################
### Liftover human enhancer to mouse enhancer and vice, identifying human-only, mouse-only, and human-mouse-shared enhancers

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu



####################################
### Hum Obs / Mus Obs for three groups: 
### Mus only (lifted but not overlap with hum enhancer)
### Share
### Hum only (lifted but not overlap with mus enhancer)


################################################################
### Step-1: liftover is performed between mm10 and hg19 directly


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
options(scipen = 999)

celltype_list = c("Hepatocytes","Erythroid_cells")
mamm_list = c("Homo_sapiens_hg19", "Mus_musculus")

for(mamm in mamm_list){
    dat = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/prediction_mammals/prediction_", mamm, "_trf/window_list_10Mb.core_region.phred_score.bed"))
    for(celltype in celltype_list){
        dat_sub = dat[dat$V5 == celltype,]
        dat_sub = dat_sub[,c(1,2,3)]
        dat_sub$V4 = paste0(mamm, "_", celltype, "_", 1:nrow(dat_sub))
        write.table(dat_sub, paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/", mamm, "_", celltype, ".bed"), row.names=F, col.names=F, sep="\t", quote=F)
    }
}



work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap
script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install

celltype_list=(Hepatocytes Erythroid_cells)

for celltype in "${celltype_list[@]}"; do
    "$script_path"/liftOver -minMatch=0.1 \
        "$work_path"/Mus_musculus_"$celltype".bed \
        "$script_path"/liftOver_chain/mm10ToHg19.over.chain.gz \
        "$work_path"/liftover_Mus_musculus_Homo_sapiens_hg19_"$celltype".bed \
        "$work_path"/liftover_Mus_musculus_Homo_sapiens_hg19_"$celltype".unmapped.bed
done

for celltype in "${celltype_list[@]}"; do
    "$script_path"/liftOver -minMatch=0.1 \
        "$work_path"/Homo_sapiens_hg19_"$celltype".bed \
        "$script_path"/liftOver_chain/hg19ToMm10.over.chain.gz \
        "$work_path"/liftover_Homo_sapiens_hg19_Mus_musculus_"$celltype".bed \
        "$work_path"/liftover_Homo_sapiens_hg19_Mus_musculus_"$celltype".unmapped.bed
done

rm *.unmapped.bed

for celltype in "${celltype_list[@]}"; do
    bedtools intersect -a liftover_Mus_musculus_Homo_sapiens_hg19_"$celltype".bed \
    -b Homo_sapiens_hg19_"$celltype".bed \
    -wa | uniq > overlap_Mus_musculus_Homo_sapiens_hg19_"$celltype".bed
done

for celltype in "${celltype_list[@]}"; do
    bedtools intersect -a liftover_Homo_sapiens_hg19_Mus_musculus_"$celltype".bed \
    -b Mus_musculus_"$celltype".bed \
    -wa | uniq > overlap_Homo_sapiens_hg19_Mus_musculus_"$celltype".bed
done



############################
### Step-2: extract obs data


work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/human_mouse_enhancer
web_path=/net/shendure/vol10/www/content/members/cxqiu/public/nobackup/jax_atac_L2_mouse_raw_data/mm10
script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install

celltype=Hepatocytes

"$script_path"/bigWigAverageOverBed \
"$work_path"/human_fatal_data/liver_hepatoblasts.bw \
"$work_path"/human_mouse_enhancer_overlap/Homo_sapiens_hg19_"$celltype".bed \
"$work_path"/human_mouse_enhancer_overlap/Homo_sapiens_hg19_"$celltype".tab

"$script_path"/bigWigAverageOverBed \
"$work_path"/human_fatal_data/liver_hepatoblasts.bw \
"$work_path"/human_mouse_enhancer_overlap/liftover_Mus_musculus_Homo_sapiens_hg19_"$celltype".bed \
"$work_path"/human_mouse_enhancer_overlap/liftover_Mus_musculus_Homo_sapiens_hg19_"$celltype".tab

"$script_path"/bigWigAverageOverBed \
"$web_path"/Hepatocytes.bw \
"$work_path"/human_mouse_enhancer_overlap/Mus_musculus_"$celltype".bed \
"$work_path"/human_mouse_enhancer_overlap/Mus_musculus_"$celltype".tab

"$script_path"/bigWigAverageOverBed \
"$web_path"/Hepatocytes.bw \
"$work_path"/human_mouse_enhancer_overlap/liftover_Homo_sapiens_hg19_Mus_musculus_"$celltype".bed \
"$work_path"/human_mouse_enhancer_overlap/liftover_Homo_sapiens_hg19_Mus_musculus_"$celltype".tab


celltype=Erythroid_cells

"$script_path"/bigWigAverageOverBed \
"$work_path"/human_fatal_data/liver_erythroblasts.bw \
"$work_path"/human_mouse_enhancer_overlap/Homo_sapiens_hg19_"$celltype".bed \
"$work_path"/human_mouse_enhancer_overlap/Homo_sapiens_hg19_"$celltype".tab

"$script_path"/bigWigAverageOverBed \
"$work_path"/human_fatal_data/liver_erythroblasts.bw \
"$work_path"/human_mouse_enhancer_overlap/liftover_Mus_musculus_Homo_sapiens_hg19_"$celltype".bed \
"$work_path"/human_mouse_enhancer_overlap/liftover_Mus_musculus_Homo_sapiens_hg19_"$celltype".tab

"$script_path"/bigWigAverageOverBed \
"$web_path"/Definitive_erythroid.bw \
"$work_path"/human_mouse_enhancer_overlap/Mus_musculus_"$celltype".bed \
"$work_path"/human_mouse_enhancer_overlap/Mus_musculus_"$celltype".tab

"$script_path"/bigWigAverageOverBed \
"$web_path"/Definitive_erythroid.bw \
"$work_path"/human_mouse_enhancer_overlap/liftover_Homo_sapiens_hg19_Mus_musculus_"$celltype".bed \
"$work_path"/human_mouse_enhancer_overlap/liftover_Homo_sapiens_hg19_Mus_musculus_"$celltype".tab



############################
### Step-3: summaring result

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
options(scipen = 999)

celltype_list = c("Hepatocytes","Erythroid_cells")
mamm_list = c("Homo_sapiens_hg19", "Mus_musculus")

celltype = "Erythroid_cells"

dat_1 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/liftover_Homo_sapiens_hg19_Mus_musculus_", celltype, ".bed"))
dat_2 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/overlap_Homo_sapiens_hg19_Mus_musculus_", celltype, ".bed"))

human_lift = as.vector(dat_1$V4)
human_share = as.vector(dat_2$V4)
human_only = human_lift[!human_lift %in% human_share]

dat_1 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/liftover_Mus_musculus_Homo_sapiens_hg19_", celltype, ".bed"))
dat_2 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/overlap_Mus_musculus_Homo_sapiens_hg19_", celltype, ".bed"))

mouse_lift = as.vector(dat_1$V4)
mouse_share = as.vector(dat_2$V4)
mouse_only = mouse_lift[!mouse_lift %in% mouse_share]


dat_1 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/Homo_sapiens_hg19_", celltype, ".tab"))
dat_2 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/liftover_Homo_sapiens_hg19_Mus_musculus_", celltype, ".tab"))
mat_1 = dat_1[,c(1,5)] %>% left_join(dat_2[,c(1,5)], by = "V1") %>% filter(!is.na(V5.y))
colnames(mat_1) = c("enhancer_id", "human_obs", "mouse_obs")


dat_1 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/Mus_musculus_", celltype, ".tab"))
dat_2 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/liftover_Mus_musculus_Homo_sapiens_hg19_", celltype, ".tab"))
mat_2 = dat_2[,c(1,5)] %>% left_join(dat_1[,c(1,5)], by = "V1") %>% filter(!is.na(V5.y))
colnames(mat_2) = c("enhancer_id", "human_obs", "mouse_obs")

mat = rbind(mat_1, mat_2)
mat$category = "human_mouse_share"
mat$category[mat$enhancer_id %in% human_only] = "human_only"
mat$category[mat$enhancer_id %in% mouse_only] = "mouse_only"

mat$log2_human_mouse_obs_ratio = log2(mat$human_obs/mat$mouse_obs)
mat = mat[is.finite(mat$log2_human_mouse_obs_ratio),]
mat$category = factor(mat$category, levels = c("mouse_only", "human_mouse_share", "human_only"))

table(mat$category)

### Hepatocytes
     mouse_only human_mouse_share        human_only
             4948              2597              5790

### Erythroid_cells
       mouse_only human_mouse_share        human_only
             3492              1045              4353

p = ggplot(mat, aes(x = category, y = log2_human_mouse_obs_ratio)) +
  #geom_boxplot(fill = "#7aa457") + ## Hepatocytes
  geom_boxplot(fill = "#cb6751") +  ## Erythroid_cells
  labs(x="", y="Log2 (Hum Obs / Mus Obs)", title="") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave(paste0("~/share/", celltype, ".human_mouse_ratio.pdf"), p, width = 5, height = 5)


mat %>% group_by(category) %>% summarize(mean_log2_human_mouse_obs_ratio = mean(log2_human_mouse_obs_ratio))

### Hepatocytes
log2(a) - log2(b) = 1.88 - (-1.26)
then a/b = 2^(1.88 + 1.26) = 8.8

## Erythroid_cells
log2(a) - log2(b) = 1.68 - (-1.40)
then a/b = 2^(1.68 + 1.40) = 8.5




















