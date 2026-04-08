
########################################################################################
### Enhancers which are exclusively in human (even not have lifted over region in mouse)

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


####################################################
### Step-1: Enhancers which are exclusively in human (even not have lifted over region in mouse)


work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/human_mouse_enhancer
web_path=/net/shendure/vol10/www/content/members/cxqiu/public/nobackup/jax_atac_L2_mouse_raw_data/mm10
script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install

cat "$work_path"/human_mouse_enhancer_overlap/Homo_sapiens_hg19_Erythroid_cells.bed \
"$work_path"/human_mouse_enhancer_overlap/Homo_sapiens_hg19_Hepatocytes.bed \
> "$work_path"/human_only/Homo_sapiens_hg19.bed

"$script_path"/bigWigAverageOverBed \
"$work_path"/human_fatal_data/liver_hepatoblasts.bw \
"$work_path"/human_only/Homo_sapiens_hg19.bed \
"$work_path"/human_only/Homo_sapiens_hg19.Hepatocytes.tab

"$script_path"/bigWigAverageOverBed \
"$work_path"/human_fatal_data/liver_erythroblasts.bw \
"$work_path"/human_only/Homo_sapiens_hg19.bed \
"$work_path"/human_only/Homo_sapiens_hg19.Erythroid_cells.tab







###############################
### Step-2: Plotting the result

model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
options(scipen = 999)

celltype_list = c("Hepatocytes","Erythroid_cells")
mamm_list = c("Homo_sapiens_hg19", "Mus_musculus")

dat_1 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/liftover_Homo_sapiens_hg19_Mus_musculus_Hepatocytes.bed"))
dat_2 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_overlap/liftover_Homo_sapiens_hg19_Mus_musculus_Erythroid_cells.bed"))

mat_1 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_only/Homo_sapiens_hg19.Hepatocytes.tab"))
mat_2 = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_only/Homo_sapiens_hg19.Erythroid_cells.tab"))
sum(mat_1$V1 != mat_2$V1)

mat = data.frame(enhancer_id = mat_1[,1], Hepatocytes = mat_1[,5], Erythroid_cells = mat_2[,5])

mat = mat[!mat$enhancer_id %in% c(dat_1$V4, dat_2$V4),]

mat$log2_Hepatocytes = log2(mat$Hepatocytes)
mat$log2_Erythroid_cells = log2(mat$Erythroid_cells)

mat$celltype = sub("_[0-9]+$", "", mat$enhancer_id)
mat$celltype = sub("Homo_sapiens_hg19_", "", mat$celltype)

mat = mat[is.finite(mat$log2_Hepatocytes) & is.finite(mat$log2_Erythroid_cells),]

print(table(mat$celltype))

Erythroid_cells     Hepatocytes
           3842            3239


p = ggplot(mat[sample(1:nrow(mat)),], aes(x = log2_Erythroid_cells, y = log2_Hepatocytes, color = celltype)) +
  geom_point(alpha = 0.6) +
  labs(x="Log2(Erythroid Obs)", y="Log2(Hepatocyte Obs)", title="") +
  theme_classic(base_size = 10) +
  scale_color_manual(values=c("Hepatocytes" = "#7aa457", "Erythroid_cells" = "#cb6751")) +
  coord_cartesian(xlim = c(-15, -2), ylim = c(-15, -2)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave(paste0("~/share/human_only_hepatocyte_erythroid_ratio.pdf"), p, width = 5, height = 5)


mat$log2_fc = log2(mat$Hepatocytes/mat$Erythroid_cells)
mat$celltype = factor(mat$celltype, levels = c("Hepatocytes", "Erythroid_cells"))

p = ggplot(mat, aes(x = celltype, y = log2_fc, fill = celltype)) +
geom_boxplot() +
labs(x="", y="Log2(Hepatocytes/Erythroid_cells)", title="") +
scale_fill_manual(values=c("Hepatocytes" = "#7aa457", "Erythroid_cells" = "#cb6751")) +
theme_classic(base_size = 10) +
theme(plot.title = element_text(hjust = 0.5)) +
theme(legend.position="none") +
theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave(paste0("~/share/human_only_hepatocyte_erythroid_ratio_boxplot.pdf"), p, width = 5, height = 5)


mat %>% group_by(celltype) %>% summarize(mean(log2_fc))

2^(2.13 + 0.699)



