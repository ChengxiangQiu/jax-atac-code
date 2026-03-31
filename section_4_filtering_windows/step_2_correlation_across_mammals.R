
###########################################################################
### Performing correlation between mouse windows and their liftover regions

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


###################################################################################################
### Step-1: filtering 100 bp windows which have been able to liftover to at least 120 other mammals

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

celltype_list = read.table("celltype_id.txt")
celltype_list = as.vector(celltype_list$V1)

mamm = "Mus_musculus"

peak_list_sub = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/window_list_merged.txt"))
colnames(peak_list_sub) = c("chr", "start", "end", "window_ID")
peak_list_sub$window_id = paste0(peak_list_sub$chr, "_", peak_list_sub$start, "_", peak_list_sub$end)

peak_list_sub_filter = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/peak_list_sub.txt"))
colnames(peak_list_sub_filter) = c("chr", "start", "end")
peak_list_sub_filter$window_id = paste0(peak_list_sub_filter$chr, "_", peak_list_sub_filter$start, "_", peak_list_sub_filter$end)

write.table(peak_list_sub[peak_list_sub$window_id %in% peak_list_sub_filter$window_id,],
            paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_", mamm, "/peak_list_filter.txt"), row.names=F, col.names=F, sep="\t", quote=F)


#################################################################
### Step-2: creating the list of candidate windows in each mammal

./corr_mammals/step_1_create_windows.py

########################
### Step-3: predicting each mammal

./corr_mammals/step_2_prediction.py

#####################################
### Step-4: merging prediction result

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

mamm_list = read.table(paste0(work_path, "/mamm_list.txt"))
mamm_list = as.vector(mamm_list$V1)
mamm_list = mamm_list[mamm_list != "Mus_musculus"]

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1]) ### 1:240
mamm = mamm_list[kk]

file_num = length(list.files(
    path = paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm),
    pattern = "\\.loc\\.txt\\.gz$"
))

dat = NULL
loc = NULL
for(i in 1:file_num){
    print(paste0(i, "/", file_num))
    dat_i = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/batch_", i,".txt.gz"))
    loc_i = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/batch_", i,".loc.txt.gz"))
    non_i = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/batch_", i,".nonzero.txt.gz"))
    dat_i = as.matrix(dat_i[,order(c(1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 36, 4, 5, 6, 7, 8, 9))])
    
    keep = non_i$V1 >= 2000
    dat = rbind(dat, dat_i[keep,])
    loc = rbind(loc, loc_i[keep,])
}

saveRDS(dat, paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/dat.rds"))
saveRDS(loc, paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/loc.rds"))


###########################################################################
### Step-5: checking correlation coefficient in a cell-type-specific manner

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

celltype_list = read.table("celltype_id.txt")
celltype_list = as.vector(celltype_list$V1)

dat = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/dat_sub.txt"))
peak_list = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/peak_list_sub.txt"))
colnames(peak_list) = c("chr", "start", "end")
peak_list$window_id = paste0(peak_list$chr, "_", peak_list$start, "_", peak_list$end)

peak_list_sub = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/peak_list_filter.txt"))
colnames(peak_list_sub) = c("chr", "start", "end", "window_ID", "window_id")

print(sum(peak_list_sub$window_id == peak_list$window_id))
peak_list = peak_list_sub
dat = as.matrix(dat)
rownames(dat) = as.vector(peak_list$window_ID)

mamm_list = read.table(paste0(work_path, "/14_crested/genome/mamm_list.txt"))
mamm_list = as.vector(mamm_list$V1)
mamm_list = mamm_list[mamm_list != "Mus_musculus"]

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1]) ### 1:240
mamm = mamm_list[kk]

dat_i = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/dat.rds"))
loc_i = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/loc.rds"))
colnames(loc_i) = c("chr", "start", "end", "window_ID")

library(data.table)
dt = as.data.table(dat_i)
dt[, group := as.vector(loc_i$window_ID)]
group_means = dt[, lapply(.SD, mean), by = group]
mat_result = as.matrix(group_means[, -"group"])
rownames(mat_result) = as.vector(group_means$group)

dat_sub = dat[rownames(mat_result),]

block_size = 10000
block_num = floor(nrow(dat_sub)/block_size) + 1

result = NULL
for(i in 1:block_num){
    print(i)
    keep = (block_size*(i-1) + 1):min((block_size*i), nrow(dat_sub))
    corr = cor(t(dat_sub[keep,]), t(mat_result[keep,]))
    result = c(result, diag(corr))
}

saveRDS(result, paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/result_corr.rds"))


##################################
### Step-6: summarizing the result

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

peak_list_sub = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/peak_list_filter.txt"))
colnames(peak_list_sub) = c("chr", "start", "end", "window_ID", "window_id")
### n = 1,465,915

mamm_list = read.table(paste0(work_path, "/14_crested/genome/mamm_list.txt"))
mamm_list = as.vector(mamm_list$V1)
mamm_list = mamm_list[mamm_list != "Mus_musculus"]

res = NULL
for(kk in 1:length(mamm_list)){
    mamm = mamm_list[kk]
    print(paste0(kk, "/", mamm))
    res_i = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/", mamm, "/result_corr.rds"))
    x = data.frame(window_ID = names(res_i),
                   corr = as.vector(res_i),
                   species = mamm)
    res = rbind(res, x)
}
saveRDS(res, paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/prediction/result_corr.rds"))

res_num = res %>% group_by(window_ID) %>% tally()
res = res[res$window_ID %in% as.vector(res_num$window_ID[res_num$n > 120]),]

### globally
res_x = res %>% group_by(window_ID) %>% summarize(mean_corr = mean(corr), median_corr = median(corr))
### n = 1,439,611

### total number of candidate windows for prediction
x = read.table('/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/genome/candidate_windows_sub_region_num.txt')
sum(x$V1)

p1 = ggplot() +
    geom_histogram(data = res_x, aes(x = mean_corr), bins = 100) +
    labs(x="Mean Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

p2 = ggplot() +
    geom_histogram(data = res_x, aes(x = median_corr), bins = 100) +
    labs(x="Median Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave("~/share/hist_corr_across_mammals_global.pdf", p1 + p2, width = 10, height = 5)

p2 = ggplot() +
    geom_histogram(data = res_x, aes(x = median_corr), bins = 50) +
    labs(x="Median Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    geom_vline(xintercept = 0.6) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 

ggsave("hist_corr_across_mammals_global_median_only.pdf", p2, width = 6, height = 4.5)


### separate by cell types
celltype_list = read.table("/net/shendure/vol2/projects/cxqiu/JAX_atac/Novaseq/call_peaks_celltype_L2/code/celltype_id.txt")
celltype_list = as.vector(celltype_list$V1)

peak_sig = NULL
for(i in 1:36){
    print(i)
    tmp = read.table(paste0(work_path, 
                            "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/celltype_", i, "_enhancer.bed"))
    tmp$V4 = celltype_list[i]
    colnames(tmp) = c("chr", "start", "end", "celltype")
    peak_sig = rbind(peak_sig, tmp)
}
peak_sig$window_id = paste0(peak_sig$chr, "_", peak_sig$start, "_", peak_sig$end)

peak_list_sub = peak_list_sub %>% left_join(peak_sig %>% select(window_id, celltype) %>% unique(), by = "window_id") %>% as.data.frame()

res_x = res_x %>% left_join(peak_list_sub[,c("window_ID", "celltype")], by = "window_ID")


p = ggplot() +
    geom_histogram(data = res_x, aes(x = mean_corr), bins = 100) +
    labs(x="Mean Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    facet_wrap(~ celltype, nrow = 6, ncol = 6) 
ggsave("~/share/hist_mean_corr_across_mammals_split_celltype.pdf", p, width = 15, height = 15)


p = ggplot() +
    geom_histogram(data = res_x, aes(x = median_corr), bins = 100) +
    labs(x="Median Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    facet_wrap(~ celltype, nrow = 6, ncol = 6) 
ggsave("~/share/hist_median_corr_across_mammals_split_celltype.pdf", p, width = 15, height = 15)


### check overlapping with promoters
bedtools intersect -a peak_list_sub.txt -b ../../../TSS_2500.bed -wa | uniq > peak_list_sub.overlap_promoter.txt

peak_list_sub_overlap_promoter = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/peak_list_sub.overlap_promoter.txt"))
colnames(peak_list_sub_overlap_promoter) = c("chr", "start", "end")
peak_list_sub_overlap_promoter$window_id = paste0(peak_list_sub_overlap_promoter$chr, "_",
                                                  peak_list_sub_overlap_promoter$start, "_",
                                                  peak_list_sub_overlap_promoter$end)
peak_list_sub$promoter = if_else(peak_list_sub$window_id %in% as.vector(peak_list_sub_overlap_promoter$window_id), "yes", "no")

res_x = res %>% group_by(window_ID) %>% summarize(mean_corr = mean(corr), median_corr = median(corr))
res_x = res_x %>% left_join(peak_list_sub[,c("window_ID", "promoter")], by = "window_ID")
res_x = unique(res_x)

p = ggplot() +
    geom_histogram(data = res_x, aes(x = mean_corr), bins = 100) +
    labs(x="Mean Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    facet_wrap(~ promoter, nrow = 1, ncol = 2) 
ggsave("~/share/hist_mean_corr_across_mammals_split_promoter.pdf", p, width = 10, height = 5)


p = ggplot() +
    geom_histogram(data = res_x, aes(x = median_corr), bins = 100) +
    labs(x="Median Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    facet_wrap(~ promoter, nrow = 1, ncol = 2) 
ggsave("~/share/hist_median_corr_across_mammals_split_promoter.pdf", p, width = 10, height = 5)

### subset 100 bp windows, n = 547,317
res_y = res_x[res_x$median_corr > 0.6,]
res_y = res_y %>% left_join(peak_list_sub[,c("window_ID","chr","start","end","celltype")], by = "window_ID")
res_y = res_y[,c("chr","start","end","promoter","celltype")]
write.table(res_y, "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_atac/14_crested/100bp_windows_pass_corr_0.6.txt", row.names=F, sep="\t", quote=F)


### check overlapping with exon
bedtools intersect -a peak_list_sub.txt -b /net/gs/vol1/home/cxqiu/work/tome/code/mouse.v12.exonID.txt -wa | uniq > peak_list_sub.overlap_exon.txt

peak_list_sub_overlap_exon = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/peak_list_sub.overlap_exon.txt"))
colnames(peak_list_sub_overlap_exon) = c("chr", "start", "end")
peak_list_sub_overlap_exon$window_id = paste0(peak_list_sub_overlap_exon$chr, "_",
                                              peak_list_sub_overlap_exon$start, "_",
                                              peak_list_sub_overlap_exon$end)
peak_list_sub$exon = if_else(peak_list_sub$window_id %in% as.vector(peak_list_sub_overlap_exon$window_id), "yes", "no")

res_x$promoter = NULL
res_x = res_x %>% left_join(peak_list_sub[,c("window_ID", "exon")], by = "window_ID")
res_x = unique(res_x)

p = ggplot() +
    geom_histogram(data = res_x, aes(x = mean_corr), bins = 100) +
    labs(x="Mean Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    facet_wrap(~ exon, nrow = 1, ncol = 2) 
ggsave("~/share/hist_mean_corr_across_mammals_split_exon.pdf", p, width = 10, height = 5)


p = ggplot() +
    geom_histogram(data = res_x, aes(x = median_corr), bins = 100) +
    labs(x="Median Corr Across Mammals", y="Counts", title="") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    facet_wrap(~ exon, nrow = 1, ncol = 2) 
ggsave("~/share/hist_median_corr_across_mammals_split_exon.pdf", p, width = 10, height = 5)

### subset 100 bp windows, n = 547,317
res_x$exon = NULL
res_y = res_x[res_x$median_corr > 0.6,]
res_y = res_y %>% left_join(peak_list_sub[,c("window_ID","chr","start","end","celltype","promoter","exon")], by = "window_ID")
res_y = res_y[,c("chr","start","end","median_corr","mean_corr","promoter","exon","celltype")]
res_y$median_corr = round(res_y$median_corr, 3)
res_y$mean_corr = round(res_y$mean_corr, 3)
write.table(res_y, "/net/shendure/vol10/www/content/members/cxqiu/private/nobackup/jax_atac/14_crested/100bp_windows_pass_corr_0.6.txt", row.names=F, sep="\t", quote=F)

### save those windows passed 0.6 for downstream analysis
### n = 547,317 (unique 100 bp windows)
res_y = res_x[res_x$median_corr > 0.6,]
res_y = res_y %>% left_join(peak_list_sub, by = "window_ID")
res_y = res_y[,c("chr","start","end","median_corr","mean_corr","promoter","exon","window_id","celltype")]
res_y$median_corr = round(res_y$median_corr, 3)
res_y$mean_corr = round(res_y$mean_corr, 3)
saveRDS(as.data.frame(res_y), paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/window_list.rds"))
saveRDS(as.data.frame(unique(res_y[,colnames(res_y) != 'celltype'])), paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/corr_Mus_musculus/window_list_uniq.rds"))





