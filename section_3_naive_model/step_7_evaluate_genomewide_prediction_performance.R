

######################################################################################
### Evaluating the model prediction perfermance on either held-out set or through chr9

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu



###########################################################
### Step-1: Distribution of predicted score and specificity

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

mamm = "Mus_musculus"

### prediction outputs
dat = readRDS(paste0(work_path, "/dat.txt.rds"))
peak_list = readRDS(paste0(work_path, "/dat_loc.txt.rds"))

dat_row_max = apply(dat, 1, max)
dat_exp = exp(dat - dat_row_max)
dat_x = dat_exp/apply(dat_exp,1,sum)
rm(dat_row_max, dat_exp)

keep = sample(1:nrow(dat), size = round(nrow(dat) * 0.1))
dat = dat[keep,]
dat_x = dat_x[keep,]

celltype_list = read.table(paste0(web_path, "/cell_class_id.txt"))[,1]

df = data.frame(score = c(as.matrix(dat)),
                celltype = rep(celltype_list, each = nrow(dat)))
df$log2_score = log2(df$score + 1)

df_x = data.frame(score = c(as.matrix(dat_x)),
                celltype = rep(celltype_list, each = nrow(dat_x)))

p1 <- df %>%
    ggplot(aes(x = log2_score, color = celltype)) +
    geom_histogram(
        aes(y = after_stat(count + 1)),
        bins = 50,
        position = "identity",
        fill = NA,          # no fill
        linewidth = 0.8     # line thickness
    ) +
    scale_color_manual(values = celltype_L2_color_plate) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))
    ) +
    labs(
        x = "Log2 (Predicted Score + 1)",
        y = "Log2 (Count + 1)"
    ) +
    theme_classic(base_size = 15) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
    )

p2 = df %>%
    ggplot(aes(x = log2_score, color = celltype)) +
    geom_density(
        position = "identity",
        alpha = 0.6,
        adjust = 1
    ) +
    scale_color_manual(values = celltype_L2_color_plate) +
    labs(
        x = "Log2 (Predicted Score + 1)",
        y = "Density"
    ) +
    theme_classic(base_size = 15) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
    )
ggsave("Hist_Density_predicted_score.pdf", p1 + p2, width = 10, height = 5)


p1 <- df_x %>%
    ggplot(aes(x = score, color = celltype)) +
    geom_histogram(
        aes(y = after_stat(count + 1)),
        bins = 50,
        position = "identity",
        fill = NA,          # no fill
        linewidth = 0.8     # line thickness
    ) +
    scale_color_manual(values = celltype_L2_color_plate) +
    scale_y_continuous(
        trans = log2_trans(),
        breaks = trans_breaks("log2", function(x) 2^x),
        labels = trans_format("log2", math_format(2^.x))
    ) +
    labs(
        x = "Softmax (Predicted Score)",
        y = "Log2 (Count + 1)"
    ) +
    theme_classic(base_size = 15) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
    )

p2 = df_x %>%
    ggplot(aes(x = score, color = celltype)) +
    geom_density(
        position = "identity",
        alpha = 0.6,
        adjust = 1
    ) +
    scale_color_manual(values = celltype_L2_color_plate) +
    labs(
        x = "Softmax (Predicted Score)",
        y = "Density"
    ) +
    theme_classic(base_size = 15) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
    )
ggsave("Hist_Density_specificity.pdf", p1 + p2, width = 10, height = 5)



###################################################################################
### Step-2: Enrichment of the peaks used for training vs. the peaks from prediction

work_path=XXX
for i in $(seq 1 36); do
    echo $i
    sort -k1,1V -k2,2n -k3,3n "$work_path"/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/celltype_"$i"_enhancer.bed | bedtools merge -i - > "$work_path"/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/overlap_with_peaks_used_for_training/celltype_"$i"_enhancer.merged.bed
done

for i in $(seq 1 36); do
    echo $i
    for j in $(seq 1 36); do
        bedtools intersect -a "$work_path"/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/overlap_with_peaks_used_for_training/celltype_"$i"_enhancer.merged.bed \
        -b "$work_path"/14_crested/celltype_L2_cut_norm/peaks_used_for_training/celltype_"$j".10K.bed | \
        sort -k1,1V -k2,2n -k3,3n | \
        bedtools merge > "$work_path"/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/overlap_with_peaks_used_for_training/celltype_"$i"_"$j".10K.bed
    done
done

for i in $(seq 1 36); do
    echo $i
    for j in $(seq 1 36); do
        bedtools intersect -a "$work_path"/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/overlap_with_peaks_used_for_training/celltype_"$i"_enhancer.merged.bed \
        -b "$work_path"/14_crested/celltype_L2_cut_norm/peaks_used_for_training/celltype_"$j".3K.bed | \
        sort -k1,1V -k2,2n -k3,3n | \
        bedtools merge > "$work_path"/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/overlap_with_peaks_used_for_training/celltype_"$i"_"$j".3K.bed
    done
done

>>> summarize results (in R)

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

mamm = "Mus_musculus"

celltype_list = read.table(paste0(web_path, "/cell_class_id.txt"))

dat_background = NULL
for(i in 1:36){
    dat_x = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/peaks_used_for_training/celltype_", i, ".10K.bed"))
    dat_background = rbind(dat_background, data.frame(train = i, bp_train = sum(dat_x$V3 - dat_x$V2)))
}

dat = NULL
for(i in 1:36){
    print(i)
    for(j in 1:36){
        dat_x = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/overlap_with_peaks_used_for_training/celltype_", i, "_", j, ".10K.bed"))
        dat = rbind(dat, data.frame(pred = i, train = j, bp = sum(dat_x$V3 - dat_x$V2)))
    }
}

df = dat %>% left_join(dat_background, by = "train") %>% mutate(norm_bp = bp/bp_train)
colnames(celltype_list) = c("celltype_name_pred", "pred")
df = df %>% left_join(celltype_list, by = "pred")
colnames(celltype_list) = c("celltype_name_train", "train")
df = df %>% left_join(celltype_list, by = "train")

df_mat = dcast(df[,c("celltype_name_pred", "celltype_name_train", "norm_bp")], celltype_name_pred~celltype_name_train)
rownames(df_mat) = df_mat[,1]
df_mat = as.matrix(df_mat[,-1])

df_mat = df_mat[celltype_list[,1], celltype_list[,1]]

library("gplots")
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("Heatmap_pred_train_comparing.pdf", 8, 5)
heatmap.2(as.matrix(df_mat), 
          col=Colors, 
          #col=viridis(100), 
          scale="row", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(10,5))
dev.off()

celltype_order = c(13,15,26,23,35,34,21,6,36,30,16,20,4,25,19,11,24,17,29,22,32,12,31,28,33,27,14,7,5,3,2,9,18,10,8,1)
celltype_order = as.vector(celltype_list$celltype_name_train)[celltype_order]

df_mat = as.matrix(df_mat[celltype_order, celltype_order])

library("gplots")
library(RColorBrewer)
library(viridis)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(120)
pdf("Heatmap_pred_train_comparing_reorder.pdf", 8, 5)
heatmap.2(df_mat, 
          #col=Colors, 
          col=viridis(100), 
          scale="none", 
          Rowv = F, 
          Colv = F, 
          key=T, 
          density.info="none", 
          trace="none", 
          cexRow = 0.5, 
          cexCol = 0.5,
          margins = c(10,5))
dev.off()



##################################################################################
### Step-3: Comparing prediction and observation across all 100 bp windows on chr9

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

mamm = "Mus_musculus"

chr = "chr9"

dat = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_", mamm, "/dat_", chr, "_obs.txt.gz"))
peak_list = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_", mamm, "/dat_", chr, ".bed"))

### reordering cell types
dat_obs = as.matrix(dat[,order(c(1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 36, 4, 5, 6, 7, 8, 9))])

dat_pre = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_", mamm, "/dat.txt.rds"))
dat_loc = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_", mamm, "/dat_loc.txt.rds"))
keep = dat_loc$chr == chr
dat_pre = dat_pre[keep,]
dat_loc = dat_loc[keep,]

df = data.frame(pre = c(as.matrix(dat_pre)),
                obs = c(as.matrix(dat_obs)),
                chr = rep(dat_loc$chr, ncol(dat_pre)),
                start = rep(dat_loc$start, ncol(dat_pre)),
                end = rep(dat_loc$end, ncol(dat_pre)))
df$window_id = paste0(df$chr, "_", df$start, "_", df$end)
df$log_pre = log(df$pre + 1)
df$log_obs = log(df$obs + 1)
print(cor(df$log_pre, df$log_obs)) 

### chr9: 0.7082918

p = ggplot(data = df, aes(x=log_obs, y=log_pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, linetype = "dashed", color = "red") +
    geom_abline(intercept = -3, slope = 1, linetype = "dashed", color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title="") +
    #labs(x="Log (observation + 1)", y="Log (prediction + 1)", title=paste0("Corr = ", round(cor(df$log_obs, df$log_pre), 2))) +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("chr9_windows_naive.pdf"), p, height = 5, width = 5)


df_over_pred = df %>% filter(log_pre > log_obs + 2 | log_pre < log_obs - 3)
print(length(unique(df_over_pred$window_id))) ### 4566 windows
print(nrow(df_over_pred)/nrow(df)) ### 0.0002131177
print(length(unique(df_over_pred$window_id))/length(unique(df$window_id))) ### 0.003664759

df_over_pred = df %>% filter(log_pre > log_obs + 3)
print(nrow(unique(df_over_pred[,c("chr","start","end")])))
print(nrow(unique(df_over_pred[,c("chr","start","end")]))/nrow(dat_pre))
write.table(unique(df_over_pred[,c("chr","start","end")]), paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_Mus_musculus/", chr, "_window_list_over_pred.naive.txt"), row.names=F, col.names=F, sep="\t", quote=F)

df_under_pred = df %>% filter(log_pre < log_obs - 3)
print(nrow(unique(df_under_pred[,c("chr","start","end")])))
print(nrow(unique(df_under_pred[,c("chr","start","end")]))/nrow(dat_pre))

write.table(unique(df_under_pred[,c("chr","start","end")]), paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_Mus_musculus/", chr, "_window_list_over_pred.naive.txt"), row.names=F, col.names=F, sep="\t", quote=F)

df_over_pred = df %>% filter(log_pre > 3, log_obs > 3)
print(nrow(unique(df_over_pred[,c("chr","start","end")])))
print(nrow(unique(df_over_pred[,c("chr","start","end")]))/nrow(dat_pre))
write.table(unique(df_over_pred[,c("chr","start","end")]), paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_Mus_musculus/", chr, "_window_list_over_pred.naive.txt"), row.names=F, col.names=F, sep="\t", quote=F)


#################################################################
### Step-4: Color freq of windows overlapped with promoter or TRF

window_overlap_promoter = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_Mus_musculus/windows_overlap_with_TRF_promoters/chr9_windows_overlap_promoter.bed"))
colnames(window_overlap_promoter) = c("chr", "start", "end")
window_overlap_promoter$window_id = paste0(window_overlap_promoter$chr, "_", window_overlap_promoter$start, "_", window_overlap_promoter$end)

window_overlap_trf = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_Mus_musculus/windows_overlap_with_TRF_promoters/chr9_windows_overlap_trf.bed"))
colnames(window_overlap_trf) = c("chr", "start", "end")
window_overlap_trf$window_id = paste0(window_overlap_trf$chr, "_", window_overlap_trf$start, "_", window_overlap_trf$end)

df$overlap_promoter = df$window_id %in% window_overlap_promoter$window_id
df$overlap_trf = df$window_id %in% window_overlap_trf$window_id

df_long <- bind_rows(
  df %>% filter(overlap_promoter) %>% mutate(category = "promoter"),
  df %>% filter(overlap_trf) %>% mutate(category = "trf"),
  df %>% filter(!overlap_promoter & !overlap_trf) %>% mutate(category = "neither")
)

hex_df <- df_long %>%
  mutate(
    hx = cut(log_obs, 70),
    hy = cut(log_pre, 70)
  ) %>%
  group_by(hx, hy) %>%
  summarise(
    total = n(),
    promoter_n = sum(category == "promoter"),
    trf_n = sum(category == "trf"),
    neither_n = sum(category == "neither"),
    .groups = "drop"
  ) %>%
  mutate(
    hx_mid = (as.numeric(sub("\\((.+),(.+)\\]", "\\1", hx)) + as.numeric(sub("\\((.+),(.+)\\]", "\\2", hx))) / 2,
    hy_mid = (as.numeric(sub("\\((.+),(.+)\\]", "\\1", hy)) + as.numeric(sub("\\((.+),(.+)\\]", "\\2", hy))) / 2,
    score = (trf_n - promoter_n) / total
  )

p <- ggplot(hex_df, aes(x = hx_mid, y = hy_mid, color = score)) +
  geom_point(size = 1) +
  scale_color_gradientn(
    colors = Colors,
    limits = c(-1, 1),
    name = "TRF vs. Promoters"
  ) +
  geom_abline(intercept = 2, slope = 1, linetype = "dashed", color = "red") +
  geom_abline(intercept = -3, slope = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
  scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
  labs(
    x = "Log (observation + 1)",
    y = "Log (prediction + 1)",
    title = ""
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position="none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(color="black"),
    axis.text.y = element_text(color="black")
  )

ggsave(paste0("chr9_windows_naive_overlap_promoter_trf.pdf"), p, height = 5, width = 5)

p1 = ggplot(data = df[df$overlap_promoter,], aes(x=log_obs, y=log_pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, linetype = "dashed", color = "red") +
    geom_abline(intercept = -3, slope = 1, linetype = "dashed", color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title = "overlap with promoters") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
p2 = ggplot(data = df[df$overlap_trf,], aes(x=log_obs, y=log_pre)) +
    geom_hex(bins = 70) +
    geom_abline(intercept = 2, slope = 1, linetype = "dashed", color = "red") +
    geom_abline(intercept = -3, slope = 1, linetype = "dashed", color = "red") +
    scale_x_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    scale_y_continuous(limits = c(-0.3, 5.25), expand = expansion(mult = 0, add = 0)) +
    labs(x="Log (observation + 1)", y="Log (prediction + 1)", title = "overlap with TRFs") +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
ggsave(paste0("tmp.pdf"), p1 + p2, height = 5, width = 12)





######################################################################################################
### Step-5: investigate whether those over-predicted regions are enriched at some specific annotations

### merging windows

work_path=XXX
model_id=celltype_L2_cut_norm
model=naive
bedtools sort -i "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/chr9_window_list_over_pred."$model".txt \
| bedtools merge -i - > \
"$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/chr9_window_list_over_pred."$model".merged.txt


### using bedtools to perform intersect

work_path=XXX
model_id=celltype_L2_cut_norm
model=naive
repeats=(SINE Simple_repeat LINE LTR Low_complexity DNA Satellite)
for j in "${repeats[@]}"; do
    echo $j
    bedtools intersect -a "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/chr9_window_list_over_pred."$model".merged.txt \
    -b "$work_path"/5_sites_anno/sites_anno/RepeatMasker_subset/"$j".bed | \
    sort -k1,1V -k2,2n -k3,3n | \
    bedtools merge > "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/overlap_chr9_over_predicted_windows/dat_chr9_false_positive_"$j"_overlap."$model".bed
done

echo $model
bedtools intersect -a "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/chr9_window_list_over_pred."$model".merged.txt -b "$work_path"/5_sites_anno/sites_anno/promoter_filter.bed | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/overlap_chr9_over_predicted_windows/dat_chr9_false_positive_promoter_overlap."$model".bed

bedtools intersect -a "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/chr9_window_list_over_pred."$model".merged.txt -b "$work_path"/5_sites_anno/sites_anno/intergenic_filter.bed | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/overlap_chr9_over_predicted_windows/dat_chr9_false_positive_intergenic_overlap."$model".bed

bedtools intersect -a "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/chr9_window_list_over_pred."$model".merged.txt -b "$work_path"/14_crested/genome/Mus_musculus/mm10_SimpleRepeats.merged.bed | \
sort -k1,1V -k2,2n -k3,3n | \
bedtools merge > "$work_path"/14_crested/"$model_id"/prediction_mammals/prediction_Mus_musculus/overlap_chr9_over_predicted_windows/dat_chr9_false_positive_SimpleRepeats_overlap."$model".bed


### summarizing result using R

model_id = "celltype_L2_cut_norm"
model = "naive"
mamm = "Mus_musculus"
result = "both_high_chr9_windows"

dat = read.table(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "/", result, "/chr9_window_list_over_pred.", model, ".merged.txt"))
fpset_size = sum(dat$V3 - dat$V2)
repeat_list = c("SINE", "Simple_repeat", "LINE", "LTR", "Low_complexity", "DNA", "Satellite", "promoter", "intergenic", "SimpleRepeats")
chr9_size = 124595110
res = NULL
for(repeat_i in repeat_list){
    print(repeat_i)
    if (repeat_i %in% c("promoter", "intergenic")){
        dat_x = read.table(paste0(work_path, "/5_sites_anno/sites_anno/", repeat_i, "_filter.bed"))
    } else if (repeat_i %in% c("SimpleRepeats")){
        dat_x = read.table(paste0(work_path, "/14_crested/genome/Mus_musculus/mm10_", repeat_i, ".merged.bed"))
    } else {
        dat_x = read.table(paste0(work_path, "/5_sites_anno/sites_anno/RepeatMasker_subset/", repeat_i, ".bed"))
    }
    colnames(dat_x) = c("chr", "start", "end")
    dat_x = dat_x[dat_x$chr == "chr9",]
    gr <- GRanges(
        seqnames = dat_x$chr,
        ranges   = IRanges(start = dat_x$start, end = dat_x$end)
    )
    gr_merged <- data.frame(reduce(gr))
    category_size = sum(gr_merged$end - gr_merged$start)
    file_path = paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_", mamm, "/", result, "/overlap_chr9_over_predicted_windows/dat_chr9_false_positive_", repeat_i, "_overlap.", model, ".bed")
    if (file.exists(file_path) && file.info(file_path)$size > 0) {
      dat_y <- read.table(file_path)
      colnames(dat_y) = c("chr", "start", "end")
      overlap_size = sum(dat_y$end - dat_y$start)
    } else {
      overlap_size = 0
    }
    res = rbind(res, data.frame(category = repeat_i, overlap_size, fpset_size, category_size, chr9_size))
}

res$foldchange = (res$overlap_size/res$fpset_size)/(res$category_size/res$chr9_size)

res$category = factor(res$category, levels = rev(repeat_list))
p = ggplot(data=res, aes(x=category, y=foldchange)) +
    geom_bar(stat="identity") +
    coord_flip() +
    scale_y_continuous(limits = c(0, 13)) +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) 
    
ggsave(paste0("~/share/Corr_chr9_over_prediction_enrichment.", model, ".pdf"), p, width = 6, height = 5)





