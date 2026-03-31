
################################################################
### Identifying candidate top 100-bp windows for each cell class

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


######################################################################
### Step-1: After prediction, identifying top candidate 100-bp windows

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

mamm = "Mus_musculus"

dat = read.table(paste0(work_path, "/dat.txt.gz"))

peak_list = read.table(paste0(work_path, "/dat_loc.txt.gz"))
colnames(peak_list) = c("chr", "start")
peak_list$end = peak_list$start + 100

### reordering cell types
dat = as.matrix(dat[,order(c(1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 36, 4, 5, 6, 7, 8, 9))])

saveRDS(dat, paste0(work_path, "/dat.txt.rds"))
saveRDS(peak_list, paste0(work_path, "/dat_loc.txt.rds"))
write.table(peak_list, paste0(work_path, "/dat_loc.bed"), row.names=F, col.names=F, sep="\t", quote=F)

### using Z (predicted score) > 3 and Z (specificity) > 3 to nominate candidate 100 bp windows

dat_row_max = apply(dat, 1, max)
dat_exp = exp(dat - dat_row_max)
dat_x = dat_exp/apply(dat_exp,1,sum)

keep = !is.nan(dat_x[,1])
if(sum(!keep) != 0){
    dat = dat[keep,]
    dat_x = dat_x[keep,]
    peak_list = peak_list[keep,]
}

dat_scale = scale(dat)
dat_scale = unclass(dat_scale)
dat_scale = as.matrix(dat_scale, nrow(dat), ncol(dat))

dat_x_scale = scale(dat_x)
dat_x_scale = unclass(dat_x_scale)
dat_x_scale = as.matrix(dat_x_scale, nrow(dat_x), ncol(dat_x))

for(i in 1:ncol(dat_scale)){
    print(i)
    keep = dat_scale[,i] >= 3 & dat_x_scale[,i] >= 3
    
    peak_list_sub = data.frame(peak_list[keep,])
    if(nrow(peak_list_sub) > 0){
        peak_list_sub$ID = paste0("ID", 1:nrow(peak_list_sub))
        write.table(peak_list_sub, 
                    paste0(work_path, "/celltype_", i, "_enhancer.bed"), 
                    row.names=F, col.names=F, sep="\t", quote=F)
    }
}



####################################################
### Step-2: Creating a merged list of 100-bp windows

peak_sig = NULL
for(i in 1:36){
    print(i)
    tmp = read.table(paste0(work_path, 
                            "/celltype_", i, "_enhancer.bed"))
    tmp = tmp[,c(1:3)]
    colnames(tmp) = c("chr", "start", "end")
    peak_sig = rbind(peak_sig, tmp)
}

peak_sig_uniq = unique(peak_sig[,c(1:3)])
peak_sig_uniq$window_id = paste0(peak_sig_uniq$chr, "_", peak_sig_uniq$start, "_", peak_sig_uniq$end)
peak_list$window_id = paste0(peak_list$chr, "_", peak_list$start, "_", peak_list$end)
peak_list_sub = peak_list[peak_list$window_id %in% peak_sig_uniq$window_id,]
peak_list_sub$window_id = paste0("I", c(1:nrow(peak_list_sub)))

write.table(peak_list_sub, paste0(work_path, 
                              "/window_list_merged.txt"),
            row.names=F, col.names=F, sep="\t", quote=F)

# This list can be downloaded from:
# https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/Supplementary_File_1_100bp_Windows_Specific_to_36_Level2_Cell_Classes.bed.gz














