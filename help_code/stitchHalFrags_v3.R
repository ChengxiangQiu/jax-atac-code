suppressMessages(suppressWarnings(library(GenomicFeatures)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(dplyr)))

args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("At least two argument must be supplied (original bed file) (hal lift bed file) (output bed file name) (minMatch) (max fraction).", call.=FALSE)
}


orig.bed <- args[1]
lift.bed <- args[2]
if(!file.exists(orig.bed)){
    stop(paste(orig.bed,"does not exist."))
}

if(!file.exists(lift.bed)){
    stop(paste(lift.bed,"does not exist."))
}

output.bed <- args[3]
minMatch <- as.numeric(args[4])
if (args[4] == "") {
  minMatch <- 0.5
}

max.frac <- as.numeric(args[5])
if (args[5] == "") {
  max.frac <- 1.5
}

#orig.bed = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/window_list_merged.txt"
#lift.bed = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_Mus_musculus/liftover/liftover_Acomys_cahirinus.bed"

orig.bed <- read.delim(orig.bed, sep ='\t', header = F)
orig.bed <- orig.bed %>% mutate(width = V3 - V2) %>% mutate(width_cutoff = width * minMatch, width_frac = width * max.frac)
lift.bed <- read.delim(lift.bed, sep = '\t', header = F)
filter.bed <- lift.bed %>% group_by(V4) %>% summarize(width = sum(V3 - V2)) %>%
    left_join(orig.bed[,c("V4", "width_cutoff")], by = "V4") %>% filter(width >= width_cutoff)
lift.bed <- lift.bed[lift.bed$V4 %in% as.vector(filter.bed$V4),]
print(paste("Filtered", nrow(orig.bed) - nrow(filter.bed), "interval(s) from original bed file with minMatch =", minMatch))

### consolidate fragments by chromosome and name
hal.join <- lift.bed %>% group_by(V1,V4) %>% summarise(start = min(V2), end = max(V3)) %>% mutate(width = end - start) 
hal.filter <- hal.join %>% group_by(V4) %>% summarise(width = sum(width))
hal.filter <- hal.filter %>% left_join(orig.bed[,c("V4", "width_frac")], by = "V4") %>% filter(width <= width_frac) 
hal.join <- hal.join[hal.join$V4 %in% as.vector(hal.filter$V4),]
print(paste("Filtered", nrow(filter.bed) - nrow(hal.join), "interval(s) from original bed file with and maxFrac =",max.frac))
print(paste("Total interval(s) filtered =", nrow(orig.bed) - nrow(hal.join)))

if(nrow(hal.join) > 0){
    colnames(hal.join) <- c('seqnames','name','start','end','width')
    hal.join <- hal.join[,c('seqnames','start','end','name')]
    hal.join <- GRanges(hal.join)    
    hal.join <- data.frame(hal.join)[,-4] ## remove width column
    hal.join <- data.frame(hal.join)[,-4] ## remove strand column
    write.table(hal.join, output.bed, sep = '\t', col.names = F, row.names = F, quote = F)
} else{
    write.table(hal.join, output.bed, sep = '\t', col.names = F, row.names = F, quote = F)
    print(paste("no regions recovered in",lift.bed))
}
