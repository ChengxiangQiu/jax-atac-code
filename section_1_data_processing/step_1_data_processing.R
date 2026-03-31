
###################################

### The data generated in this study can be downloaded in raw and processed forms from:
### GSE325776
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE325776

### Additional data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


##########################################
### Step-1: processing the sequencing data

### The data generation was using sci-ATAC-seq3 protocol:
### https://www.protocols.io/view/sci-atac-seq3-ewov18xn7gr2/v1

### The data processing was following the previous pipeline:
### https://github.com/shendurelab/human-atac

### Note: This script is only used for reprocessing raw data and performing quality control.
### The processed and annotated dataset can be downloaded directly from the GEO (GSE325776)


####################################################
### Step-2: reading peak x cell matrix of 36 samples

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

blacklist_mm10 = readBed(paste0(web_path, "/mm10-blacklist.bed"))

### read sample name of 36 samples
sample_list_file = read.table(paste0(web_path, "/JAX_ATAC_sample_list.txt"))
experiment_id_list = sample_list_file[,1] ### experiment label
sample_id_list = sample_list_file[,2]     ### sample ID label
sample_list = sample_list_file[,3]        ### time point label

for(kk in 1:length(sample_list)){
    sample_i = sample_list[kk]
    counts = Read10X(paste0(work_path, sample_i), gene.column = 1) ### the peak x cell data could be downloaded from GEO site (GSE325776; peak_matrix.mtx, peak_matrix.row.txt, peak_matrix.column.txt)
    metadata = read.table(paste0(web_path, "/sample_sheet/", sample_id_list[kk], ".sample_info.txt"), header=T, row.names=2, as.is=T)
    
    if(sum(colnames(counts)!=rownames(metadata))==0){
        print("no mistake on barcode and cell name matching.")
        chrom_assay = CreateChromatinAssay(
            counts = counts,
            sep = c("_","_"),
            genome = "mm10",
            fragments = paste0(work_path, sample_i, ".fragments.txt.gz"), ### the fragment files could be downloaded from GEO site (GSE325776)
            min.cells = 0,
            min.features = 0
        )
        pbmc = CreateSeuratObject(counts=chrom_assay, assay='peaks', meta.data=metadata)
    }
    
    annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevelsStyle(annotations) <- 'UCSC'
    genome(annotations) <- "mm10"
    Annotation(pbmc) = annotations

    pbmc <- NucleosomeSignal(object = pbmc)
    pbmc$pct_reads_in_peaks <- pbmc$reads_peak / pbmc$reads * 100
    pbmc$pct_reads_in_tss <- pbmc$reads_tss / pbmc$reads * 100
    pbmc$blacklist_ratio <- FractionCountsInRegion(
        object = pbmc,
        assay = 'peaks',
        regions = blacklist_mm10
    )

    p1 <- VlnPlot(pbmc, c('pct_reads_in_peaks', 'pct_reads_in_tss', 'reads_peak'), pt.size = 0.1)
    p2 <- VlnPlot(pbmc, c('blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1) & scale_y_log10()
    ggsave(paste0(work_path, "/", sample_i, ".QC_1.pdf"), p1 | p2, width = 13, height = 4)

    pd = data.frame(pbmc[[]])

    p1 = ggplot(pd, aes(pct_reads_in_peaks)) + geom_histogram(aes(y=..density..), binwidth=0.25) + geom_density(alpha=.2) + theme_classic(base_size = 5) + scale_x_continuous(breaks = seq(1,round(max(pd$pct_reads_in_peaks)), by = 1))
    p2 = ggplot(pd, aes(pct_reads_in_peaks)) + geom_histogram(binwidth=0.25) + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) + theme_classic(base_size = 5)  + scale_x_continuous(breaks = seq(1,round(max(pd$pct_reads_in_peaks)), by = 1))
    p3 = ggplot(pd, aes(pct_reads_in_tss)) + geom_histogram(aes(y=..density..), binwidth=0.25) + geom_density(alpha=.2) + theme_classic(base_size = 5) + scale_x_continuous(breaks = seq(1,round(max(pd$pct_reads_in_tss)), by = 1))
    p4 = ggplot(pd, aes(pct_reads_in_tss)) + geom_histogram(binwidth=0.25) + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) + theme_classic(base_size = 5)  + scale_x_continuous(breaks = seq(1,round(max(pd$pct_reads_in_peaks)), by = 1))
    ggsave(paste0(work_path, "/", sample_i, ".QC_2.pdf"), (p1 | p2) / (p3 | p4), width = 12, height = 6)

    saveRDS(pd, paste0(work_path, "/", sample_i, ".pre_filter_pd.rds"))
}


#################################################################################
### Step-3: filtering cells based on quality distributions for individual samples

### Unique reads: >1,000
### Unique reads in peaks: 200 - 200,000
### Fraction of reads in peaks: above a sample-specific cutoff (15 - 37%)
### Fraction of reads in TSSs (+/-1 kb): above a sample-specific cutoff (8 - 13%)
### Fraction of reads in ENCODE blacklist regions: 0.001 - 0.1
### Nucleosome signal (calculated using the NucleosomeSignal function in Signac): <10

### read sample name of 36 samples
sample_list_file = read.table(paste0(web_path, "/JAX_ATAC_sample_list.txt"))
experiment_id_list = sample_list_file[,1] ### experiment label
sample_id_list = sample_list_file[,2]     ### sample ID label
sample_list = sample_list_file[,3]        ### time point label

pct_reads_in_peaks_cutoff = c("JAX_E10.0-L3-01"  = 26,
                              "JAX_E10.5_L4-07"  = 25,
                              "JAX_E10.75_L4-04" = 26,
                              "JAX_E10.75_L5-03" = 24,
                              "JAX_E11.0_L4-02"  = 26,
                              "JAX_E11.5_L4-02"  = 18,
                              "JAX_E11.75_L4-06" = 27,
                              "JAX_E11.75_L5-03" = 26,
                              "JAX_E12.25_L2-03" = 28,
                              "JAX_E12.5-L3-05"  = 28,
                              
                              "JAX_E12.75_L7-02" = 26,
                              "JAX_E12.75_L7-08" = 27,
                              "JAX_E13.0_L1-02"  = 25,
                              "JAX_E13.5_L11_07" = 34,
                              "JAX_E13.5_L13_02" = 27,
                              "JAX_E13.5_L3-06"  = 27,
                              "JAX_E13.5_L5_01"  = 28,
                              "JAX_E14.5_L4-05"  = 35,
                              "JAX_E14.5_L5-07"  = 34,
                              "JAX_E14.5_L6B-09" = 34,
                              
                              "JAX_E14.5_L7-07"  = 37,
                              "JAX_E15.25_L1-03" = 28,
                              "JAX_E15.5_L2-02"  = 30,
                              "JAX_E15.75_L1-02" = 30,
                              "JAX_E16.0_L2-03"  = 25,
                              "JAX_E16.25_L1-03" = 23,
                              "JAX_E16.5_L1-07"  = 23,
                              "JAX_E16.75_L1-03" = 23,
                              "JAX_E17.0_L1-02"  = 24,
                              "JAX_E17.25_L1-05" = 23,
                              
                              "JAX_E17.5_L1-05"  = 23,
                              "JAX_E18.0_L1-02"  = 23,
                              "JAX_E18.25_L1-05" = 22,
                              "JAX_E18.50_L1-07" = 23,
                              "JAX_E18.75_L1-07" = 22,
                              "JAX_P0_L1-08"     = 15)

pct_reads_in_tss_cutoff = c("JAX_E10.0-L3-01"  = 9,
                            "JAX_E10.5_L4-07"  = 9,
                            "JAX_E10.75_L4-04" = 9,
                            "JAX_E10.75_L5-03" = 9,
                            "JAX_E11.0_L4-02"  = 9,
                            "JAX_E11.5_L4-02"  = 8,
                            "JAX_E11.75_L4-06" = 9,
                            "JAX_E11.75_L5-03" = 9,
                            "JAX_E12.25_L2-03" = 10,
                            "JAX_E12.5-L3-05"  = 9,
                            
                            "JAX_E12.75_L7-02" = 9,
                            "JAX_E12.75_L7-08" = 9,
                            "JAX_E13.0_L1-02"  = 9,
                            "JAX_E13.5_L11_07" = 11,
                            "JAX_E13.5_L13_02" = 9,
                            "JAX_E13.5_L3-06"  = 9,
                            "JAX_E13.5_L5_01"  = 9,
                            "JAX_E14.5_L4-05"  = 13,
                            "JAX_E14.5_L5-07"  = 12,
                            "JAX_E14.5_L6B-09" = 12,
                            
                            "JAX_E14.5_L7-07"  = 12,
                            "JAX_E15.25_L1-03" = 10,
                            "JAX_E15.5_L2-02"  = 10,
                            "JAX_E15.75_L1-02" = 10,
                            "JAX_E16.0_L2-03"  = 8,
                            "JAX_E16.25_L1-03" = 8,
                            "JAX_E16.5_L1-07"  = 8,
                            "JAX_E16.75_L1-03" = 8,
                            "JAX_E17.0_L1-02"  = 8,
                            "JAX_E17.25_L1-05" = 8,
                            
                            "JAX_E17.5_L1-05"  = 8,
                            "JAX_E18.0_L1-02"  = 8,
                            "JAX_E18.25_L1-05" = 8,
                            "JAX_E18.50_L1-07" = 8,
                            "JAX_E18.75_L1-07" = 8,
                            "JAX_P0_L1-08"     = 8)

for(kk in 1:length(sample_list)){
    sample_i = sample_list[kk]
    pd = readRDS(paste0(work_path, "/", sample_i, ".pre_filter_pd.rds"))
    
    keep = pd$reads > 1000 &
        pd$reads_peak > 200  &
        pd$reads_peak < 2e5  &
        pd$pct_reads_in_peaks > pct_reads_in_peaks_cutoff[sample_id_list[kk]] &
        pd$pct_reads_in_tss > pct_reads_in_tss_cutoff[sample_id_list[kk]] &
        pd$blacklist_ratio < 0.1 &
        pd$blacklist_ratio > 1e-3 &
        pd$nucleosome_signal < 10
    print(paste0(sample_i, " : ", pct_reads_in_peaks_cutoff[sample_id_list[kk]], " : ", pct_reads_in_tss_cutoff[sample_id_list[kk]], " : ", sum(keep), " : ", nrow(pd), " : ", round(sum(keep)/nrow(pd),3)))
    pd$keep = keep
    saveRDS(pd, paste0(work_path, "/", sample_i, ".pre_filter_pd.rds"))
}


##########################################################
### Step-4: merging 36 samples and creating BPCells object

### BPCells: https://github.com/bnprks/BPCells

pd_all = NULL
for(kk in 1:length(sample_list)){
    sample_i = sample_list[kk]; print(sample_i)
    pd = readRDS(paste0(work_path, "/", sample_i, ".pre_filter_pd.rds"))
    pd_all = rbind(pd_all, pd)
}
print(table(pd_all$keep))

pd_all$barcode = unlist(lapply(as.vector(pd_all$id), function(x) strsplit(x,"[_]")[[1]][2])) 
pd_all$sample_id = unlist(lapply(as.vector(pd_all$id), function(x) strsplit(x,"[_]")[[1]][1])) 
rownames(pd_all) = NULL
pd_all$cell_id = as.vector(pd_all$id)
pd_all$id = NULL

sex = rep("Female", nrow(pd_all))
sex[pd_all$sex == "M"] = "Male"
pd_all$sex = as.vector(sex)

pd_all_uniq = unique(pd_all[,c("sample","sample_id")])
pd_all_uniq$sample_name = pd_all_uniq$sample

### Reading the doublet detection result from AMULET
### AMULET: https://github.com/UcarLab/AMULET

amulet_result = readRDS(paste0(web_path, "/AMULET_res.rds"))
amulet_result = amulet_result %>%
    left_join(pd_all_uniq[,c("sample_id","sample_name")], by = "sample_name")
amulet_result$cell_id = paste0(amulet_result$sample_id, "_", amulet_result$cell_id)
amulet_result$amulet_pval = amulet_result$p_val
amulet_result$amulet_qval = amulet_result$q_val

pd_all = pd_all %>%
    left_join(amulet_result[,c("cell_id", "amulet_pval", "amulet_qval")], by = "cell_id")
pd_all$orig.ident = NULL
saveRDS(pd_all, paste0(work_path, "/pd_all.rds"))

rownames(pd_all) = as.vector(pd_all$cell_id)
pd_all = pd_all[pd_all$keep,]
print(dim(pd_all))
### n = 4,344,905 cells, this is before filtering doublets

pd_x = unique(pd_all[,c("sample_id","sample")])
sample_list = as.vector(pd_x$sample)

for(sample_i in sample_list){
    x = pd_x$sample_id[pd_x$sample == sample_i]

    counts = Read10X(paste0(work_path, sample_i), gene.column = 1) ### the peak x cell data could be downloaded from GEO site (GSE325776; peak_matrix.mtx, peak_matrix.row.txt, peak_matrix.column.txt)
    colnames(counts) = paste0(x, "_", colnames(counts))

    counts_sub = counts[,colnames(counts) %in% as.vector(pd_all$cell_id)]
    counts_sub@x[counts_sub@x > 0] = 1
    
    pd_sub = pd_all[colnames(counts_sub),]
    print(dim(pd_sub))
    
    Matrix::writeMM(t(counts_sub), paste0(work_path, "/gene_count_", x,".mtx"))
    write.csv(pd_sub, paste0(work_path, "/df_cell_", x, ".csv"))
    write.csv(rownames(counts_sub), paste0(work_path, "/df_gene_", x, ".csv"), row.names = F)

}

### Using scanpy to merge samples and create a combined h5ad file, followed by reading in BPCells
### python help_code/create_scanpy_object.py













