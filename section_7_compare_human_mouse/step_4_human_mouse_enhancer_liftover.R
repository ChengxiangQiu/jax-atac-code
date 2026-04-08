
#######################################################
### Liftover human enhancer to mouse enhancer and vice, identifying which enhancers are shared, lifted-over but not overlapped, not even lifted over

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


############################
### Step-1: perform liftover

model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
options(scipen = 999)

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")
mamm_list = c("Homo_sapiens", "Mus_musculus")

for(mamm in mamm_list){
    dat = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/prediction_mammals/prediction_", mamm, "_trf/window_list_10Mb.core_region.phred_score.bed"))
    for(celltype in celltype_list){
        dat_sub = dat[dat$V5 == celltype,]
        dat_sub = dat_sub[,c(1,2,3)]
        dat_sub$V4 = paste0(mamm, "_", celltype, "_", 1:nrow(dat_sub))
        write.table(dat_sub, paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_liftover/", mamm, "_", celltype, ".bed"), row.names=F, col.names=F, sep="\t", quote=F)
    }
}


########################################
### Step-2: liftover from mouse to human

script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install/cactus-bin-v2.9.9/bin
data_path=/net/shendure/vol8/projects/tli/ucsc_cactus/
work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/human_mouse_enhancer

mamm=Mus_musculus
target_mamm=Homo_sapiens

celltype_list=(Adipocyte_cells Adipocyte_cells_Cyp2e1 B_cells Brain_capillary_endothelial_cells CNS_neurons Cardiomyocytes Corticofugal_neurons Endocardial_cells Endothelium Epithelial_cells Erythroid_cells Eye Glia Glomerular_endothelial_cells Gut_epithelial_cells Hepatocytes Intermediate_neuronal_progenitors Kidney Lateral_plate_and_intermediate_mesoderm Liver_sinusoidal_endothelial_cells Lung_and_airway Lymphatic_vessel_endothelial_cells Melanocyte_cells Mesoderm Neural_crest_PNS_neurons Neuroectoderm_and_glia Olfactory_ensheathing_cells Olfactory_neurons Oligodendrocytes Skeletal_muscle_cells T_cells White_blood_cells)
celltype="${celltype_list[$SGE_TASK_ID - 1]}"
echo $celltype

"$script_path"/halLiftover --noDupes --bedType 4 \
"$data_path/241-mammalian-2020v2.hal" \
"${mamm}" "$work_path/human_mouse_enhancer_liftover/${mamm}_${celltype}.bed" \
"${target_mamm}" "$work_path/human_mouse_enhancer_liftover/liftover_${mamm}_${celltype}_to_${target_mamm}.bed"



######## After lifting over, filtering and stitching fragments

script_path=/net/gs/vol1/home/cxqiu/bin/python_script
work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/human_mouse_enhancer

mamm=Mus_musculus
target_mamm=Homo_sapiens

celltype_list=(Adipocyte_cells Adipocyte_cells_Cyp2e1 B_cells Brain_capillary_endothelial_cells CNS_neurons Cardiomyocytes Corticofugal_neurons Endocardial_cells Endothelium Epithelial_cells Erythroid_cells Eye Glia Glomerular_endothelial_cells Gut_epithelial_cells Hepatocytes Intermediate_neuronal_progenitors Kidney Lateral_plate_and_intermediate_mesoderm Liver_sinusoidal_endothelial_cells Lung_and_airway Lymphatic_vessel_endothelial_cells Melanocyte_cells Mesoderm Neural_crest_PNS_neurons Neuroectoderm_and_glia Olfactory_ensheathing_cells Olfactory_neurons Oligodendrocytes Skeletal_muscle_cells T_cells White_blood_cells)
celltype="${celltype_list[$SGE_TASK_ID - 1]}"
echo $celltype

Rscript "$script_path/stitchHalFrags_v3.R" \
"$work_path/human_mouse_enhancer_liftover/${mamm}_${celltype}.bed" \
"$work_path/human_mouse_enhancer_liftover/liftover_${mamm}_${celltype}_to_${target_mamm}.bed" \
"$work_path/human_mouse_enhancer_liftover/stitchHalFrags_${mamm}_${celltype}_to_${target_mamm}.bed" \
0.1 2.0


############################################################
### Step-3: Overlap liftover regions with the target regions

work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/human_mouse_enhancer

celltype_list=(Adipocyte_cells Adipocyte_cells_Cyp2e1 B_cells Brain_capillary_endothelial_cells CNS_neurons Cardiomyocytes Corticofugal_neurons Endocardial_cells Endothelium Epithelial_cells Erythroid_cells Eye Glia Glomerular_endothelial_cells Gut_epithelial_cells Hepatocytes Intermediate_neuronal_progenitors Kidney Lateral_plate_and_intermediate_mesoderm Liver_sinusoidal_endothelial_cells Lung_and_airway Lymphatic_vessel_endothelial_cells Melanocyte_cells Mesoderm Neural_crest_PNS_neurons Neuroectoderm_and_glia Olfactory_ensheathing_cells Olfactory_neurons Oligodendrocytes Skeletal_muscle_cells T_cells White_blood_cells)

mamm=Mus_musculus
target_mamm=Homo_sapiens

for celltype in "${celltype_list[@]}"; do
    bedtools intersect \
        -a "$work_path/human_mouse_enhancer_liftover/stitchHalFrags_${mamm}_${celltype}_to_${target_mamm}.bed" \
        -b "$work_path/human_mouse_enhancer_liftover/${target_mamm}_${celltype}.bed" \
        -wa -wb \
        > "$work_path/human_mouse_enhancer_liftover/overlap_${mamm}_${target_mamm}_${celltype}.bed"
done

mamm=Homo_sapiens
target_mamm=Mus_musculus

for celltype in "${celltype_list[@]}"; do
    bedtools intersect \
        -a "$work_path/human_mouse_enhancer_liftover/stitchHalFrags_${mamm}_${celltype}_to_${target_mamm}.bed" \
        -b "$work_path/human_mouse_enhancer_liftover/${target_mamm}_${celltype}.bed" \
        -wa -wb \
        > "$work_path/human_mouse_enhancer_liftover/overlap_${mamm}_${target_mamm}_${celltype}.bed"
done



############################################################################################
### Step-4: Asking how many enhancers are species-specific vs. share between human and mouse


model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

mamm = "Mus_musculus"
target_mamm = "Homo_sapiens"

df = NULL

for(i in 1:length(celltype_list)){
    celltype = celltype_list[i]
    x = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_liftover/", mamm, "_", celltype, ".bed"))
    y = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_liftover/stitchHalFrags_", mamm, "_", celltype, "_to_", target_mamm, ".bed"))
    z = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_liftover/overlap_", mamm, "_", target_mamm, "_", celltype, ".bed"))
    x_num = length(unique(x[,4]))
    y_num = length(unique(y[,4]))
    z_num = length(unique(z[,4]))
    df = rbind(df, data.frame(celltype = celltype, num_enhancer = x_num, num_liftover = y_num, num_overlap = z_num))
}

x = round(100*sum(df$num_overlap)/sum(df$num_enhancer))
print(paste0(x, "% of the mouse enhancers had a human ortholog called as an enhancer for the same cell type"))
x = round(100*(sum(df$num_liftover) - sum(df$num_overlap))/sum(df$num_enhancer))
print(paste0(x, "% had a human ortholog that did not meet this criterion"))
x = round(100*(sum(df$num_enhancer) - sum(df$num_liftover))/sum(df$num_enhancer))
print(paste0(x, "% could not be lifted over"))


### Following bidirectional liftover, 23% of mouse enhancers had a human ortholog called as an enhancer for the same cell type, 49% had a human ortholog not matched to the same cell type, and 28% could not be lifted over

mamm = "Homo_sapiens"
target_mamm = "Mus_musculus"

df = NULL

for(i in 1:length(celltype_list)){
    celltype = celltype_list[i]
    x = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_liftover/", mamm, "_", celltype, ".bed"))
    y = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_liftover/stitchHalFrags_", mamm, "_", celltype, "_to_", target_mamm, ".bed"))
    z = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/human_mouse_enhancer/human_mouse_enhancer_liftover/overlap_", mamm, "_", target_mamm, "_", celltype, ".bed"))
    x_num = length(unique(x[,4]))
    y_num = length(unique(y[,4]))
    z_num = length(unique(z[,4]))
    df = rbind(df, data.frame(celltype = celltype, num_enhancer = x_num, num_liftover = y_num, num_overlap = z_num))
}

x = round(100*sum(df$num_overlap)/sum(df$num_enhancer))
print(paste0(x, "% of the human enhancers had a mouse ortholog called as an enhancer for the same cell type"))
x = round(100*(sum(df$num_liftover) - sum(df$num_overlap))/sum(df$num_enhancer))
print(paste0(x, "% had a human ortholog that did not meet this criterion"))
x = round(100*(sum(df$num_enhancer) - sum(df$num_liftover))/sum(df$num_enhancer))
print(paste0(x, "% could not be lifted over"))

### a distribution nearly identical in the human-to-mouse direction (24%, 49%, 27%). 

