

##########################################################
### Focusing on Afp TSS region predicted by STEAM-v1 model


### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu




##############################
### Step-1: perform prediction

import sys, os
import anndata as ad
import crested
import numpy as np
import pandas as pd
import gzip
import keras
import pysam
from collections import defaultdict

work_path = ""
model_id = ""

model_path = f"STEAM-v1-model.keras"
model = keras.models.load_model(model_path, compile=False)

#celltype_list = ["Hepatocytes", "Erythroid_cells"]
celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]

adata = ad.read_h5ad(f"{work_path}/data_window_cluster_mouse.top3K.h5ad")

with open(f"{work_path}/{model_id}/Afp_TSS/candidate_liftover_regions.txt") as f:
    mamm_list = [line.rstrip().split("\t")[3] for line in f if line.rstrip().split("\t")[3] != "species"]

for mamm in mamm_list:
    print(mamm)
    regions = []
    candidate_regions = []
    with open(f"{work_path}/{model_id}/Afp_TSS/prediction_multispecies/{mamm}.txt") as f:
        for line in f:
            chr, start, end = line.rstrip().split('\t')
            regions.append((chr, int(start) + 557, int(end) - 557))
            candidate_regions.append((chr, int(start), int(end)))
    predictions = []
    for iter in range(1,11):
        fasta = pysam.FastaFile(f"{work_path}/genome/fasta_trf_Afp/{mamm}_trf_{iter}.fa")
        candidate_sequences = []
        for region in candidate_regions:
            candidate_sequences.append(fasta.fetch(region[0], region[1], region[2]))
        predictions.append(crested.tl.predict(input=candidate_sequences, model=model))
    predictions_arr = np.stack(predictions, axis = 0)
    predictions_mean = np.mean(predictions_arr, axis = 0)
    for celltype in celltype_list:
        celltype_id = np.where(adata.obs.index == celltype)[0][0]
        scores = defaultdict(list)
        for i, (chrom, start, end) in enumerate(regions):
            for pos in range(start, end, 100):
                scores[(chrom, pos)].append(predictions_mean[i, celltype_id])
        rows = []
        for (chrom, pos), vals in sorted(scores.items()):
            rows.append((chrom, pos, pos + 100, np.mean(vals)))
        df = pd.DataFrame(rows, columns=["chr", "start", "end", "score"])
        df.to_csv(f"{work_path}/{model_id}/Afp_TSS/prediction_multispecies/{mamm}_100bp_{celltype}.tsv", sep="\t", index=False, header=False)



######################################################
### Step-2: counting how many elements in each species


library(tidyr)
library(GenomicRanges)
library(dplyr)

dat_phred = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/dat.txt.rds"))

mamm_list = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/species_order.txt"))
mamm_list = as.vector(mamm_list$V1)

region = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/candidate_liftover_regions.txt"), header=T)

celltype_list=c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

dat = NULL
for(celltype in celltype_list){
    print(celltype)
    values = dat_phred[, colnames(dat_phred) == celltype]
    sorted = sort(values)
    N = length(sorted)

    Q_from_x <- function(x) {
        r <- findInterval(x, sorted)
        q <- (r - 0.5) / N
        q <- pmin(q, 1 - 1e-12)
        -10 * log10(1 - q)
    }

    for(mamm in mamm_list){
        dat_i = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_multispecies/", mamm, "_100bp_", celltype, ".tsv"))
        colnames(dat_i) = c("chr", "start", "end", "score")
        dat_i$phred_score = round(Q_from_x(dat_i$score), 3)
        dat_i$species = mamm
        dat_i$celltype = celltype
        dat = rbind(dat, dat_i[,c("chr", "start", "end", "phred_score", "species", "celltype")])
    }
}

saveRDS(dat, paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_multispecies_res.rds"))





###################################################
### Step-3: Identify core region for Afp TSS region 


library(tidyr)
library(GenomicRanges)
library(dplyr)
options(scipen = 999)

dat = readRDS(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_multispecies_res.rds"))

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

for(celltype in celltype_list){
    print(celltype)

    dat_sub = dat[dat$phred_score > 24.5 & dat$celltype == celltype,]

    gr <- makeGRangesFromDataFrame(dat_sub[,c(1,2,3,5)],
                                    keep.extra.columns = TRUE,
                                    seqnames.field = "chr",
                                    start.field = "start",
                                    end.field = "end")

    gr_split <- split(gr, gr$species)
    merged   <- reduce(gr_split)

    result <- as.data.frame(merged) %>%
      rename(chr = seqnames, species = group_name) %>%
      select(-width, -strand, -group)

    mamm_list = unique(result$species)
    for(mamm in mamm_list){
        print(mamm)
        res = result[result$species == mamm,c(2:4)]
        write.table(res, paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/core_region/", mamm, ".", celltype, ".bed"), row.names=F, col.names=F, sep="\t", quote=F)
    }
}


##############################################################################
### Identifying the core regions

import random
import sys, os
import numpy as np
import gzip
import pysam
from collections import Counter

import anndata as ad
import crested
import keras

work_path = ""

########################################
# Constants
########################################

BASES = np.array(['A','C','G','T'])
BASE_TO_INT = {b:i for i,b in enumerate(BASES)}

model_id = "mouse_fake_track_15"
step_size = 20
seq_length = 2114
n_reps = 10

########################################
# FASTA + model
########################################

model_path = f"STEAM-v1-model.keras"
model = keras.models.load_model(model_path, compile=False)

adata = ad.read_h5ad(f"{work_path}/{model_id}/data_window_cluster_mouse.top3K.h5ad")

celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]
celltype = celltype_list[int(sys.argv[1]) - 1]
celltype_index = np.where(np.isin(adata.obs.index, celltype))[0]

########################################
# Dinucleotide model
########################################

def prepare_dinuc_probabilities(seq):
    seq = seq.upper()
    counts = Counter()
    total = 0
    for a, b in zip(seq[:-1], seq[1:]):
        if a in BASE_TO_INT and b in BASE_TO_INT:
            counts[a+b] += 1
            total += 1
    dinuc_freq = {
        a+b: counts.get(a+b, 0) / total
        for a in BASES for b in BASES
    }
    mono = np.zeros(4)
    for i, a in enumerate(BASES):
        mono[i] = sum(dinuc_freq[a+b] for b in BASES)
    mono /= mono.sum()
    cond = np.zeros((4,4))
    for i,a in enumerate(BASES):
        row = np.array([dinuc_freq[a+b] for b in BASES])
        cond[i] = row / row.sum() if row.sum() > 0 else 0.25
    return mono, cond

########################################
# FAST dinucleotide sequence generator
########################################

def generate_dinuc_sequence(length, first_probs, cond_mat):
    if length == 0:
        return ""
    seq = np.empty(length, dtype=np.int8)
    seq[0] = np.random.choice(4, p=first_probs)
    r = np.random.rand(length - 1)
    cdf = np.cumsum(cond_mat, axis=1)
    for i in range(1, length):
        seq[i] = np.searchsorted(cdf[seq[i-1]], r[i-1])
    return ''.join(BASES[seq])

########################################
# Precompute replacement schedule
########################################

replacement_schedule = []
left = right = 0
while left + right < seq_length:
    if left == right:
        left += step_size
    else:
        right += step_size
    if left + right > seq_length:
        break
    replacement_schedule.append((left, right))

########################################
# Firt step of trimming
########################################

mamm_list = ["Otolemur_garnettii","Daubentonia_madagascariensis","Propithecus_coquereli","Cheirogaleus_medius","Microcebus_murinus","Eulemur_flavifrons","Pithecia_pithecia","Saguinus_imperator","Callithrix_jacchus","Aotus_nancymaae","Saimiri_boliviensis","Cebus_albifrons","Cebus_capucinus","Mandrillus_leucophaeus","Cercocebus_atys","Papio_anubis","Macaca_mulatta","Macaca_fascicularis","Macaca_nemestrina","Chlorocebus_sabaeus","Rhinopithecus_bieti","Piliocolobus_tephrosceles","Colobus_angolensis","Nomascus_leucogenys","Gorilla_gorilla","Pan_troglodytes","Pan_paniscus","Homo_sapiens","Pongo_abelii","Tupaia_chinensis","Oryctolagus_cuniculus","Ochotona_princeps","Ctenodactylus_gundi","Heterocephalus_glaber","Fukomys_damarensis","Cavia_porcellus","Cavia_aperea","Octodon_degus","Myocastor_coypus","Rattus_norvegicus","Mus_musculus","Mus_spretus","Mus_caroli","Mus_pahari","Meriones_unguiculatus","Mesocricetus_auratus","Cricetulus_griseus","Microtus_ochrogaster","Peromyscus_maniculatus","Cricetomys_gambianus","Nannospalax_galili","Jaculus_jaculus","Dipodomys_ordii","Spermophilus_dauricus","Ictidomys_tridecemlineatus","Marmota_marmota","Erinaceus_europaeus","Sorex_araneus","Condylura_cristata","Megaderma_lyra","Hipposideros_armiger","Rhinolophus_sinicus","Macroglossus_sobrinus","Pteropus_vampyrus","Pteropus_alecto","Anoura_caudifer","Micronycteris_hirsuta","Desmodus_rotundus","Eptesicus_fuscus","Myotis_brandtii","Myotis_lucifugus","Myotis_davidii","Miniopterus_natalensis","Vicugna_pacos","Camelus_dromedarius","Camelus_ferus","Camelus_bactrianus","Eubalaena_japonica","Balaenoptera_acutorostrata","Eschrichtius_robustus","Lipotes_vexillifer","Neophocaena_asiaeorientalis","Delphinapterus_leucas","Monodon_monoceros","Tursiops_truncatus","Orcinus_orca","Tragulus_javanicus","Bos_taurus","Bos_indicus","Bos_mutus","Bison_bison","Hemitragus_hylocrius","Capra_hircus","Capra_aegagrus","Ovis_aries","Ovis_canadensis","Antilocapra_americana","Odocoileus_virginianus","Rangifer_tarandus","Elaphurus_davidianus","Sus_scrofa","Suricata_suricatta","Paradoxurus_hermaphroditus","Panthera_tigris","Panthera_pardus","Panthera_onca","Felis_catus","Puma_concolor","Acinonyx_jubatus","Lycaon_pictus","Canis_lupus","Canis_lupus_familiaris","Vulpes_lagopus","Enhydra_lutris","Mustela_putorius","Ailurus_fulgens","Odobenus_rosmarus","Leptonychotes_weddellii","Neomonachus_schauinslandi","Ailuropoda_melanoleuca","Ursus_maritimus","Manis_javanica","Equus_caballus","Equus_przewalskii","Equus_asinus","Tapirus_indicus","Tapirus_terrestris","Dicerorhinus_sumatrensis","Diceros_bicornis","Ceratotherium_simum","Dasypus_novemcinctus","Choloepus_hoffmanni","Trichechus_manatus","Loxodonta_africana","Echinops_telfairi","Elephantulus_edwardii"]

for mamm in mamm_list:

    if os.path.exists(f"{work_path}/{model_id}/Afp_TSS/core_region/{mamm}.{celltype}.bed"):

        fasta = pysam.FastaFile(f"{work_path}/genome/{mamm}/{mamm}.fa")

        chr_size = {}
        with open(f"{work_path}/genome/{mamm}/{mamm}.chrom.sizes.update") as f:
            for line in f:
                c, s = line.rstrip().split('\t')
                chr_size[c] = int(s)

        all_regions = []
        with open(f"{work_path}/{model_id}/Afp_TSS/core_region/{mamm}.{celltype}.bed") as f:
            for line in f:
                l = line.rstrip().split('\t')
                loc = (int(l[1]) + int(l[2])) // 2
                if loc-1057 > 0 and loc+1057 < chr_size[l[0]]:
                    all_regions.append((l[0], int(l[1]), int(l[2])))

        result = []
        for chrom, s, e in all_regions:
            loc = (s + e) // 2
            local_seq = fasta.fetch(chrom, max(1, loc - 2500), min(chr_size[chrom], loc + 2500))
            first_probs, cond_mat = prepare_dinuc_probabilities(local_seq)
            orig_seq = fasta.fetch(chrom, loc-1057, loc+1057)
            orig_score = float(crested.tl.predict(input=orig_seq, model=model)[:,celltype_index].item())
            result.append((chrom, s, e, 0, 0, orig_score))
            for left, right in replacement_schedule:
                middle = orig_seq[left:seq_length-right]
                candidate_sequences = []
                for _ in range(n_reps):
                    lp = generate_dinuc_sequence(left, first_probs, cond_mat)
                    rp = generate_dinuc_sequence(right, first_probs, cond_mat)
                    candidate_sequences.append(lp + middle + rp)
                current_score = float(np.median(crested.tl.predict(input=candidate_sequences, model=model)[:,celltype_index]))
                result.append((chrom, s, e, left, right, current_score))
                if current_score < orig_score * 0.9:
                    if left == right:
                        right = right - 20
                        rp = generate_dinuc_sequence(right, first_probs, cond_mat)
                        for x in range(left + step_size, seq_length - right, step_size):
                            middle = orig_seq[x:seq_length-right]
                            candidate_sequences = []
                            for _ in range(n_reps):
                                lp = generate_dinuc_sequence(x, first_probs, cond_mat)
                                candidate_sequences.append(lp + middle + rp)
                            current_score = float(np.median(crested.tl.predict(input=candidate_sequences, model=model)[:,celltype_index]))
                            result.append((chrom, s, e, x, right, current_score))
                            if current_score < orig_score * 0.9:
                                break
                    else:
                        left = left - 20
                        lp = generate_dinuc_sequence(left, first_probs, cond_mat)
                        for x in range(right + step_size, seq_length - left, step_size):
                            middle = orig_seq[left:seq_length-x]
                            candidate_sequences = []
                            for _ in range(n_reps):
                                rp = generate_dinuc_sequence(x, first_probs, cond_mat)
                                candidate_sequences.append(lp + middle + rp)
                            current_score = float(np.median(crested.tl.predict(input=candidate_sequences, model=model)[:,celltype_index]))
                            result.append((chrom, s, e, left, x, current_score))
                            if current_score < orig_score * 0.9:
                                break
                    break

        with open(f"{work_path}/{model_id}/Afp_TSS/core_region/{mamm}.{celltype}.core_region.txt", "w") as out:
            for i in result:
                out.write('\t'.join(map(str, i)) + '\n')


##################################
#### Identify core regions



celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

mamm_list = c("Otolemur_garnettii","Daubentonia_madagascariensis","Propithecus_coquereli","Cheirogaleus_medius","Microcebus_murinus","Eulemur_flavifrons","Pithecia_pithecia","Saguinus_imperator","Callithrix_jacchus","Aotus_nancymaae","Saimiri_boliviensis","Cebus_albifrons","Cebus_capucinus","Mandrillus_leucophaeus","Cercocebus_atys","Papio_anubis","Macaca_mulatta","Macaca_fascicularis","Macaca_nemestrina","Chlorocebus_sabaeus","Rhinopithecus_bieti","Piliocolobus_tephrosceles","Colobus_angolensis","Nomascus_leucogenys","Gorilla_gorilla","Pan_troglodytes","Pan_paniscus","Homo_sapiens","Pongo_abelii","Tupaia_chinensis","Oryctolagus_cuniculus","Ochotona_princeps","Ctenodactylus_gundi","Heterocephalus_glaber","Fukomys_damarensis","Cavia_porcellus","Cavia_aperea","Octodon_degus","Myocastor_coypus","Rattus_norvegicus","Mus_musculus","Mus_spretus","Mus_caroli","Mus_pahari","Meriones_unguiculatus","Mesocricetus_auratus","Cricetulus_griseus","Microtus_ochrogaster","Peromyscus_maniculatus","Cricetomys_gambianus","Nannospalax_galili","Jaculus_jaculus","Dipodomys_ordii","Spermophilus_dauricus","Ictidomys_tridecemlineatus","Marmota_marmota","Erinaceus_europaeus","Sorex_araneus","Condylura_cristata","Megaderma_lyra","Hipposideros_armiger","Rhinolophus_sinicus","Macroglossus_sobrinus","Pteropus_vampyrus","Pteropus_alecto","Anoura_caudifer","Micronycteris_hirsuta","Desmodus_rotundus","Eptesicus_fuscus","Myotis_brandtii","Myotis_lucifugus","Myotis_davidii","Miniopterus_natalensis","Vicugna_pacos","Camelus_dromedarius","Camelus_ferus","Camelus_bactrianus","Eubalaena_japonica","Balaenoptera_acutorostrata","Eschrichtius_robustus","Lipotes_vexillifer","Neophocaena_asiaeorientalis","Delphinapterus_leucas","Monodon_monoceros","Tursiops_truncatus","Orcinus_orca","Tragulus_javanicus","Bos_taurus","Bos_indicus","Bos_mutus","Bison_bison","Hemitragus_hylocrius","Capra_hircus","Capra_aegagrus","Ovis_aries","Ovis_canadensis","Antilocapra_americana","Odocoileus_virginianus","Rangifer_tarandus","Elaphurus_davidianus","Sus_scrofa","Suricata_suricatta","Paradoxurus_hermaphroditus","Panthera_tigris","Panthera_pardus","Panthera_onca","Felis_catus","Puma_concolor","Acinonyx_jubatus","Lycaon_pictus","Canis_lupus","Canis_lupus_familiaris","Vulpes_lagopus","Enhydra_lutris","Mustela_putorius","Ailurus_fulgens","Odobenus_rosmarus","Leptonychotes_weddellii","Neomonachus_schauinslandi","Ailuropoda_melanoleuca","Ursus_maritimus","Manis_javanica","Equus_caballus","Equus_przewalskii","Equus_asinus","Tapirus_indicus","Tapirus_terrestris","Dicerorhinus_sumatrensis","Diceros_bicornis","Ceratotherium_simum","Dasypus_novemcinctus","Choloepus_hoffmanni","Trichechus_manatus","Loxodonta_africana","Echinops_telfairi","Elephantulus_edwardii")

dat_all = NULL
for(celltype in celltype_list){
    print(celltype)
    for(mamm in mamm_list){
        file_name = paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/core_region/", mamm, ".", celltype, ".core_region.txt")
        if(file.exists(file_name) && file.info(file_name)$size > 0){
            dat = read.table(file_name)
            colnames(dat) = c("chr", "start", "end", "left", "right", "score")
            dat_orig = dat %>% filter(left == 0, right == 0) %>% rename(orig_score = score)
            dat = dat %>% left_join(dat_orig[,c("chr", "start", "end", "orig_score")], by = c("chr", "start", "end"))
            dat = dat %>% filter(score > orig_score * 0.9)
            dat$min_len = dat$left + dat$right
            dat = dat %>% group_by(chr, start, end) %>% slice_max(order_by = min_len, n = 1)
            dat$mid = (dat$start + dat$end)/2
            dat$core_start = dat$mid - 1057 + dat$left
            dat$core_end = dat$mid + 1057 - dat$right
            dat = dat[,c("chr", "start", "end", "mid","core_start", "core_end", "score")]
            dat$species = mamm
            dat$celltype = celltype
            dat_all = rbind(dat_all, dat)
            if (celltype == "Hepatocytes"){
                dat = dat[,c("chr", "core_start", "core_end")]
                dat$species = paste0(mamm, "_", 1:nrow(dat))
                write.table(dat, paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/core_region_2/", mamm, ".", celltype, ".bed"), row.names=F, col.names=F, quote=F, sep="\t")
            }
        }
    }
}

saveRDS(dat_all, paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_multispecies.rds"))










