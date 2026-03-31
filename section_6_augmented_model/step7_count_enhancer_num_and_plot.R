

#################################################################
### count the number of enhancers in each cell type and each mamm

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
library(tidyr)
options(scipen = 999)

celltype_list = c("Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells")

mamm_list = c("Otolemur_garnettii","Daubentonia_madagascariensis","Propithecus_coquereli","Cheirogaleus_medius","Microcebus_murinus","Eulemur_flavifrons","Pithecia_pithecia","Saguinus_imperator","Callithrix_jacchus","Aotus_nancymaae","Saimiri_boliviensis","Cebus_albifrons","Cebus_capucinus","Mandrillus_leucophaeus","Cercocebus_atys","Papio_anubis","Macaca_mulatta","Macaca_fascicularis","Macaca_nemestrina","Chlorocebus_sabaeus","Rhinopithecus_bieti","Piliocolobus_tephrosceles","Colobus_angolensis","Nomascus_leucogenys","Gorilla_gorilla","Pan_troglodytes","Pan_paniscus","Homo_sapiens","Pongo_abelii","Tupaia_chinensis","Oryctolagus_cuniculus","Ochotona_princeps","Ctenodactylus_gundi","Heterocephalus_glaber","Fukomys_damarensis","Cavia_porcellus","Cavia_aperea","Octodon_degus","Myocastor_coypus","Rattus_norvegicus","Mus_musculus","Mus_spretus","Mus_caroli","Mus_pahari","Meriones_unguiculatus","Mesocricetus_auratus","Cricetulus_griseus","Microtus_ochrogaster","Peromyscus_maniculatus","Cricetomys_gambianus","Nannospalax_galili","Jaculus_jaculus","Dipodomys_ordii","Spermophilus_dauricus","Ictidomys_tridecemlineatus","Marmota_marmota","Erinaceus_europaeus","Sorex_araneus","Condylura_cristata","Megaderma_lyra","Hipposideros_armiger","Rhinolophus_sinicus","Macroglossus_sobrinus","Pteropus_vampyrus","Pteropus_alecto","Anoura_caudifer","Micronycteris_hirsuta","Desmodus_rotundus","Eptesicus_fuscus","Myotis_brandtii","Myotis_lucifugus","Myotis_davidii","Miniopterus_natalensis","Vicugna_pacos","Camelus_dromedarius","Camelus_ferus","Camelus_bactrianus","Eubalaena_japonica","Balaenoptera_acutorostrata","Eschrichtius_robustus","Lipotes_vexillifer","Neophocaena_asiaeorientalis","Delphinapterus_leucas","Monodon_monoceros","Tursiops_truncatus","Orcinus_orca","Tragulus_javanicus","Bos_taurus","Bos_indicus","Bos_mutus","Bison_bison","Hemitragus_hylocrius","Capra_hircus","Capra_aegagrus","Ovis_aries","Ovis_canadensis","Antilocapra_americana","Odocoileus_virginianus","Rangifer_tarandus","Elaphurus_davidianus","Sus_scrofa","Suricata_suricatta","Paradoxurus_hermaphroditus","Panthera_tigris","Panthera_pardus","Panthera_onca","Felis_catus","Puma_concolor","Acinonyx_jubatus","Lycaon_pictus","Canis_lupus","Canis_lupus_familiaris","Vulpes_lagopus","Enhydra_lutris","Mustela_putorius","Ailurus_fulgens","Odobenus_rosmarus","Leptonychotes_weddellii","Neomonachus_schauinslandi","Ailuropoda_melanoleuca","Ursus_maritimus","Manis_javanica","Equus_caballus","Equus_przewalskii","Equus_asinus","Tapirus_indicus","Tapirus_terrestris","Dicerorhinus_sumatrensis","Diceros_bicornis","Ceratotherium_simum","Dasypus_novemcinctus","Choloepus_hoffmanni","Trichechus_manatus","Loxodonta_africana","Echinops_telfairi","Elephantulus_edwardii")

### aware model
dat = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_aware.rds"))

df = dat %>% group_by(species, celltype) %>% tally() %>% ungroup() %>%
  complete(species = mamm_list, celltype = celltype_list, fill = list(n = 0))

df_x = df[df$celltype == "Hepatocytes",]
print(paste0("Hepatocytes: mean = ", round(mean(df_x$n),1), ", sd = ", round(sd(df_x$n),1)))

df_x = df[df$celltype != "Hepatocytes",]
print(paste0("Others: mean = ", round(mean(df_x$n),1), ", sd = ", round(sd(df_x$n),1)))


### multispecies model
dat = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_multispecies.rds"))

df = dat %>% group_by(species, celltype) %>% tally() %>% ungroup() %>%
  complete(species = mamm_list, celltype = celltype_list, fill = list(n = 0))

df_x = df[df$celltype == "Hepatocytes",]
print(paste0("Hepatocytes: mean = ", round(mean(df_x$n),1), ", sd = ", round(sd(df_x$n),1)))

df_x = df[df$celltype != "Hepatocytes",]
print(paste0("Others: mean = ", round(mean(df_x$n),1), ", sd = ", round(sd(df_x$n),1)))


############## focusing on "Rodentia" only

node_category = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/species_category.txt"))
colnames(node_category) = c("species", "species_category")

dat = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_aware.rds"))

df = dat %>% group_by(species, celltype) %>% tally() %>% ungroup() %>%
  complete(species = mamm_list, celltype = celltype_list, fill = list(n = 0)) %>% 
  left_join(node_category, by = "species")

df_x = df[df$celltype == "Hepatocytes" & df$species_category == "Rodentia",]
print(paste0("Hepatocytes: mean = ", round(mean(df_x$n),1), ", sd = ", round(sd(df_x$n),1)))


dat = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_multispecies.rds"))

df = dat %>% group_by(species, celltype) %>% tally() %>% ungroup() %>%
  complete(species = mamm_list, celltype = celltype_list, fill = list(n = 0)) %>% 
  left_join(node_category, by = "species")

df_x = df[df$celltype == "Hepatocytes" & df$species_category == "Rodentia",]
print(paste0("Hepatocytes: mean = ", round(mean(df_x$n),1), ", sd = ", round(sd(df_x$n),1)))



###########################################################
### subset the phylogenetic tree to 136 species and plot it

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
import matplotlib.pyplot as plt
import io

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
model_id = "mouse_fake_track_15"

mamm_list = []
with open(f"{work_path}/{model_id}/Afp_TSS/candidate_liftover_regions.txt") as f:
    line = f.readline()
    line = f.readline()
    while line:
        l = line.rstrip().split('\t')
#        if (int(l[4]) - int(l[1])) >= 80000 and (int(l[2]) - int(l[4])) >= 20000:
#            mamm_list.append(l[3])
        mamm_list.append(l[3])
        line = f.readline()

# Load the full tree
tree = Phylo.read(f"{work_path}/{model_id}/Afp_TSS/241-mammalian-2020v2.phast-242.nh", "newick")
species = [c.name for c in tree.get_terminals()]

to_remove = [c.name for c in tree.get_terminals() if c.name not in mamm_list]
for species in to_remove:
    tree.prune(species)

for node in tree.get_nonterminals():
    node.name = None

# Plot vertically (top to bottom)
fig, ax = plt.subplots(figsize=(5, len(mamm_list) * 0.15))
Phylo.draw(tree, axes=ax, do_show=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xticks([])
ax.set_xlabel("")
ax.set_yticks([])
ax.set_ylabel("")
plt.tight_layout()
plt.savefig("/net/gs/vol1/home/cxqiu/share/phylo_tree.pdf")

# Save top-to-bottom order of species
order = [c.name for c in tree.get_terminals()]
with open(f"{work_path}/{model_id}/Afp_TSS/species_order.txt", "w") as f:
    for s in order:
        f.write(s + "\n")

len(order)


for node in tree.get_terminals():
    node.name = None

# Plot vertically (top to bottom)
fig, ax = plt.subplots(figsize=(5, len(mamm_list) * 0.15))
Phylo.draw(tree, axes=ax, do_show=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xticks([])
ax.set_xlabel("")
ax.set_yticks([])
ax.set_ylabel("")
plt.tight_layout()
plt.savefig("/net/gs/vol1/home/cxqiu/share/phylo_tree.no_label.pdf")

def equalize_depths(tree):
    max_depth = max(tree.distance(tip) for tip in tree.get_terminals())
    for tip in tree.get_terminals():
        depth = tree.distance(tip)
        tip.branch_length += (max_depth - depth)

equalize_depths(tree)

# Plot vertically (top to bottom)
fig, ax = plt.subplots(figsize=(5, len(mamm_list) * 0.15))
Phylo.draw(tree, axes=ax, do_show=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xticks([])
ax.set_xlabel("")
ax.set_yticks([])
ax.set_ylabel("")
plt.tight_layout()
plt.savefig("/net/gs/vol1/home/cxqiu/share/phylo_tree.no_label_2.pdf")



##############################################################################
### check if the liftover regions are flipped from the original mouse sequence

from Bio import pairwise2
from Bio.Seq import Seq
import pysam
import os, sys

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested"
model_id = "mouse_fake_track_15"

mamm_list = ["Acinonyx_jubatus","Acomys_cahirinus","Ailuropoda_melanoleuca","Ailurus_fulgens","Allactaga_bullata","Alouatta_palliata","Ammotragus_lervia","Anoura_caudifer","Antilocapra_americana","Aotus_nancymaae","Aplodontia_rufa","Artibeus_jamaicensis","Ateles_geoffroyi","Balaenoptera_acutorostrata","Balaenoptera_bonaerensis","Beatragus_hunteri","Bison_bison","Bos_indicus","Bos_mutus","Bos_taurus","Bubalus_bubalis","Callicebus_donacophilus","Callithrix_jacchus","Camelus_bactrianus","Camelus_dromedarius","Camelus_ferus","Canis_lupus","Canis_lupus_familiaris","Capra_aegagrus","Capra_hircus","Capromys_pilorides","Carollia_perspicillata","Castor_canadensis","Catagonus_wagneri","Cavia_aperea","Cavia_porcellus","Cavia_tschudii","Cebus_albifrons","Cebus_capucinus","Ceratotherium_simum","Ceratotherium_simum_cottoni","Cercocebus_atys","Cercopithecus_neglectus","Chaetophractus_vellerosus","Cheirogaleus_medius","Chinchilla_lanigera","Chlorocebus_sabaeus","Choloepus_didactylus","Choloepus_hoffmanni","Chrysochloris_asiatica","Colobus_angolensis","Condylura_cristata","Craseonycteris_thonglongyai","Cricetomys_gambianus","Cricetulus_griseus","Crocidura_indochinensis","Cryptoprocta_ferox","Ctenodactylus_gundi","Ctenomys_sociabilis","Cuniculus_paca","Dasyprocta_punctata","Dasypus_novemcinctus","Daubentonia_madagascariensis","Delphinapterus_leucas","Desmodus_rotundus","Dicerorhinus_sumatrensis","Diceros_bicornis","Dinomys_branickii","Dipodomys_ordii","Dipodomys_stephensi","Dolichotis_patagonum","Echinops_telfairi","Eidolon_helvum","Elaphurus_davidianus","Elephantulus_edwardii","Ellobius_lutescens","Ellobius_talpinus","Enhydra_lutris","Eptesicus_fuscus","Equus_asinus","Equus_caballus","Equus_przewalskii","Erinaceus_europaeus","Erythrocebus_patas","Eschrichtius_robustus","Eubalaena_japonica","Eulemur_flavifrons","Eulemur_fulvus","Felis_catus","Felis_nigripes","Fukomys_damarensis","Galeopterus_variegatus","Giraffa_tippelskirchi","Glis_glis","Gorilla_gorilla","Graphiurus_murinus","Helogale_parvula","Hemitragus_hylocrius","Heterocephalus_glaber","Heterohyrax_brucei","Hippopotamus_amphibius","Hipposideros_armiger","Hipposideros_galeritus","Homo_sapiens","Hyaena_hyaena","Hydrochoerus_hydrochaeris","Hystrix_cristata","Ictidomys_tridecemlineatus","Indri_indri","Inia_geoffrensis","Jaculus_jaculus","Kogia_breviceps","Lasiurus_borealis","Lemur_catta","Leptonychotes_weddellii","Lepus_americanus","Lipotes_vexillifer","Loxodonta_africana","Lycaon_pictus","Macaca_fascicularis","Macaca_mulatta","Macaca_nemestrina","Macroglossus_sobrinus","Mandrillus_leucophaeus","Manis_javanica","Manis_pentadactyla","Marmota_marmota","Megaderma_lyra","Mellivora_capensis","Meriones_unguiculatus","Mesocricetus_auratus","Mesoplodon_bidens","Microcebus_murinus","Microgale_talazaci","Micronycteris_hirsuta","Microtus_ochrogaster","Miniopterus_natalensis","Miniopterus_schreibersii","Mirounga_angustirostris","Mirza_coquereli","Monodon_monoceros","Mormoops_blainvillei","Moschus_moschiferus","Mungos_mungo","Murina_feae","Mus_caroli","Mus_pahari","Mus_spretus","Muscardinus_avellanarius","Mustela_putorius","Myocastor_coypus","Myotis_brandtii","Myotis_davidii","Myotis_lucifugus","Myotis_myotis","Myrmecophaga_tridactyla","Nannospalax_galili","Nasalis_larvatus","Neomonachus_schauinslandi","Neophocaena_asiaeorientalis","Noctilio_leporinus","Nomascus_leucogenys","Nycticebus_coucang","Ochotona_princeps","Octodon_degus","Odobenus_rosmarus","Odocoileus_virginianus","Okapia_johnstoni","Ondatra_zibethicus","Onychomys_torridus","Orcinus_orca","Orycteropus_afer","Oryctolagus_cuniculus","Otolemur_garnettii","Ovis_aries","Ovis_canadensis","Pan_paniscus","Pan_troglodytes","Panthera_onca","Panthera_pardus","Panthera_tigris","Pantholops_hodgsonii","Papio_anubis","Paradoxurus_hermaphroditus","Perognathus_longimembris","Peromyscus_maniculatus","Petromus_typicus","Phocoena_phocoena","Piliocolobus_tephrosceles","Pipistrellus_pipistrellus","Pithecia_pithecia","Platanista_gangetica","Pongo_abelii","Procavia_capensis","Propithecus_coquereli","Psammomys_obesus","Pteronotus_parnellii","Pteronura_brasiliensis","Pteropus_alecto","Pteropus_vampyrus","Puma_concolor","Pygathrix_nemaeus","Rangifer_tarandus","Rattus_norvegicus","Rhinolophus_sinicus","Rhinopithecus_bieti","Rhinopithecus_roxellana","Rousettus_aegyptiacus","Saguinus_imperator","Saiga_tatarica","Saimiri_boliviensis","Scalopus_aquaticus","Semnopithecus_entellus","Sigmodon_hispidus","Solenodon_paradoxus","Sorex_araneus","Spermophilus_dauricus","Spilogale_gracilis","Suricata_suricatta","Sus_scrofa","Tadarida_brasiliensis","Tamandua_tetradactyla","Tapirus_indicus","Tapirus_terrestris","Thryonomys_swinderianus","Tolypeutes_matacus","Tonatia_saurophila","Tragulus_javanicus","Trichechus_manatus","Tupaia_chinensis","Tupaia_tana","Tursiops_truncatus","Uropsilus_gracilis","Ursus_maritimus","Vicugna_pacos","Vulpes_lagopus","Xerus_inauris","Zalophus_californianus","Zapus_hudsonius","Ziphius_cavirostris"]

def get_identity(seq1, seq2):
    aln = pairwise2.align.globalms(seq1, seq2, 2, -1, -2, -0.5)[0]
    matches = sum(a == b for a, b in zip(aln[0], aln[1]))
    return matches / len(aln[0])

mamm = "Mus_musculus"
fa = pysam.FastaFile(f"{work_path}/genome/{mamm}/{mamm}.fa")
with open(f"{work_path}/{model_id}/Afp_TSS/liftover/mouse_Afp_TSS.bed") as file:
    chr, start, end, name = file.readline().rstrip().split('\t')
    mouse_seq = fa.fetch(chr, int(start), int(end))

with open(f"{work_path}/{model_id}/Afp_TSS/liftover_flip_check.txt", "w") as out:
    for mamm in mamm_list:
        print(mamm)
        fa = pysam.FastaFile(f"{work_path}/genome/{mamm}/{mamm}.fa")
        with open(f"{work_path}/{model_id}/Afp_TSS/liftover/stitchHalFrags_{mamm}.bed") as file:
            for line in file:
                chr, start, end, name = line.rstrip().split('\t')
                target_seq = fa.fetch(chr, int(start), int(end))
                fwd_id = get_identity(mouse_seq, target_seq)
                rc_id  = get_identity(mouse_seq, str(Seq(target_seq).reverse_complement()))
                if fwd_id >= rc_id:
                    out.write(f"{chr}\t{start}\t{end}\t{mamm}\tsame\t{fwd_id}\n")
                else:
                    out.write(f"{chr}\t{start}\t{end}\t{mamm}\topposite\t{rc_id}\n")





###########################################
### plot the big scatter plot (aware model)

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
library(tidyr)
library(GenomicRanges)

dat_phred = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_14/prediction_mammals/prediction_Mus_musculus_trf/dat.txt.rds"))

mamm_list = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/species_order.txt"))
mamm_list = as.vector(mamm_list$V1)

region = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/candidate_liftover_regions.txt"), header=T)

check_flip = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/liftover_flip_check.txt"), header=F)
check_flip = check_flip[check_flip$V4 %in% mamm_list,]
check_flip = check_flip[,c(4,5)]
colnames(check_flip) = c("species", "strand")
check_flip = rbind(check_flip, data.frame(species = "Mus_musculus", strand = "same"))

for(celltype in c("Hepatocytes", "Erythroid_cells")){

    values = dat_phred[, colnames(dat_phred) == celltype]
    sorted = sort(values)
    N = length(sorted)

    Q_from_x <- function(x) {
        r <- findInterval(x, sorted)
        q <- (r - 0.5) / N
        q <- pmin(q, 1 - 1e-12)
        -10 * log10(1 - q)
    }

    dat_all = list()
    cnt = 1
    for(mamm in rev(mamm_list)){
        print(mamm)
        dat_i = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_aware/", mamm, "_100bp_", celltype, ".tsv"))
        colnames(dat_i) = c("chr", "start", "end", "score")
        mid = region$mid[region$species == mamm]
        strand = check_flip$strand[check_flip$species == mamm]
        dat_i$tmp = 1:nrow(dat_i)
        dat_j = dat_i[dat_i$start < mid & dat_i$end > mid,]
        if (strand == "same"){
          dat_i$x_axis = dat_i$tmp - dat_j$tmp
        } else {
          dat_i$x_axis = dat_j$tmp - dat_i$tmp
        }
        dat_i$y_axis = cnt
        dat_i$phred_score = Q_from_x(dat_i$score)
        dat_all[[mamm]] = dat_i[,c("x_axis", "y_axis", "phred_score")]
        cnt = cnt + 1
    }

      df = do.call(rbind, dat_all)
      df$phred_score[df$phred_score < 24.5] = 0

      p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_void() +
      theme(legend.position="none")
      ggsave(paste0("~/share/heatmap_", celltype, "_aware.png"), p, width = 5, height = 5, dpi = 300)

      p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE)
      ggsave(paste0("~/share/heatmap_", celltype, "_aware.pdf"), p, width = 5, height = 5)

}






#################################################
### plot the big scatter plot (multispecies model)

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
library(tidyr)
library(GenomicRanges)

dat_phred = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/dat.txt.rds"))

mamm_list = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/species_order.txt"))
mamm_list = as.vector(mamm_list$V1)

region = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/candidate_liftover_regions.txt"), header=T)

check_flip = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/liftover_flip_check.txt"), header=F)
check_flip = check_flip[check_flip$V4 %in% mamm_list,]
check_flip = check_flip[,c(4,5)]
colnames(check_flip) = c("species", "strand")
check_flip = rbind(check_flip, data.frame(species = "Mus_musculus", strand = "same"))

for(celltype in c("Hepatocytes", "Erythroid_cells")){

    values = dat_phred[, colnames(dat_phred) == celltype]
    sorted = sort(values)
    N = length(sorted)

    Q_from_x <- function(x) {
        r <- findInterval(x, sorted)
        q <- (r - 0.5) / N
        q <- pmin(q, 1 - 1e-12)
        -10 * log10(1 - q)
    }

    dat_all = list()
    cnt = 1
    for(mamm in rev(mamm_list)){
        print(mamm)
        dat_i = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_multispecies/", mamm, "_100bp_", celltype, ".tsv"))
        colnames(dat_i) = c("chr", "start", "end", "score")
        mid = region$mid[region$species == mamm]
        strand = check_flip$strand[check_flip$species == mamm]
        dat_i$tmp = 1:nrow(dat_i)
        dat_j = dat_i[dat_i$start < mid & dat_i$end > mid,]
        if (strand == "same"){
          dat_i$x_axis = dat_i$tmp - dat_j$tmp
        } else {
          dat_i$x_axis = dat_j$tmp - dat_i$tmp
        }
        dat_i$y_axis = cnt
        dat_i$phred_score = Q_from_x(dat_i$score)
        dat_all[[mamm]] = dat_i[,c("x_axis", "y_axis", "phred_score")]
        cnt = cnt + 1
    }

      df = do.call(rbind, dat_all)
      df$phred_score[df$phred_score < 24.5] = 0

      p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_void() +
      theme(legend.position="none")
      ggsave(paste0("~/share/heatmap_", celltype, "_multispecies.png"), p, width = 5, height = 5, dpi = 300)

      p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE)
      ggsave(paste0("~/share/heatmap_", celltype, "_multispecies.pdf"), p, width = 5, height = 5)

}





#################################################################################################
### plot the big scatter plot (multispecies model) + original score without Phred like converting

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
library(tidyr)
library(GenomicRanges)

mamm_list = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/species_order.txt"))
mamm_list = as.vector(mamm_list$V1)

region = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/candidate_liftover_regions.txt"), header=T)

check_flip = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/liftover_flip_check.txt"), header=F)
check_flip = check_flip[check_flip$V4 %in% mamm_list,]
check_flip = check_flip[,c(4,5)]
colnames(check_flip) = c("species", "strand")
check_flip = rbind(check_flip, data.frame(species = "Mus_musculus", strand = "same"))

celltype = "Hepatocytes"

dat_all = list()
cnt = 1
for(mamm in rev(mamm_list)){
    print(mamm)
    dat_i = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_multispecies/", mamm, "_100bp_", celltype, ".tsv"))
    colnames(dat_i) = c("chr", "start", "end", "score")
    mid = region$mid[region$species == mamm]
    strand = check_flip$strand[check_flip$species == mamm]
    dat_i$tmp = 1:nrow(dat_i)
    dat_j = dat_i[dat_i$start < mid & dat_i$end > mid,]
    if (strand == "same"){
      dat_i$x_axis = dat_i$tmp - dat_j$tmp
    } else {
      dat_i$x_axis = dat_j$tmp - dat_i$tmp
    }
    dat_i$y_axis = cnt
    dat_all[[mamm]] = dat_i[,c("x_axis", "y_axis", "score")]
    cnt = cnt + 1
}

  df = do.call(rbind, dat_all)

  p = ggplot(df, aes(x_axis, y_axis, fill= score)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position="none",
    panel.background = element_rect(fill = "#440154", color = NA))
  ggsave(paste0("~/share/heatmap_", celltype, "_multispecies_raw.png"), p, width = 5, height = 5, dpi = 300)

  p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE)
  ggsave(paste0("~/share/heatmap_", celltype, "_multispecies.pdf"), p, width = 5, height = 5)

df_sub = df[df$x_axis > (-100) & df$x_axis < 20,]

  p = ggplot(df_sub, aes(x_axis, y_axis, fill= score)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  scale_x_continuous(limits = c(-100, 20), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position="none",
    panel.background = element_rect(fill = "#440154", color = NA))
  ggsave(paste0("~/share/heatmap_", celltype, "_multispecies_raw_zoom_in.png"), p, width = 5, height = 5, dpi = 300)








##################################################################
### BACKUP: trimming to core regions, but its too thin for display, so give up


###########################################
### plot the big scatter plot (aware model)

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
library(tidyr)
library(GenomicRanges)

find_trimmed_overlaps <- function(dat_i, filter_region, return_gr = FALSE) {
  gr_A <- GRanges(seqnames = dat_i$chr, ranges = IRanges(start = dat_i$start, end = dat_i$end), score = dat_i$score, phred_score = dat_i$phred_score)
  gr_B <- GRanges(seqnames = filter_region$chr, ranges = IRanges(start = filter_region$core_start, end = filter_region$core_end))
  hits <- findOverlaps(gr_A, gr_B)
  if (length(hits) == 0) {
    if (return_gr) return(GRanges())
    return(data.frame())
  }
  overlap_gr <- pintersect(
    gr_A[queryHits(hits)],
    gr_B[subjectHits(hits)]
  )
  if (return_gr) return(overlap_gr)
  result <- data.frame(
    chr = as.character(seqnames(overlap_gr)),
    start = start(overlap_gr),
    end = end(overlap_gr),
    score = mcols(overlap_gr)$score,
    phred_score = mcols(overlap_gr)$phred_score
  )
  return(result)
}

dat_phred = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_14/prediction_mammals/prediction_Mus_musculus_trf/dat.txt.rds"))

mamm_list = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/species_order.txt"))
mamm_list = as.vector(mamm_list$V1)

region = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/candidate_liftover_regions.txt"), header=T)

check_flip = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/liftover_flip_check.txt"), header=F)
check_flip = check_flip[check_flip$V4 %in% mamm_list,]
check_flip = check_flip[,c(4,5)]
colnames(check_flip) = c("species", "strand")
check_flip = rbind(check_flip, data.frame(species = "Mus_musculus", strand = "same"))

core_region = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_aware.rds"))

for(celltype in c("Hepatocytes", "Erythroid_cells")){

    values = dat_phred[, colnames(dat_phred) == celltype]
    sorted = sort(values)
    N = length(sorted)

    Q_from_x <- function(x) {
        r <- findInterval(x, sorted)
        q <- (r - 0.5) / N
        q <- pmin(q, 1 - 1e-12)
        -10 * log10(1 - q)
    }

    dat_all = list()
    cnt = 1
    for(mamm in rev(mamm_list)){
        print(mamm)
        dat_i = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_aware/", mamm, "_100bp_", celltype, ".tsv"))
        colnames(dat_i) = c("chr", "start", "end", "score")
        dat_i$phred_score = round(Q_from_x(dat_i$score), 3)
        mid = region$mid[region$species == mamm]
        strand = check_flip$strand[check_flip$species == mamm]
        filter_region = core_region[core_region$species == mamm & core_region$celltype == celltype,c("chr","core_start","core_end")]
        dat_i_sub = find_trimmed_overlaps(dat_i, filter_region)
        if(nrow(dat_i_sub) > 0){
          xxx = NULL
          for(i in 1:nrow(dat_i_sub)){
            for(j in dat_i_sub$start[i]:(dat_i_sub$end[i]-1)){
              xxx = rbind(xxx, data.frame(x_axis = j, phred_score = dat_i_sub$phred_score[i]))
            }
          }
          xxx$y_axis = cnt
          if(strand == "same"){
            xxx$x_axis = xxx$x_axis - mid
          } else {
            xxx$x_axis = mid - xxx$x_axis
          }
          dat_all[[mamm]] = xxx
        }
        cnt = cnt + 1
      }

      df = do.call(rbind, dat_all)
      df = rbind(df, data.frame(x_axis = 100000, y_axis = 1, phred_score = 0))

      p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_void() +
      theme(legend.position="none",
        panel.background = element_rect(fill = "#440154", color = NA))
      ggsave(paste0("~/share/heatmap_", celltype, "_aware.png"), p, width = 5, height = 5, dpi = 300)

      p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE)
      ggsave(paste0("~/share/heatmap_", celltype, "_aware.pdf"), p, width = 5, height = 5)
}








##################################################
### plot the big scatter plot (multispecies model)


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
library(tidyr)
library(GenomicRanges)

find_trimmed_overlaps <- function(dat_i, filter_region, return_gr = FALSE) {
  gr_A <- GRanges(seqnames = dat_i$chr, ranges = IRanges(start = dat_i$start, end = dat_i$end), score = dat_i$score, phred_score = dat_i$phred_score)
  gr_B <- GRanges(seqnames = filter_region$chr, ranges = IRanges(start = filter_region$core_start, end = filter_region$core_end))
  hits <- findOverlaps(gr_A, gr_B)
  if (length(hits) == 0) {
    if (return_gr) return(GRanges())
    return(data.frame())
  }
  overlap_gr <- pintersect(
    gr_A[queryHits(hits)],
    gr_B[subjectHits(hits)]
  )
  if (return_gr) return(overlap_gr)
  result <- data.frame(
    chr = as.character(seqnames(overlap_gr)),
    start = start(overlap_gr),
    end = end(overlap_gr),
    score = mcols(overlap_gr)$score,
    phred_score = mcols(overlap_gr)$phred_score
  )
  return(result)
}

dat_phred = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/dat.txt.rds"))

mamm_list = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/species_order.txt"))
mamm_list = as.vector(mamm_list$V1)

region = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/candidate_liftover_regions.txt"), header=T)

check_flip = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/liftover_flip_check.txt"), header=F)
check_flip = check_flip[check_flip$V4 %in% mamm_list,]
check_flip = check_flip[,c(4,5)]
colnames(check_flip) = c("species", "strand")
check_flip = rbind(check_flip, data.frame(species = "Mus_musculus", strand = "same"))

core_region = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_multispecies.rds"))

for(celltype in c("Hepatocytes", "Erythroid_cells")){

    values = dat_phred[, colnames(dat_phred) == celltype]
    sorted = sort(values)
    N = length(sorted)

    Q_from_x <- function(x) {
        r <- findInterval(x, sorted)
        q <- (r - 0.5) / N
        q <- pmin(q, 1 - 1e-12)
        -10 * log10(1 - q)
    }

    dat_all = list()
    cnt = 1
    for(mamm in rev(mamm_list)){
        print(mamm)
        dat_i = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_multispecies/", mamm, "_100bp_", celltype, ".tsv"))
        colnames(dat_i) = c("chr", "start", "end", "score")
        dat_i$phred_score = round(Q_from_x(dat_i$score), 3)
        mid = region$mid[region$species == mamm]
        strand = check_flip$strand[check_flip$species == mamm]
        filter_region = core_region[core_region$species == mamm & core_region$celltype == celltype,c("chr","core_start","core_end")]
        dat_i_sub = find_trimmed_overlaps(dat_i, filter_region)
        if(nrow(dat_i_sub) > 0){
          xxx = NULL
          for(i in 1:nrow(dat_i_sub)){
            for(j in dat_i_sub$start[i]:(dat_i_sub$end[i]-1)){
              xxx = rbind(xxx, data.frame(x_axis = j, phred_score = dat_i_sub$phred_score[i]))
            }
          }
          xxx$y_axis = cnt
          if(strand == "same"){
            xxx$x_axis = xxx$x_axis - mid
          } else {
            xxx$x_axis = mid - xxx$x_axis
          }
          dat_all[[mamm]] = xxx
        }
        cnt = cnt + 1
      }

      df = do.call(rbind, dat_all)
      df = rbind(df, data.frame(x_axis = 100000, y_axis = 1, phred_score = 0))

      p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_void() +
      theme(legend.position="none",
        panel.background = element_rect(fill = "#440154", color = NA))
      ggsave(paste0("~/share/heatmap_", celltype, "_multispecies.png"), p, width = 5, height = 5, dpi = 300)

      p = ggplot(df, aes(x_axis, y_axis, fill= phred_score)) + 
      geom_tile() +
      scale_fill_viridis(discrete=FALSE)
      ggsave(paste0("~/share/heatmap_", celltype, "_multispecies.pdf"), p, width = 5, height = 5)
}






