
##############################################
### liftover between 135 species x 135 species

mamm_list=(Acinonyx_jubatus Ailuropoda_melanoleuca Ailurus_fulgens Anoura_caudifer Antilocapra_americana Aotus_nancymaae Balaenoptera_acutorostrata Bison_bison Bos_indicus Bos_mutus Bos_taurus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Cavia_aperea Cavia_porcellus Cebus_albifrons Cebus_capucinus Ceratotherium_simum Cercocebus_atys Cheirogaleus_medius Chlorocebus_sabaeus Choloepus_hoffmanni Colobus_angolensis Condylura_cristata Cricetomys_gambianus Cricetulus_griseus Ctenodactylus_gundi Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dipodomys_ordii Echinops_telfairi Elaphurus_davidianus Elephantulus_edwardii Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Felis_catus Fukomys_damarensis Gorilla_gorilla Hemitragus_hylocrius Heterocephalus_glaber Hipposideros_armiger Homo_sapiens Ictidomys_tridecemlineatus Jaculus_jaculus Leptonychotes_weddellii Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Marmota_marmota Megaderma_lyra Meriones_unguiculatus Mesocricetus_auratus Microcebus_murinus Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Monodon_monoceros Mus_caroli Mus_musculus Mus_pahari Mus_spretus Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Nannospalax_galili Neomonachus_schauinslandi Neophocaena_asiaeorientalis Nomascus_leucogenys Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Orcinus_orca Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Papio_anubis Paradoxurus_hermaphroditus Peromyscus_maniculatus Piliocolobus_tephrosceles Pithecia_pithecia Pongo_abelii Propithecus_coquereli Pteropus_alecto Pteropus_vampyrus Puma_concolor Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Saguinus_imperator Saimiri_boliviensis Sorex_araneus Spermophilus_dauricus Suricata_suricatta Sus_scrofa Tapirus_indicus Tapirus_terrestris Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tursiops_truncatus Ursus_maritimus Vicugna_pacos Vulpes_lagopus)
target_mamm_list=(Acinonyx_jubatus Ailuropoda_melanoleuca Ailurus_fulgens Anoura_caudifer Antilocapra_americana Aotus_nancymaae Balaenoptera_acutorostrata Bison_bison Bos_indicus Bos_mutus Bos_taurus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Cavia_aperea Cavia_porcellus Cebus_albifrons Cebus_capucinus Ceratotherium_simum Cercocebus_atys Cheirogaleus_medius Chlorocebus_sabaeus Choloepus_hoffmanni Colobus_angolensis Condylura_cristata Cricetomys_gambianus Cricetulus_griseus Ctenodactylus_gundi Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dipodomys_ordii Echinops_telfairi Elaphurus_davidianus Elephantulus_edwardii Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Erinaceus_europaeus Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Felis_catus Fukomys_damarensis Gorilla_gorilla Hemitragus_hylocrius Heterocephalus_glaber Hipposideros_armiger Homo_sapiens Ictidomys_tridecemlineatus Jaculus_jaculus Leptonychotes_weddellii Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Marmota_marmota Megaderma_lyra Meriones_unguiculatus Mesocricetus_auratus Microcebus_murinus Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Monodon_monoceros Mus_caroli Mus_musculus Mus_pahari Mus_spretus Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Nannospalax_galili Neomonachus_schauinslandi Neophocaena_asiaeorientalis Nomascus_leucogenys Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Orcinus_orca Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Papio_anubis Paradoxurus_hermaphroditus Peromyscus_maniculatus Piliocolobus_tephrosceles Pithecia_pithecia Pongo_abelii Propithecus_coquereli Pteropus_alecto Pteropus_vampyrus Puma_concolor Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Saguinus_imperator Saimiri_boliviensis Sorex_araneus Spermophilus_dauricus Suricata_suricatta Sus_scrofa Tapirus_indicus Tapirus_terrestris Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tursiops_truncatus Ursus_maritimus Vicugna_pacos Vulpes_lagopus)
mamm="${mamm_list[$SGE_TASK_ID - 1]}"
echo $mamm

script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install/cactus-bin-v2.9.9/bin
data_path=/net/shendure/vol8/projects/tli/ucsc_cactus/
work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15
celltype=Hepatocytes

mkdir -p "$work_path/Afp_TSS_cross_species/liftover/$mamm"

for target_mamm in "${target_mamm_list[@]}"; do
    "$script_path"/halLiftover --noDupes --bedType 4 \
    "$data_path/241-mammalian-2020v2.hal" \
    "$mamm" "$work_path/Afp_TSS/core_region_2/$mamm.$celltype.bed" \
    "$target_mamm" "$work_path/Afp_TSS_cross_species/liftover/$mamm/liftover_$target_mamm.bed"
done


##############################################################
######## After lifting over, filtering and stitching fragments

mamm_list=(Acinonyx_jubatus Ailuropoda_melanoleuca Ailurus_fulgens Anoura_caudifer Antilocapra_americana Aotus_nancymaae Balaenoptera_acutorostrata Bison_bison Bos_indicus Bos_mutus Bos_taurus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Cavia_aperea Cavia_porcellus Cebus_albifrons Cebus_capucinus Ceratotherium_simum Cercocebus_atys Cheirogaleus_medius Chlorocebus_sabaeus Choloepus_hoffmanni Colobus_angolensis Condylura_cristata Cricetomys_gambianus Cricetulus_griseus Ctenodactylus_gundi Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dipodomys_ordii Echinops_telfairi Elaphurus_davidianus Elephantulus_edwardii Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Felis_catus Fukomys_damarensis Gorilla_gorilla Hemitragus_hylocrius Heterocephalus_glaber Hipposideros_armiger Homo_sapiens Ictidomys_tridecemlineatus Jaculus_jaculus Leptonychotes_weddellii Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Marmota_marmota Megaderma_lyra Meriones_unguiculatus Mesocricetus_auratus Microcebus_murinus Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Monodon_monoceros Mus_caroli Mus_musculus Mus_pahari Mus_spretus Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Nannospalax_galili Neomonachus_schauinslandi Neophocaena_asiaeorientalis Nomascus_leucogenys Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Orcinus_orca Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Papio_anubis Paradoxurus_hermaphroditus Peromyscus_maniculatus Piliocolobus_tephrosceles Pithecia_pithecia Pongo_abelii Propithecus_coquereli Pteropus_alecto Pteropus_vampyrus Puma_concolor Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Saguinus_imperator Saimiri_boliviensis Sorex_araneus Spermophilus_dauricus Suricata_suricatta Sus_scrofa Tapirus_indicus Tapirus_terrestris Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tursiops_truncatus Ursus_maritimus Vicugna_pacos Vulpes_lagopus)
target_mamm_list=(Acinonyx_jubatus Ailuropoda_melanoleuca Ailurus_fulgens Anoura_caudifer Antilocapra_americana Aotus_nancymaae Balaenoptera_acutorostrata Bison_bison Bos_indicus Bos_mutus Bos_taurus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Cavia_aperea Cavia_porcellus Cebus_albifrons Cebus_capucinus Ceratotherium_simum Cercocebus_atys Cheirogaleus_medius Chlorocebus_sabaeus Choloepus_hoffmanni Colobus_angolensis Condylura_cristata Cricetomys_gambianus Cricetulus_griseus Ctenodactylus_gundi Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dipodomys_ordii Echinops_telfairi Elaphurus_davidianus Elephantulus_edwardii Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Erinaceus_europaeus Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Felis_catus Fukomys_damarensis Gorilla_gorilla Hemitragus_hylocrius Heterocephalus_glaber Hipposideros_armiger Homo_sapiens Ictidomys_tridecemlineatus Jaculus_jaculus Leptonychotes_weddellii Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Marmota_marmota Megaderma_lyra Meriones_unguiculatus Mesocricetus_auratus Microcebus_murinus Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Monodon_monoceros Mus_caroli Mus_musculus Mus_pahari Mus_spretus Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Nannospalax_galili Neomonachus_schauinslandi Neophocaena_asiaeorientalis Nomascus_leucogenys Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Orcinus_orca Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Papio_anubis Paradoxurus_hermaphroditus Peromyscus_maniculatus Piliocolobus_tephrosceles Pithecia_pithecia Pongo_abelii Propithecus_coquereli Pteropus_alecto Pteropus_vampyrus Puma_concolor Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Saguinus_imperator Saimiri_boliviensis Sorex_araneus Spermophilus_dauricus Suricata_suricatta Sus_scrofa Tapirus_indicus Tapirus_terrestris Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tursiops_truncatus Ursus_maritimus Vicugna_pacos Vulpes_lagopus)
mamm="${mamm_list[$SGE_TASK_ID - 1]}"
echo $mamm

script_path=/net/gs/vol1/home/cxqiu/bin/python_script
work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15
celltype=Hepatocytes

for target_mamm in "${mamm_list[@]}"; do
    Rscript "$script_path"/stitchHalFrags_v3.R \
    "$work_path/Afp_TSS/core_region_2/$mamm.$celltype.bed" \
    "$work_path/Afp_TSS_cross_species/liftover/$mamm/liftover_$target_mamm.bed" \
    "$work_path/Afp_TSS_cross_species/liftover/$mamm/stitchHalFrags_$target_mamm.bed" 0.1 2.0
done


#########################################
######## Overlap enhancers across species

mamm_list=(Acinonyx_jubatus Ailuropoda_melanoleuca Ailurus_fulgens Anoura_caudifer Antilocapra_americana Aotus_nancymaae Balaenoptera_acutorostrata Bison_bison Bos_indicus Bos_mutus Bos_taurus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Cavia_aperea Cavia_porcellus Cebus_albifrons Cebus_capucinus Ceratotherium_simum Cercocebus_atys Cheirogaleus_medius Chlorocebus_sabaeus Choloepus_hoffmanni Colobus_angolensis Condylura_cristata Cricetomys_gambianus Cricetulus_griseus Ctenodactylus_gundi Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dipodomys_ordii Echinops_telfairi Elaphurus_davidianus Elephantulus_edwardii Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Felis_catus Fukomys_damarensis Gorilla_gorilla Hemitragus_hylocrius Heterocephalus_glaber Hipposideros_armiger Homo_sapiens Ictidomys_tridecemlineatus Jaculus_jaculus Leptonychotes_weddellii Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Marmota_marmota Megaderma_lyra Meriones_unguiculatus Mesocricetus_auratus Microcebus_murinus Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Monodon_monoceros Mus_caroli Mus_musculus Mus_pahari Mus_spretus Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Nannospalax_galili Neomonachus_schauinslandi Neophocaena_asiaeorientalis Nomascus_leucogenys Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Orcinus_orca Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Papio_anubis Paradoxurus_hermaphroditus Peromyscus_maniculatus Piliocolobus_tephrosceles Pithecia_pithecia Pongo_abelii Propithecus_coquereli Pteropus_alecto Pteropus_vampyrus Puma_concolor Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Saguinus_imperator Saimiri_boliviensis Sorex_araneus Spermophilus_dauricus Suricata_suricatta Sus_scrofa Tapirus_indicus Tapirus_terrestris Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tursiops_truncatus Ursus_maritimus Vicugna_pacos Vulpes_lagopus)
mamm="${mamm_list[$SGE_TASK_ID - 1]}"
echo $mamm

work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15
celltype=Hepatocytes

mkdir -p "$work_path/Afp_TSS_cross_species/overlap/$mamm"

for target_mamm in "${mamm_list[@]}"; do
    bedtools intersect \
    -a "$work_path/Afp_TSS/core_region_2/$mamm.$celltype.bed" \
    -b "$work_path/Afp_TSS_cross_species/liftover/$target_mamm/stitchHalFrags_$mamm.bed" \
    -wa -wb > "$work_path/Afp_TSS_cross_species/overlap/$mamm/overlap_$target_mamm.bed"
done


##################################################
######## adding results together to create network

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
celltype = "Hepatocytes"
source("~/work/scripts/utils.R")

mamm_list = c("Acinonyx_jubatus","Ailuropoda_melanoleuca","Ailurus_fulgens","Anoura_caudifer","Antilocapra_americana","Aotus_nancymaae","Balaenoptera_acutorostrata","Bison_bison","Bos_indicus","Bos_mutus","Bos_taurus","Callithrix_jacchus","Camelus_bactrianus","Camelus_dromedarius","Camelus_ferus","Canis_lupus","Canis_lupus_familiaris","Capra_aegagrus","Capra_hircus","Cavia_aperea","Cavia_porcellus","Cebus_albifrons","Cebus_capucinus","Ceratotherium_simum","Cercocebus_atys","Cheirogaleus_medius","Chlorocebus_sabaeus","Choloepus_hoffmanni","Colobus_angolensis","Condylura_cristata","Cricetomys_gambianus","Cricetulus_griseus","Ctenodactylus_gundi","Dasypus_novemcinctus","Daubentonia_madagascariensis","Delphinapterus_leucas","Desmodus_rotundus","Dicerorhinus_sumatrensis","Diceros_bicornis","Dipodomys_ordii","Echinops_telfairi","Elaphurus_davidianus","Elephantulus_edwardii","Enhydra_lutris","Eptesicus_fuscus","Equus_asinus","Equus_caballus","Equus_przewalskii","Eschrichtius_robustus","Eubalaena_japonica","Eulemur_flavifrons","Felis_catus","Fukomys_damarensis","Gorilla_gorilla","Hemitragus_hylocrius","Heterocephalus_glaber","Hipposideros_armiger","Homo_sapiens","Ictidomys_tridecemlineatus","Jaculus_jaculus","Leptonychotes_weddellii","Lipotes_vexillifer","Loxodonta_africana","Lycaon_pictus","Macaca_fascicularis","Macaca_mulatta","Macaca_nemestrina","Macroglossus_sobrinus","Mandrillus_leucophaeus","Manis_javanica","Marmota_marmota","Megaderma_lyra","Meriones_unguiculatus","Mesocricetus_auratus","Microcebus_murinus","Micronycteris_hirsuta","Microtus_ochrogaster","Miniopterus_natalensis","Monodon_monoceros","Mus_caroli","Mus_musculus","Mus_pahari","Mus_spretus","Mustela_putorius","Myocastor_coypus","Myotis_brandtii","Myotis_davidii","Myotis_lucifugus","Nannospalax_galili","Neomonachus_schauinslandi","Neophocaena_asiaeorientalis","Nomascus_leucogenys","Ochotona_princeps","Octodon_degus","Odobenus_rosmarus","Odocoileus_virginianus","Orcinus_orca","Oryctolagus_cuniculus","Otolemur_garnettii","Ovis_aries","Ovis_canadensis","Pan_paniscus","Pan_troglodytes","Panthera_onca","Panthera_pardus","Panthera_tigris","Papio_anubis","Paradoxurus_hermaphroditus","Peromyscus_maniculatus","Piliocolobus_tephrosceles","Pithecia_pithecia","Pongo_abelii","Propithecus_coquereli","Pteropus_alecto","Pteropus_vampyrus","Puma_concolor","Rangifer_tarandus","Rattus_norvegicus","Rhinolophus_sinicus","Rhinopithecus_bieti","Saguinus_imperator","Saimiri_boliviensis","Sorex_araneus","Spermophilus_dauricus","Suricata_suricatta","Sus_scrofa","Tapirus_indicus","Tapirus_terrestris","Tragulus_javanicus","Trichechus_manatus","Tupaia_chinensis","Tursiops_truncatus","Ursus_maritimus","Vicugna_pacos","Vulpes_lagopus")

dat = NULL

for(mamm_i in mamm_list){
    print(mamm_i)
    for(mamm_j in mamm_list){
        if(mamm_i != mamm_j){
            file = paste0(work_path, "/14_crested/", model_id, "/Afp_TSS_cross_species/overlap/", mamm_i, "/overlap_", mamm_j, ".bed")
            if (file.info(file)$size != 0){
                dat_i = read.table(file)
                dat = rbind(dat, dat_i[,c(4,8)])
            }
        }
    }
}

df = dat
key = apply(df, 1, function(row) paste(sort(row), collapse = "_"))
edge = df[!duplicated(key), ]
node = unique(c(edge$V4, edge$V8))
node = data.frame(node_id = node,
    species = sub("_[0-9]+$", "", node))
colnames(edge) = c("node_1", "node_2")

library(igraph)
g <- graph_from_data_frame(d = edge, vertices = node, directed = FALSE)
cl <- cluster_louvain(g)
membership <- membership(cl) 

node$cluster = paste0("cluster_",membership[as.vector(node$node_id)])

target_nodes <- node$node_id[node$cluster == "cluster_1"]
subg <- induced_subgraph(g, vids = as.character(target_nodes))
sub_cl <- cluster_louvain(subg)
sub_membership <- membership(sub_cl)
node$subcluster <- NA
sub_node_ids <- names(sub_membership)
node$cluster[node$node_id %in% sub_node_ids] <- 
  paste0("cluster_1_sub_", sub_membership[as.character(node$node_id[node$node_id %in% sub_node_ids])])

tmp = node %>% group_by(cluster) %>% tally() %>% arrange(-n) %>% filter(n > 3)
tmp$cluster_update = paste0("cluster_", 1:11)
node = node %>% left_join(tmp[,(1:3)], by = "cluster")
node$cluster_update[is.na(node$cluster_update)] = "cluster_0"
node$cluster = node$cluster_update
node = node[,c("node_id", "species", "cluster")]

node_loc = NULL
for(mamm_i in mamm_list){
    print(mamm_i)
    dat_i = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_2/", mamm_i, ".", celltype, ".bed"))
    node_loc = rbind(node_loc, dat_i)
}
colnames(node_loc) = c("chr", "start", "end", "node_id")

node_add = node_loc$node_id[!node_loc$node_id %in% node$node_id]
for(x in node_add){
    node = rbind(node, data.frame(node_id = x, species = sub("_[0-9]+$", "", x), cluster = "cluster_0"))
    edge = rbind(edge, data.frame(node_1 = x, node_2 = x))
}

node_category = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/species_category.txt"))
colnames(node_category) = c("species", "species_category")

node = node %>% left_join(node_loc, by = "node_id") %>% 
    left_join(node_category, by = "species") %>% as.data.frame()

print(sum(node$cluster != 'cluster_0')/nrow(node))
### 591/621 = 95%

write.table(node, paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/node.txt"), row.names=F, col.names=T, sep="\t", quote=F)
write.table(edge, paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/edge.txt"), row.names=F, col.names=T, sep="\t", quote=F)






##########################################################
### Can we plot enhancer clusters in the big scatter plot?


work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
source("~/work/scripts/utils.R")
library(tidyr)

mamm_list = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/species_order.txt"))
mamm_list = as.vector(mamm_list$V1)

region = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/candidate_liftover_regions.txt"), header=T)

check_flip = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/liftover_flip_check.txt"), header=F)
check_flip = check_flip[check_flip$V4 %in% mamm_list,]
check_flip = check_flip[,c(4,5)]
colnames(check_flip) = c("species", "strand")
check_flip = rbind(check_flip, data.frame(species = "Mus_musculus", strand = "same"))

node = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/node.txt"), header=T)

dat_list = list()
cnt = 1
for(mamm in rev(mamm_list)){
    print(mamm)
    mid = region$mid[region$species == mamm]
    strand = check_flip$strand[check_flip$species == mamm]
    node_sub = node[node$species == mamm,]
    if(nrow(node_sub) > 0){
        for(i in 1:nrow(node_sub)){
            dat = NULL
            for(j in node_sub$start[i]:node_sub$end[i]){
                dat = rbind(dat, data.frame(x_axis = j - mid))
            }
            dat$y_axis = cnt
            dat$cluster = node_sub$cluster[i]
            if(strand != 'same'){
                dat$x_axis = 0 - dat$x_axis
            }
            dat_list[[node_sub$node_id[i]]] = dat
        }
    }
    cnt = cnt + 1
}

res = do.call(rbind, dat_list)

cluster_color_plate = c(
"cluster_0" =  "#7f7f7f",
"cluster_1" =  "#ff7f0e",
"cluster_2" =  "#2ca02c",
"cluster_3" =  "#d62728",
"cluster_4" =  "#9467bd",
"cluster_5" =  "#8c564b",
"cluster_6" =  "#e377c2",
"cluster_7" =  "#ff9896",
"cluster_8" =  "#bcbd22",
"cluster_9" =  "#17becf",
"cluster_10" = "#1f77b4",
"cluster_11" = "#ffbb78")

p = ggplot(res, aes(x_axis, y_axis, fill= cluster)) + 
      geom_tile() +
      scale_fill_manual(values=cluster_color_plate) +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_void() +
      theme(legend.position="none")
ggsave(paste0("~/share/heatmap_", celltype, "_multi_enhancer_cluster.png"), p, width = 5, height = 5, dpi = 300)


### zoom in -10k to +2k
res_sub = res[res$x_axis > (-10000) & res$x_axis < 2000,]

p = ggplot(res_sub, aes(x_axis, y_axis, fill= cluster)) + 
      geom_tile() +
      scale_fill_manual(values=cluster_color_plate) +
      scale_x_continuous(limits = c(-10000, 2000), expand = c(0, 0)) +
      scale_y_discrete(expand = c(0, 0)) +
      theme_void() +
      theme(legend.position="none")
ggsave(paste0("~/share/heatmap_", celltype, "_multi_enhancer_cluster_zoom_in.png"), p, width = 5, height = 5, dpi = 300)




#########################################################
### label each of the 9 clusters in the figure to include 
### a) number of members 
### b) number of species
### c) distance from TSS +/- stdev
### d) phred score +/- stdev
### e) mean phyloP score +/- stdev

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
celltype = "Hepatocytes"
source("~/work/scripts/utils.R")

dat = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_multispecies.rds"))
dat = as.data.frame(dat[dat$celltype == celltype,])
dat = dat[,c("chr", "core_start", "core_end", "score", "species")]
colnames(dat) = c("chr", "start", "end", "score", "species")

node = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/node.txt"), header=T)

dat = dat %>% left_join(node, by = c("chr","start","end","species"))

dat_loc = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/candidate_liftover_regions.txt"), header=T)

dat = dat %>% left_join(dat_loc[,c("species", "mid")], by = "species")

dat$distance = abs((dat$start + dat$end)/2 - dat$mid)

dat_phred = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/dat.txt.rds"))

values = dat_phred[, colnames(dat_phred) == celltype]
sorted = sort(values)
N = length(sorted)

Q_from_x <- function(x) {
    r <- findInterval(x, sorted)
    q <- (r - 0.5) / N
    q <- pmin(q, 1 - 1e-12)
    -10 * log10(1 - q)
}

dat$phred_score = Q_from_x(dat$score)

dat_sub = dat[dat$cluster != "cluster_0",]

x1 = dat_sub %>% group_by(cluster) %>% tally() %>% rename(num_enhancers = n)
x2 = dat_sub %>% group_by(cluster, species) %>% tally() %>% as.data.frame() %>%
    group_by(cluster) %>% tally() %>% rename(num_species = n)
x3 = dat_sub %>% group_by(cluster) %>% summarize(mean_dis = mean(distance), sd_dis = sd(distance))
x4 = dat_sub %>% group_by(cluster) %>% summarize(mean_phred = mean(phred_score), sd_phred = sd(phred_score))

x = x1 %>% left_join(x2, by = "cluster") %>% 
    left_join(x3, by = "cluster") %>% 
    left_join(x4, by = "cluster") 

x$dis = paste0(round(x$mean_dis/1000), "k ± ", round(x$sd_dis/1000), "k")
x$phred = paste0(round(x$mean_phred,1), " ± ", round(x$sd_phred,1))
print(x[,c("cluster", "num_enhancers", "num_species", "dis", "phred")])

for(i in 1:nrow(x)){
    print(paste0(x$cluster[i],
                 ": ", x$num_enhancers[i], " enhancers, ",
                 x$num_species[i], " species, ",
                 x$dis[i], " to TSS, ",
                 x$phred[i], " phred-score"))
}


#################################################################################
### For each of the eleven clusters, can you identify the set of syntenic orthologs 
### for ALL species (Regardless of whether an enhancer was called), and then plot 
### a histogram of the corresponding Phred-scores for each? If some are bimodal, 
### then we might actually be able to draw a conclusion

work_path = "/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq"
model_id = "mouse_fake_track_15"
celltype = "Hepatocytes"
source("~/work/scripts/utils.R")

dat_phred = readRDS(paste0(work_path, "/14_crested/", model_id, "/prediction_mammals/prediction_Mus_musculus_trf/dat.txt.rds"))

values = dat_phred[, colnames(dat_phred) == celltype]
sorted = sort(values)
N = length(sorted)

Q_from_x <- function(x) {
    r <- findInterval(x, sorted)
    q <- (r - 0.5) / N
    q <- pmin(q, 1 - 1e-12)
    -10 * log10(1 - q)
}

library(GenomicRanges)
find_trimmed_overlaps <- function(dat_i, filter_region, return_gr = FALSE) {
    gr_A <- GRanges(seqnames = dat_i$chr, ranges = IRanges(start = dat_i$start, end = dat_i$end), score = dat_i$score, phred_score = dat_i$phred_score)
    gr_B_raw <- GRanges(seqnames = filter_region$chr, ranges = IRanges(start = filter_region$start, end = filter_region$end))
    gr_B <- reduce(gr_B_raw)
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

target_mamm_list = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/species_order.txt"))[,1]
node = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/node.txt"), header=T)

node_cluster = paste0("cluster_", 1:11)

for(node_cluster_id in node_cluster){
    print(node_cluster_id)
    node_list = node$node_id[node$cluster == node_cluster_id]
    mamm_list = unique(node$species[node$cluster == node_cluster_id])

    res_all = NULL
    for(target_mamm in target_mamm_list){
        if(target_mamm %in% mamm_list){
            res = node[node$cluster == node_cluster_id &
                           node$species == target_mamm,]
            res = res[,c("chr","start","end","node_id")]
        } else {
            res = NULL
            for(mamm in mamm_list){
                file_name = paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/liftover/", mamm, "/stitchHalFrags_", target_mamm, ".bed")
                if(file.exists(file_name) && file.info(file_name)$size > 0) {
                    res_i = read.table(file_name)
                    res_i = res_i[res_i$V4 %in% node_list,]
                    res = rbind(res, res_i)
                }
            }
            colnames(res) = c("chr","start","end","node_id")
        }

        dat = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/prediction_multispecies/", target_mamm, "_100bp_", celltype, ".tsv"))
        colnames(dat) = c("chr","start","end","score")
        dat$phred_score = Q_from_x(dat$score)
        
        dat_sub = find_trimmed_overlaps(dat, res)
        if(nrow(dat_sub) > 0){
            dat_sub$species = target_mamm
            res_all = rbind(res_all, dat_sub)
        }
    }
    saveRDS(res_all, paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/syntenic_orthologs/res_", node_cluster_id, ".rds"))
}


### --- Phred score (histogram)
### --- accessibility prediction (histogram)
### --- distance from Afp TSS (log scale dot plot?)

dat_loc = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/candidate_liftover_regions.txt"), header=T)

flip_check = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/liftover_flip_check.txt"))
flip_check = unique(flip_check[,c(4,5)])
colnames(flip_check) = c("species", "strand")
flip_check = rbind(flip_check, data.frame(species = "Mus_musculus", strand = "same"))

df_all = NULL
for(node_cluster_id in node_cluster){
    df_list = list()
    print(node_cluster_id)
    dat = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/syntenic_orthologs/res_", node_cluster_id, ".rds"))
    dat = dat %>% left_join(dat_loc[,c("species", "mid")], by = "species")
    for(i in 1:nrow(dat)){
        print(paste0(i,"/",nrow(dat)))
        df_i = data.frame(loc = dat$start[i]:(dat$end[i]-1),
                          score = dat$score[i],
                          phred_score = dat$phred_score[i],
                          species = dat$species[i])
        if(flip_check$strand[flip_check$species == dat$species[i]] == "same"){
            df_i$distance = df_i$loc - dat$mid[i]
        } else {
            df_i$distance = dat$mid[i] - df_i$loc
        }
        df_list[[i]] = df_i
    }
    df = do.call(rbind, df_list)
    
    df_merged <- df %>%
        arrange(species, loc) %>%
        group_by(species) %>%
        mutate(group_id = cumsum(c(1, diff(loc) != 1))) %>% 
        group_by(species, group_id) %>%
        summarise(
            loc_start = min(loc),
            loc_end   = max(loc),
            score       = mean(score),
            phred_score = mean(phred_score),
            distance    = mean(distance),
            .groups = "drop"
        ) %>%
        select(-group_id)
    
    df_merged$cluster = node_cluster_id
    df_all = rbind(df_all, df_merged)
}

saveRDS(df_all, paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/syntenic_orthologs/df_all.rds"))


df_all$log_score = log2(df_all$score + 1)
df_all$log_phred_score = log2(df_all$phred_score + 1)

node = read.table(paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/node.txt"), header=T)
node_cluster = unique(node$cluster)
node_cluster = node_cluster[node_cluster != "cluster_0"]

df = NULL
for(node_cluster_id in node_cluster){
    df_x = df_all[df_all$cluster == node_cluster_id,]
    node_x = node[node$cluster == node_cluster_id,]
    df_x$if_enhancer_detected = if_else(df_x$species %in% node_x$species, "yes", "no")
    print(table(df_x$if_enhancer_detected))
    df = rbind(df, df_x)
}

write.table(df, paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/df.txt"), row.names=F, col.names=T, sep='\t', quote=F)



### Making the dot plot (similar to density plot but just also displaying single dots)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/mouse_fake_track_15/Afp_TSS_cross_species/df.txt", sep="\t")

# numeric sort order for clusters
df["cluster_n"] = df["cluster"].str.replace("cluster_", "").astype(int)
clusters = sorted(df["cluster_n"].unique())
cluster_to_y = {c: i for i, c in enumerate(clusters)}

# custom colors
cluster_color_plate = {
"cluster_1": "#ff7f0e",
"cluster_2": "#2ca02c",
"cluster_3": "#d62728",
"cluster_4": "#9467bd",
"cluster_5": "#8c564b",
"cluster_6": "#e377c2",
"cluster_7": "#ff9896",
"cluster_8": "#bcbd22",
"cluster_9": "#17becf",
"cluster_10": "#1f77b4",
"cluster_11": "#ffbb78"}

metrics = ["distance", "log_score", "log_phred_score"]
fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharey=True)
rng = np.random.default_rng(0)

for ax, metric in zip(axes, metrics):
    for c in clusters:
        sub = df[df["cluster_n"] == c]
        x = sub[metric].values
        y = np.full(len(x), cluster_to_y[c]) + rng.uniform(-0.35, 0.35, len(x))
        det = (sub["if_enhancer_detected"] == "yes").values
        cluster_name = f"cluster_{c}"
        col = cluster_color_plate.get(cluster_name, "grey")
        ax.scatter(x[~det], y[~det], s=14, alpha=0.5, color=col, edgecolors="none")
        ax.scatter(x[det], y[det], s=16, alpha=0.85, color=col, edgecolors="black", linewidths=0.5)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels([f"cluster_{c}" for c in clusters])
    ax.set_xlabel(metric)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    ax.set_axisbelow(True)

axes[0].set_title("Per-species ortholog scores by cluster  (n per cluster shown at right)")

for c in clusters:
    n = (df["cluster_n"] == c).sum()
    col = cluster_color_plate.get(f"cluster_{c}", "grey")
    axes[0].annotate(f"n={n}", xy=(1.01, cluster_to_y[c]),
                 xycoords=("axes fraction", "data"),
                 va="center", fontsize=8, color=col)

plt.tight_layout()
out = "/net/gs/vol1/home/cxqiu/share/df_dotplot.pdf"
plt.savefig(out, bbox_inches="tight")
print(f"wrote {out}")






