
####################################################################
### Lifting over mouse candidate 100-bp windows to other 241 mammals

### Support data can be downloaded from:
### https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download/

### Please contact Chengxiang (CX) Qiu for any questions!
### cxqiu@uw.edu or chengxiang.qiu@dartmouth.edu


#################################################################
### Step-1: Lifting over those candidate windows to other mammals

mamm_list=(Acinonyx_jubatus Acomys_cahirinus Ailuropoda_melanoleuca Ailurus_fulgens Allactaga_bullata Alouatta_palliata Ammotragus_lervia Anoura_caudifer Antilocapra_americana Aotus_nancymaae Aplodontia_rufa Artibeus_jamaicensis Ateles_geoffroyi Balaenoptera_acutorostrata Balaenoptera_bonaerensis Beatragus_hunteri Bison_bison Bos_indicus Bos_mutus Bos_taurus Bubalus_bubalis Callicebus_donacophilus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Capromys_pilorides Carollia_perspicillata Castor_canadensis Catagonus_wagneri Cavia_aperea Cavia_porcellus Cavia_tschudii Cebus_albifrons Cebus_capucinus Ceratotherium_simum Ceratotherium_simum_cottoni Cercocebus_atys Cercopithecus_neglectus Chaetophractus_vellerosus Cheirogaleus_medius Chinchilla_lanigera Chlorocebus_sabaeus Choloepus_didactylus Choloepus_hoffmanni Chrysochloris_asiatica Colobus_angolensis Condylura_cristata Craseonycteris_thonglongyai Cricetomys_gambianus Cricetulus_griseus Crocidura_indochinensis Cryptoprocta_ferox Ctenodactylus_gundi Ctenomys_sociabilis Cuniculus_paca Dasyprocta_punctata Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dinomys_branickii Dipodomys_ordii Dipodomys_stephensi Dolichotis_patagonum Echinops_telfairi Eidolon_helvum Elaphurus_davidianus Elephantulus_edwardii Ellobius_lutescens Ellobius_talpinus Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Erinaceus_europaeus Erythrocebus_patas Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Eulemur_fulvus Felis_catus Felis_nigripes Fukomys_damarensis Galeopterus_variegatus Giraffa_tippelskirchi Glis_glis Gorilla_gorilla Graphiurus_murinus Helogale_parvula Hemitragus_hylocrius Heterocephalus_glaber Heterohyrax_brucei Hippopotamus_amphibius Hipposideros_armiger Hipposideros_galeritus Homo_sapiens Hyaena_hyaena Hydrochoerus_hydrochaeris Hystrix_cristata Ictidomys_tridecemlineatus Indri_indri Inia_geoffrensis Jaculus_jaculus Kogia_breviceps Lasiurus_borealis Lemur_catta Leptonychotes_weddellii Lepus_americanus Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Manis_pentadactyla Marmota_marmota Megaderma_lyra Mellivora_capensis Meriones_unguiculatus Mesocricetus_auratus Mesoplodon_bidens Microcebus_murinus Microgale_talazaci Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Miniopterus_schreibersii Mirounga_angustirostris Mirza_coquereli Monodon_monoceros Mormoops_blainvillei Moschus_moschiferus Mungos_mungo Murina_feae Mus_caroli Mus_pahari Mus_spretus Muscardinus_avellanarius Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Myotis_myotis Myrmecophaga_tridactyla Nannospalax_galili Nasalis_larvatus Neomonachus_schauinslandi Neophocaena_asiaeorientalis Noctilio_leporinus Nomascus_leucogenys Nycticebus_coucang Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Okapia_johnstoni Ondatra_zibethicus Onychomys_torridus Orcinus_orca Orycteropus_afer Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Pantholops_hodgsonii Papio_anubis Paradoxurus_hermaphroditus Perognathus_longimembris Peromyscus_maniculatus Petromus_typicus Phocoena_phocoena Piliocolobus_tephrosceles Pipistrellus_pipistrellus Pithecia_pithecia Platanista_gangetica Pongo_abelii Procavia_capensis Propithecus_coquereli Psammomys_obesus Pteronotus_parnellii Pteronura_brasiliensis Pteropus_alecto Pteropus_vampyrus Puma_concolor Pygathrix_nemaeus Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Rhinopithecus_roxellana Rousettus_aegyptiacus Saguinus_imperator Saiga_tatarica Saimiri_boliviensis Scalopus_aquaticus Semnopithecus_entellus Sigmodon_hispidus Solenodon_paradoxus Sorex_araneus Spermophilus_dauricus Spilogale_gracilis Suricata_suricatta Sus_scrofa Tadarida_brasiliensis Tamandua_tetradactyla Tapirus_indicus Tapirus_terrestris Thryonomys_swinderianus Tolypeutes_matacus Tonatia_saurophila Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tupaia_tana Tursiops_truncatus Uropsilus_gracilis Ursus_maritimus Vicugna_pacos Vulpes_lagopus Xerus_inauris Zalophus_californianus Zapus_hudsonius Ziphius_cavirostris)
target_mamm="${mamm_list[$SGE_TASK_ID - 1]}"
echo $target_mamm

mamm=Mus_musculus
script_path=/net/shendure/vol10/projects/cxqiu/nobackup/install/cactus-bin-v2.9.9/bin
data_path=/net/shendure/vol8/projects/tli/ucsc_cactus/
work_path="/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_${mamm}"

mkdir -p "$work_path/liftover"

"$script_path"/halLiftover --noDupes --bedType 4 \
"$data_path/241-mammalian-2020v2.hal" \
"$mamm" "$work_path/window_list_merged.txt" \
"$target_mamm" "$work_path/liftover/liftover_$target_mamm.bed"


######## After lifting over, filtering and stitching fragments

mamm_list=(Acinonyx_jubatus Acomys_cahirinus Ailuropoda_melanoleuca Ailurus_fulgens Allactaga_bullata Alouatta_palliata Ammotragus_lervia Anoura_caudifer Antilocapra_americana Aotus_nancymaae Aplodontia_rufa Artibeus_jamaicensis Ateles_geoffroyi Balaenoptera_acutorostrata Balaenoptera_bonaerensis Beatragus_hunteri Bison_bison Bos_indicus Bos_mutus Bos_taurus Bubalus_bubalis Callicebus_donacophilus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Capromys_pilorides Carollia_perspicillata Castor_canadensis Catagonus_wagneri Cavia_aperea Cavia_porcellus Cavia_tschudii Cebus_albifrons Cebus_capucinus Ceratotherium_simum Ceratotherium_simum_cottoni Cercocebus_atys Cercopithecus_neglectus Chaetophractus_vellerosus Cheirogaleus_medius Chinchilla_lanigera Chlorocebus_sabaeus Choloepus_didactylus Choloepus_hoffmanni Chrysochloris_asiatica Colobus_angolensis Condylura_cristata Craseonycteris_thonglongyai Cricetomys_gambianus Cricetulus_griseus Crocidura_indochinensis Cryptoprocta_ferox Ctenodactylus_gundi Ctenomys_sociabilis Cuniculus_paca Dasyprocta_punctata Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dinomys_branickii Dipodomys_ordii Dipodomys_stephensi Dolichotis_patagonum Echinops_telfairi Eidolon_helvum Elaphurus_davidianus Elephantulus_edwardii Ellobius_lutescens Ellobius_talpinus Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Erinaceus_europaeus Erythrocebus_patas Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Eulemur_fulvus Felis_catus Felis_nigripes Fukomys_damarensis Galeopterus_variegatus Giraffa_tippelskirchi Glis_glis Gorilla_gorilla Graphiurus_murinus Helogale_parvula Hemitragus_hylocrius Heterocephalus_glaber Heterohyrax_brucei Hippopotamus_amphibius Hipposideros_armiger Hipposideros_galeritus Homo_sapiens Hyaena_hyaena Hydrochoerus_hydrochaeris Hystrix_cristata Ictidomys_tridecemlineatus Indri_indri Inia_geoffrensis Jaculus_jaculus Kogia_breviceps Lasiurus_borealis Lemur_catta Leptonychotes_weddellii Lepus_americanus Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Manis_pentadactyla Marmota_marmota Megaderma_lyra Mellivora_capensis Meriones_unguiculatus Mesocricetus_auratus Mesoplodon_bidens Microcebus_murinus Microgale_talazaci Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Miniopterus_schreibersii Mirounga_angustirostris Mirza_coquereli Monodon_monoceros Mormoops_blainvillei Moschus_moschiferus Mungos_mungo Murina_feae Mus_caroli Mus_pahari Mus_spretus Muscardinus_avellanarius Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Myotis_myotis Myrmecophaga_tridactyla Nannospalax_galili Nasalis_larvatus Neomonachus_schauinslandi Neophocaena_asiaeorientalis Noctilio_leporinus Nomascus_leucogenys Nycticebus_coucang Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Okapia_johnstoni Ondatra_zibethicus Onychomys_torridus Orcinus_orca Orycteropus_afer Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Pantholops_hodgsonii Papio_anubis Paradoxurus_hermaphroditus Perognathus_longimembris Peromyscus_maniculatus Petromus_typicus Phocoena_phocoena Piliocolobus_tephrosceles Pipistrellus_pipistrellus Pithecia_pithecia Platanista_gangetica Pongo_abelii Procavia_capensis Propithecus_coquereli Psammomys_obesus Pteronotus_parnellii Pteronura_brasiliensis Pteropus_alecto Pteropus_vampyrus Puma_concolor Pygathrix_nemaeus Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Rhinopithecus_roxellana Rousettus_aegyptiacus Saguinus_imperator Saiga_tatarica Saimiri_boliviensis Scalopus_aquaticus Semnopithecus_entellus Sigmodon_hispidus Solenodon_paradoxus Sorex_araneus Spermophilus_dauricus Spilogale_gracilis Suricata_suricatta Sus_scrofa Tadarida_brasiliensis Tamandua_tetradactyla Tapirus_indicus Tapirus_terrestris Thryonomys_swinderianus Tolypeutes_matacus Tonatia_saurophila Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tupaia_tana Tursiops_truncatus Uropsilus_gracilis Ursus_maritimus Vicugna_pacos Vulpes_lagopus Xerus_inauris Zalophus_californianus Zapus_hudsonius Ziphius_cavirostris)
target_mamm="${mamm_list[$SGE_TASK_ID - 1]}"
echo $target_mamm

mamm=Mus_musculus
script_path=/net/gs/vol1/home/cxqiu/bin/python_script
data_path=/net/shendure/vol8/projects/tli/ucsc_cactus/
work_path="/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_${mamm}"

Rscript "$script_path"/stitchHalFrags_v3.R \
"$work_path/window_list_merged.txt" \
"$work_path/liftover/liftover_$target_mamm.bed" \
"$work_path/liftover/stitchHalFrags_$target_mamm.bed" 0.5 2


#####################################################################################
### Step-2: summarizing for each window, how many mammals are consistently identified

import os

work_path = ""
mamm = "Mus_musculus"

mamm_list = []
with open(f"{work_path}/14_crested/genome/mamm_list.txt") as file:
    for line in file:
    if line.strip() != mamm:
        mamm_list.append(line.strip())

window_count = {}
for cnt, mamm_i in enumerate(mamm_list, 1):
    print(f"{cnt}/{len(mamm_list)}: {mamm_i}")
    tmp = set()
    bed_path = f"{work_path}/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_{mamm}/liftover/stitchHalFrags_{mamm_i}.bed"
    with open(bed_path) as file:
        for line in file:
            fields = line.rstrip().split('\t')
            tmp.add(fields[3])
    for i in tmp:
        window_count[i] = window_count.get(i, 0) + 1

output_path = f"{work_path}/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_{mamm}/window_list_merged.replicate_count.txt"
with open(output_path, "w") as out:
    for i in window_count:
        out.write(f"{i}\t{window_count[i]}\n")


>>> shell >>>
bedtools intersect -a window_list_merged.txt -b ../../../TSS_2500.bed -wa | uniq > window_list_merged.overlap_promoter.txt

bedtools intersect -a window_list_merged.txt -b ../../../TSS_2500_human.bed -wa | uniq > window_list_merged.overlap_promoter.txt

>>> R >>>

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

mamm = "Mus_musculus"

peak_list_sub = read.table(paste0(work_path, 
                                  "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/window_list_merged.txt"))
colnames(peak_list_sub) = c("chr", "start", "end", "window_id")

peak_list_sub_count = read.table(paste0(work_path, 
                                  "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/window_list_merged.replicate_count.txt"))
colnames(peak_list_sub_count) = c("window_id", "count")

peak_list_sub_promoter = read.table(paste0(work_path, 
                                  "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/window_list_merged.overlap_promoter.txt"))
colnames(peak_list_sub_promoter) = c("chr", "start", "end", "window_id")

df = peak_list_sub %>% left_join(peak_list_sub_count, by = "window_id") %>% as.data.frame()
df$count[is.na(df$count)] = 0
df$promoter = if_else(df$window_id %in% peak_list_sub_promoter$window_id, "yes", "no")

p = ggplot() +
    geom_histogram(data = df, aes(count), bins = 50) +
    theme_classic(base_size = 10) +
    geom_vline(xintercept = 120) 
ggsave(paste0(mamm, "_windows_replicated_in_other_mammals.pdf"), p, width = 6, height = 4.5)

p1 = ggplot() +
    geom_histogram(data = subset(df, promoter == 'yes'), aes(count), bins = 50) +
    labs(title = "Overlapped with promoter")
p2 = ggplot() +
    geom_histogram(data = subset(df, promoter != 'yes'), aes(count), bins = 50) +
    labs(title = "Non-overlapped with promoter")
ggsave(paste0(mamm, "_windows_replicated_in_other_mammals_by_promoter.pdf"), p1+p2, width = 10, height = 5)


celltype_list = read.table("celltype_id.txt")
celltype_list = as.vector(celltype_list$V1)

peak_sig = NULL
for(i in 1:36){
    print(i)
    tmp = read.table(paste0(work_path, 
                            "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/celltype_", i, "_enhancer.bed"))
    tmp = tmp[,c(1:3)]
    colnames(tmp) = c("chr", "start", "end")
    tmp$celltype = celltype_list[i]
    peak_sig = rbind(peak_sig, tmp)
}

df$x = paste0(df$chr, "_", df$start, "_", df$end)
peak_sig$x = paste0(peak_sig$chr, "_", peak_sig$start, "_", peak_sig$end)
df_x = df %>% left_join(peak_sig[,c("x", "celltype")], by = "x")

p = ggplot() +
    geom_histogram(data = df_x, aes(count), bins = 50) +
    facet_wrap(~celltype, nrow = 6, ncol = 6, scales = "free")
    #facet_wrap(~celltype, nrow = 6, ncol = 6)
ggsave(paste0(mamm, "_windows_replicated_in_other_mammals_by_celltype.pdf"), p, width = 10, height = 10)



################################################################################
### Step-3: Subset 100 bp windows which have been succ liftover to > 120 mammals

### Mus_musculus
### before filtering: n = 2,790,521
### after filtering: n = 1,465,915

work_path = ""
web_path = "https://shendure-web.gs.washington.edu/content/members/cxqiu/public/backup/jax_atac/download"
source("help_code/utils.R")

celltype_list = read.table("celltype_id.txt")
celltype_list = as.vector(celltype_list$V1)

mamm = "Mus_musculus"

dat = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_", mamm, "/dat.txt.rds"))
peak_list = readRDS(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/prediction_", mamm, "/dat_loc.txt.rds"))
peak_list$window_id = paste0(peak_list$chr, "_", peak_list$start, "_", peak_list$end)

peak_list_sub = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/window_list_merged.txt"))
colnames(peak_list_sub) = c("chr", "start", "end", "window_id")
peak_list_sub_count = read.table(paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/window_list_merged.replicate_count.txt"))
colnames(peak_list_sub_count) = c("window_id", "count")
df = peak_list_sub %>% left_join(peak_list_sub_count, by = "window_id") %>% as.data.frame()
df$count[is.na(df$count)] = 0
df_sub = df[df$count > 120,]
df_sub$window_id = paste0(df_sub$chr, "_", df_sub$start, "_", df_sub$end)

keep = peak_list$window_id %in% df_sub$window_id

write.table(dat[keep,], paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/dat_sub.txt"),
            row.names=F, col.names=F, sep="\t", quote=F)
write.table(peak_list[keep, c(1:3)], paste0(work_path, "/14_crested/celltype_L2_cut_norm/prediction_mammals/tmp_", mamm, "/peak_list_sub.txt"),
            row.names=F, col.names=F, sep="\t", quote=F)



