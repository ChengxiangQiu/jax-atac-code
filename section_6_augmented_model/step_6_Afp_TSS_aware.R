

#######################################################
### Focusing on Afp TSS region predicted by aware model

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
script_path=/cactus-bin-v2.9.9/bin
data_path=XXX
work_path=XXX

"$script_path"/halLiftover --noDupes --bedType 4 \
"$data_path/241-mammalian-2020v2.hal" \
"$mamm" "$work_path/mouse_Afp_TSS.bed" \
"$target_mamm" "$work_path/liftover_$target_mamm.bed"


######## After lifting over, filtering and stitching fragments

mamm_list=(Acinonyx_jubatus Acomys_cahirinus Ailuropoda_melanoleuca Ailurus_fulgens Allactaga_bullata Alouatta_palliata Ammotragus_lervia Anoura_caudifer Antilocapra_americana Aotus_nancymaae Aplodontia_rufa Artibeus_jamaicensis Ateles_geoffroyi Balaenoptera_acutorostrata Balaenoptera_bonaerensis Beatragus_hunteri Bison_bison Bos_indicus Bos_mutus Bos_taurus Bubalus_bubalis Callicebus_donacophilus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Capromys_pilorides Carollia_perspicillata Castor_canadensis Catagonus_wagneri Cavia_aperea Cavia_porcellus Cavia_tschudii Cebus_albifrons Cebus_capucinus Ceratotherium_simum Ceratotherium_simum_cottoni Cercocebus_atys Cercopithecus_neglectus Chaetophractus_vellerosus Cheirogaleus_medius Chinchilla_lanigera Chlorocebus_sabaeus Choloepus_didactylus Choloepus_hoffmanni Chrysochloris_asiatica Colobus_angolensis Condylura_cristata Craseonycteris_thonglongyai Cricetomys_gambianus Cricetulus_griseus Crocidura_indochinensis Cryptoprocta_ferox Ctenodactylus_gundi Ctenomys_sociabilis Cuniculus_paca Dasyprocta_punctata Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dinomys_branickii Dipodomys_ordii Dipodomys_stephensi Dolichotis_patagonum Echinops_telfairi Eidolon_helvum Elaphurus_davidianus Elephantulus_edwardii Ellobius_lutescens Ellobius_talpinus Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Erinaceus_europaeus Erythrocebus_patas Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Eulemur_fulvus Felis_catus Felis_nigripes Fukomys_damarensis Galeopterus_variegatus Giraffa_tippelskirchi Glis_glis Gorilla_gorilla Graphiurus_murinus Helogale_parvula Hemitragus_hylocrius Heterocephalus_glaber Heterohyrax_brucei Hippopotamus_amphibius Hipposideros_armiger Hipposideros_galeritus Homo_sapiens Hyaena_hyaena Hydrochoerus_hydrochaeris Hystrix_cristata Ictidomys_tridecemlineatus Indri_indri Inia_geoffrensis Jaculus_jaculus Kogia_breviceps Lasiurus_borealis Lemur_catta Leptonychotes_weddellii Lepus_americanus Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Manis_pentadactyla Marmota_marmota Megaderma_lyra Mellivora_capensis Meriones_unguiculatus Mesocricetus_auratus Mesoplodon_bidens Microcebus_murinus Microgale_talazaci Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Miniopterus_schreibersii Mirounga_angustirostris Mirza_coquereli Monodon_monoceros Mormoops_blainvillei Moschus_moschiferus Mungos_mungo Murina_feae Mus_caroli Mus_pahari Mus_spretus Muscardinus_avellanarius Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Myotis_myotis Myrmecophaga_tridactyla Nannospalax_galili Nasalis_larvatus Neomonachus_schauinslandi Neophocaena_asiaeorientalis Noctilio_leporinus Nomascus_leucogenys Nycticebus_coucang Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Okapia_johnstoni Ondatra_zibethicus Onychomys_torridus Orcinus_orca Orycteropus_afer Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Pantholops_hodgsonii Papio_anubis Paradoxurus_hermaphroditus Perognathus_longimembris Peromyscus_maniculatus Petromus_typicus Phocoena_phocoena Piliocolobus_tephrosceles Pipistrellus_pipistrellus Pithecia_pithecia Platanista_gangetica Pongo_abelii Procavia_capensis Propithecus_coquereli Psammomys_obesus Pteronotus_parnellii Pteronura_brasiliensis Pteropus_alecto Pteropus_vampyrus Puma_concolor Pygathrix_nemaeus Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Rhinopithecus_roxellana Rousettus_aegyptiacus Saguinus_imperator Saiga_tatarica Saimiri_boliviensis Scalopus_aquaticus Semnopithecus_entellus Sigmodon_hispidus Solenodon_paradoxus Sorex_araneus Spermophilus_dauricus Spilogale_gracilis Suricata_suricatta Sus_scrofa Tadarida_brasiliensis Tamandua_tetradactyla Tapirus_indicus Tapirus_terrestris Thryonomys_swinderianus Tolypeutes_matacus Tonatia_saurophila Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tupaia_tana Tursiops_truncatus Uropsilus_gracilis Ursus_maritimus Vicugna_pacos Vulpes_lagopus Xerus_inauris Zalophus_californianus Zapus_hudsonius Ziphius_cavirostris)
target_mamm="${mamm_list[$SGE_TASK_ID - 1]}"
echo $target_mamm

mamm=Mus_musculus
script_path=XXX
data_path=XXX
work_path=XXX

Rscript "$script_path"/stitchHalFrags_v3.R \
"$work_path/mouse_Afp_TSS.bed" \
"$work_path/liftover_$target_mamm.bed" \
"$work_path/stitchHalFrags_$target_mamm.bed" 0.1 2.0



####################################################################
### Step-2: selecting windows which are moving forward to prediction

work_path = ""

mamm_list = read.table(paste0(work_path, "/mamm_list.txt"))
mamm_list = as.vector(mamm_list$V1)
mamm_list = mamm_list[mamm_list != "Mus_musculus"]

dat = NULL

for(mamm in mamm_list){
    print(mamm)
    file = paste0(work_path, "/stitchHalFrags_", mamm, ".bed")
    if(file.info(file)$size != 0){
        x = read.table(file)
        x$species = mamm
        dat = rbind(dat, x)
    }
}
dat = dat[,c(1,2,3,5)]
colnames(dat) = c("chr","start","end","species")
dat$region_length = dat$end - dat$start

dat = dat[dat$region_length >= 175,]
dat = rbind(dat, data.frame(chr = "chr5", start = 90490637, end = 90490837, species = "Mus_musculus", region_length = 200))

mamm_list = unique(dat$species)
### 239 species (including mouse)

res = NULL
for(mamm in mamm_list){
    print(mamm)
    chr = dat$chr[dat$species == mamm]; start = dat$start[dat$species == mamm]; end = dat$end[dat$species == mamm]
    mid = round((start + end)/2)
    chr_size = read.table(paste0(work_path, "/", mamm, "/", mamm, ".chrom.sizes.update"))
    chr_size = chr_size[chr_size[,1] == chr, 2]
    start_x = mid - 100000
    if (start_x < 1){start_x = 1}
    end_x = mid + 100000
    if (end_x > chr_size) {end_x = chr_size}
    res = rbind(res, data.frame(chr = chr, start = start_x, end = end_x, species = mamm, mid = mid))
}

write.table(res, paste0(work_path, "/candidate_liftover_regions_all.txt"), row.names=F, col.names=T, sep="\t", quote=F)
### 239 species

res_sub = res[res$mid - res$start >= 50000 & res$end - res$mid >= 50000,]
write.table(res_sub, paste0(work_path, "/candidate_liftover_regions.txt"), row.names=F, col.names=T, sep="\t", quote=F)

### 136 species (required both sides >= 50K)

### write candidate region for each species

import pandas as pd

df = pd.read_csv("candidate_liftover_regions.txt", sep=r'\s+')

for species, group in df.groupby("species"):
    rows = []
    for _, row in group.iterrows():
        for start in range(int(row["start"]), int(row["end"]) - 2114, 100):
            rows.append((row["chr"], start, start + 2114))
    pd.DataFrame(rows, columns=["chr", "start", "end"]).to_csv(
        f"prediction/{species}.txt", sep="\t", index=False, header=False
    )


###################################################################################################################################
### Step-3: creating TRF masked genome for those candidate species (only one chrs, not the whole genome for this specific analysis)

### bedtools sort -i mm10_SimpleRepeats.bed | bedtools merge -i - > mm10_SimpleRepeats.merged.bed

mamm_list=(Acinonyx_jubatus Acomys_cahirinus Ailuropoda_melanoleuca Ailurus_fulgens Allactaga_bullata Alouatta_palliata Ammotragus_lervia Anoura_caudifer Antilocapra_americana Aotus_nancymaae Aplodontia_rufa Artibeus_jamaicensis Ateles_geoffroyi Balaenoptera_acutorostrata Balaenoptera_bonaerensis Beatragus_hunteri Bison_bison Bos_indicus Bos_mutus Bos_taurus Bubalus_bubalis Callicebus_donacophilus Callithrix_jacchus Camelus_bactrianus Camelus_dromedarius Camelus_ferus Canis_lupus Canis_lupus_familiaris Capra_aegagrus Capra_hircus Capromys_pilorides Carollia_perspicillata Castor_canadensis Catagonus_wagneri Cavia_aperea Cavia_porcellus Cavia_tschudii Cebus_albifrons Cebus_capucinus Ceratotherium_simum Ceratotherium_simum_cottoni Cercocebus_atys Cercopithecus_neglectus Chaetophractus_vellerosus Cheirogaleus_medius Chinchilla_lanigera Chlorocebus_sabaeus Choloepus_didactylus Choloepus_hoffmanni Chrysochloris_asiatica Colobus_angolensis Condylura_cristata Craseonycteris_thonglongyai Cricetomys_gambianus Cricetulus_griseus Crocidura_indochinensis Cryptoprocta_ferox Ctenodactylus_gundi Ctenomys_sociabilis Cuniculus_paca Dasyprocta_punctata Dasypus_novemcinctus Daubentonia_madagascariensis Delphinapterus_leucas Desmodus_rotundus Dicerorhinus_sumatrensis Diceros_bicornis Dinomys_branickii Dipodomys_ordii Dipodomys_stephensi Dolichotis_patagonum Echinops_telfairi Eidolon_helvum Elaphurus_davidianus Elephantulus_edwardii Ellobius_lutescens Ellobius_talpinus Enhydra_lutris Eptesicus_fuscus Equus_asinus Equus_caballus Equus_przewalskii Erinaceus_europaeus Erythrocebus_patas Eschrichtius_robustus Eubalaena_japonica Eulemur_flavifrons Eulemur_fulvus Felis_catus Felis_nigripes Fukomys_damarensis Galeopterus_variegatus Giraffa_tippelskirchi Glis_glis Gorilla_gorilla Graphiurus_murinus Helogale_parvula Hemitragus_hylocrius Heterocephalus_glaber Heterohyrax_brucei Hippopotamus_amphibius Hipposideros_armiger Hipposideros_galeritus Hyaena_hyaena Hydrochoerus_hydrochaeris Hystrix_cristata Ictidomys_tridecemlineatus Indri_indri Inia_geoffrensis Jaculus_jaculus Kogia_breviceps Lasiurus_borealis Lemur_catta Leptonychotes_weddellii Lepus_americanus Lipotes_vexillifer Loxodonta_africana Lycaon_pictus Macaca_fascicularis Macaca_mulatta Macaca_nemestrina Macroglossus_sobrinus Mandrillus_leucophaeus Manis_javanica Manis_pentadactyla Marmota_marmota Megaderma_lyra Mellivora_capensis Meriones_unguiculatus Mesocricetus_auratus Mesoplodon_bidens Microcebus_murinus Microgale_talazaci Micronycteris_hirsuta Microtus_ochrogaster Miniopterus_natalensis Miniopterus_schreibersii Mirounga_angustirostris Mirza_coquereli Monodon_monoceros Mormoops_blainvillei Moschus_moschiferus Mungos_mungo Murina_feae Mus_caroli Mus_pahari Mus_spretus Muscardinus_avellanarius Mustela_putorius Myocastor_coypus Myotis_brandtii Myotis_davidii Myotis_lucifugus Myotis_myotis Myrmecophaga_tridactyla Nannospalax_galili Nasalis_larvatus Neomonachus_schauinslandi Neophocaena_asiaeorientalis Noctilio_leporinus Nomascus_leucogenys Nycticebus_coucang Ochotona_princeps Octodon_degus Odobenus_rosmarus Odocoileus_virginianus Okapia_johnstoni Ondatra_zibethicus Onychomys_torridus Orcinus_orca Orycteropus_afer Oryctolagus_cuniculus Otolemur_garnettii Ovis_aries Ovis_canadensis Pan_paniscus Pan_troglodytes Panthera_onca Panthera_pardus Panthera_tigris Pantholops_hodgsonii Papio_anubis Paradoxurus_hermaphroditus Perognathus_longimembris Peromyscus_maniculatus Petromus_typicus Phocoena_phocoena Piliocolobus_tephrosceles Pipistrellus_pipistrellus Pithecia_pithecia Platanista_gangetica Pongo_abelii Procavia_capensis Propithecus_coquereli Psammomys_obesus Pteronotus_parnellii Pteronura_brasiliensis Pteropus_alecto Pteropus_vampyrus Puma_concolor Pygathrix_nemaeus Rangifer_tarandus Rattus_norvegicus Rhinolophus_sinicus Rhinopithecus_bieti Rhinopithecus_roxellana Rousettus_aegyptiacus Saguinus_imperator Saiga_tatarica Saimiri_boliviensis Scalopus_aquaticus Semnopithecus_entellus Sigmodon_hispidus Solenodon_paradoxus Sorex_araneus Spermophilus_dauricus Spilogale_gracilis Suricata_suricatta Sus_scrofa Tadarida_brasiliensis Tamandua_tetradactyla Tapirus_indicus Tapirus_terrestris Thryonomys_swinderianus Tolypeutes_matacus Tonatia_saurophila Tragulus_javanicus Trichechus_manatus Tupaia_chinensis Tupaia_tana Tursiops_truncatus Uropsilus_gracilis Ursus_maritimus Vicugna_pacos Vulpes_lagopus Xerus_inauris Zalophus_californianus Zapus_hudsonius Ziphius_cavirostris)
data_path=/net/shendure/vol10/projects/CRE_prediction/trf_output
work_path=/net/shendure/vol2/projects/cxqiu/work/jax/atac_seq/novaseq/14_crested/genome

for mamm in "${mamm_list[@]}"; do
    echo $mamm
    input="$data_path"/"$mamm"_trf_full_genome.bed
    [ -f "$input" ] || { echo "Missing: $mamm"; continue; }
    bedtools sort -i "$input" | \
    bedtools merge -i - > "$work_path"/"$mamm"/SimpleRepeats.merged.bed
done


### creating TRF genome (only one chromosome)

import pysam
import numpy as np
import os, sys
from collections import Counter

work_path = ""
model_id = ""

mamm_list = []
mamm_chr = {}
with open(f"{work_path}/{model_id}/candidate_liftover_regions.txt") as f:
    line = f.readline()
    line = f.readline()
    while line:
        l = line.rstrip().split('\t')
        mamm_list.append(l[3])
        mamm_chr[l[3]] = l[0]
        line = f.readline()

mamm = mamm_list[int(sys.argv[1]) - 1]
chr_list = [mamm_chr[mamm]]
print(mamm)

def dinuc_freq_from_sequence(seq):
    seq = seq.upper()
    counts = Counter()
    total = 0
    for a, b in zip(seq[:-1], seq[1:]):
        if a in "ACGT" and b in "ACGT":
            counts[a + b] += 1
            total += 1
    all_dinucs = [a + b for a in "ACGT" for b in "ACGT"]
    freqs = {d: counts.get(d, 0) / total for d in all_dinucs}
    return freqs

def prepare_dinuc_probabilities(seq):
    dinuc_freq = dinuc_freq_from_sequence(seq)
    bases = np.array(['A','C','G','T'])
    mono_counts = {b:0 for b in bases}
    for d, f in dinuc_freq.items():
        mono_counts[d[0]] += f
    first_base_probs = np.array([mono_counts[b] for b in bases])
    first_base_probs /= first_base_probs.sum()
    cond_probs = {}
    for b1 in bases:
        probs = np.array([dinuc_freq[b1+b2] for b2 in bases])
        total = probs.sum()
        if total == 0:
            probs = np.ones(4)/4
        else:
            probs /= total
        cond_probs[b1] = probs
    return first_base_probs, cond_probs

def generate_dinuc_sequence(length, first_base_probs, cond_probs, seed=None):
    if seed is not None:
        np.random.seed(seed)
    bases = np.array(['A','C','G','T'])
    seq = []
    first_base = np.random.choice(bases, p=first_base_probs)
    seq.append(first_base)
    for _ in range(length-1):
        prev_base = seq[-1]
        next_base = np.random.choice(bases, p=cond_probs[prev_base])
        seq.append(next_base)
    return "".join(seq)

fa = pysam.FastaFile(f"{work_path}/genome/{mamm}/{mamm}.fa")

seqs = {name: list(fa.fetch(name)) for name in fa.references if name in chr_list}

seqs_length = {name: len(seqs[name]) for name in fa.references if name in chr_list}

folder_path = f"{work_path}/genome/{mamm}/fasta_trf_Afp/"
os.makedirs(folder_path, exist_ok=True)

for iter in range(1,11):
    print(iter)
    with open(f"{work_path}/genome/{mamm}/SimpleRepeats.merged.bed") as bed:
        for line in bed:
            chrom, start, end = line.strip().split('\t')[:3]
            if chrom in chr_list:
                start, end = int(start), int(end)
                length = end - start
                middle = (start + end)//2
                seq_local = fa.fetch(chrom, max(1, middle - 2500), min(seqs_length[chrom], middle + 2500))
                first_base_probs, cond_probs = prepare_dinuc_probabilities(seq_local)
                rand_bases = generate_dinuc_sequence(length, first_base_probs, cond_probs)
                seqs[chrom][start:end] = rand_bases
    with open(f"{folder_path}/{mamm}_trf_{iter}.fa", "w") as out:
        for chrom in chr_list:
            out.write(f">{chrom}\n")
            seq = "".join(seqs[chrom])
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")




##############################
### Step-4: perform prediction

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

model_path = f"evoluation-aware.keras"
model = keras.models.load_model(model_path, compile=False)

celltype_list = ["Adipocyte_cells","Adipocyte_cells_Cyp2e1","B_cells","Brain_capillary_endothelial_cells","CNS_neurons","Cardiomyocytes","Corticofugal_neurons","Endocardial_cells","Endothelium","Epithelial_cells","Erythroid_cells","Eye","Glia","Glomerular_endothelial_cells","Gut_epithelial_cells","Hepatocytes","Intermediate_neuronal_progenitors","Kidney","Lateral_plate_and_intermediate_mesoderm","Liver_sinusoidal_endothelial_cells","Lung_and_airway","Lymphatic_vessel_endothelial_cells","Melanocyte_cells","Mesoderm","Neural_crest_PNS_neurons","Neuroectoderm_and_glia","Olfactory_ensheathing_cells","Olfactory_neurons","Oligodendrocytes","Skeletal_muscle_cells","T_cells","White_blood_cells"]

adata = ad.read_h5ad(f"{work_path}/data_window_cluster_top3K.h5ad")

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
        print(celltype)
        celltype_id = np.where(adata.obs.index == celltype)[0][0]
        scores = defaultdict(list)
        for i, (chrom, start, end) in enumerate(regions):
            for pos in range(start, end, 100):
                scores[(chrom, pos)].append(predictions_mean[i, celltype_id])
        rows = []
        for (chrom, pos), vals in sorted(scores.items()):
            rows.append((chrom, pos, pos + 100, np.mean(vals)))
        df = pd.DataFrame(rows, columns=["chr", "start", "end", "score"])
        df.to_csv(f"{work_path}/{model_id}/Afp_TSS/prediction_aware/{mamm}_100bp_{celltype}.tsv", sep="\t", index=False, header=False)


######################################################
### Step-5: counting how many elements in each species


library(tidyr)
library(GenomicRanges)
library(dplyr)

dat_phred = readRDS(paste0(work_path, "/14_crested/mouse_fake_track_14/prediction_mammals/prediction_Mus_musculus_trf/dat.txt.rds"))

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
        dat_i = read.table(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_aware/", mamm, "_100bp_", celltype, ".tsv"))
        colnames(dat_i) = c("chr", "start", "end", "score")
        dat_i$phred_score = round(Q_from_x(dat_i$score), 3)
        dat_i$species = mamm
        dat_i$celltype = celltype
        dat = rbind(dat, dat_i[,c("chr", "start", "end", "phred_score", "species", "celltype")])
    }
}

saveRDS(dat, paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_aware_res.rds"))






###################################################
### Step-6: Identify core region for Afp TSS region 


library(tidyr)
library(GenomicRanges)
library(dplyr)
options(scipen = 999)

dat = readRDS(paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/prediction_aware_res.rds"))

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
        res = result[result$species == mamm,c(2:4)]
        write.table(res, paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/core_region_aware/", mamm, ".", celltype, ".bed"), row.names=F, col.names=F, sep="\t", quote=F)
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

model_path = f"evoluation-aware.keras"
model = keras.models.load_model(model_path, compile=False)

adata = ad.read_h5ad(f"{work_path}/mouse_fake_track_14/data_window_cluster_top3K.h5ad")

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
    if os.path.exists(f"{work_path}/{model_id}/Afp_TSS/core_region_aware/{mamm}.{celltype}.bed"):
        fasta = pysam.FastaFile(f"{work_path}/genome/{mamm}/{mamm}.fa")
        chr_size = {}
        with open(f"{work_path}/genome/{mamm}/{mamm}.chrom.sizes.update") as f:
            for line in f:
                c, s = line.rstrip().split('\t')
                chr_size[c] = int(s)
        all_regions = []
        with open(f"{work_path}/{model_id}/Afp_TSS/core_region_aware/{mamm}.{celltype}.bed") as f:
            for line in f:
                l = line.rstrip().split('\t')
                loc = (int(l[1]) + int(l[2])) // 2
                if loc-1057 > 0 and loc+1057 < chr_size[l[0]]:
                    all_regions.append((l[0], int(l[1]), int(l[2])))
        result = []
        for chrom, s, e in all_regions:
            loc = (s + e) // 2
            local_seq = fasta.fetch(chrom, max(1, loc - 2500), min(chr_size[chrom], loc + 2500))
            if all(x == 'N' for x in local_seq):
                continue
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
        with open(f"{work_path}/{model_id}/Afp_TSS/core_region_aware/{mamm}.{celltype}.core_region.txt", "w") as out:
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
        file_name = paste0(work_path, "/14_crested/", model_id, "/Afp_TSS/core_region_aware/", mamm, ".", celltype, ".core_region.txt")
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
        }
    }
}

saveRDS(dat_all, paste0(work_path, "/14_crested/mouse_fake_track_15/Afp_TSS/core_region_aware.rds"))









