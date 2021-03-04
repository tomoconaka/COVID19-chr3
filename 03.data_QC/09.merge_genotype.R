setwd("/home/tomoko/covid19-hgi-clinical-values/bgen")
library(stringr)
library(tidyr)
library(dplyr)
library(data.table)

#tmp0.8 <- fread("/home/tomoko/11.pca/cases_controls_pca_sup_pops_0.8_probs.txt")
tmp0.5 <- fread("/home/tomoko/11.pca/cases_controls_pca_sup_pops_0.5_probs.txt")

ibd <- fread("/home/tomoko/11.pca/ibd_all_rev2.txt")
ibd  <- ibd  %>% mutate(pihat1 = str_split(ibd, pattern = ":",simplify=TRUE)[,5])
ibd  <- ibd  %>% mutate(pihat = as.numeric(gsub( "}", "",pihat1)))
ibd  <- ibd  %>% filter(pihat > 0.1875)
#ibd  <- ibd  %>% filter(pihat >= 0.0442)

path = "../"

#canada
canada_geno <- readRDS("hgi_canada/all_chr.rds") %>% mutate(anonymized_patient_id = paste0("CA_",as.numeric(str_split(FID, "_",simplify=TRUE)[,1])))
a <- canada_geno %>% select(c(-FID))
map <- fread("/home/tomoko/11.pca/canada.sampleid.studyid.map", header=F) %>% filter(grepl("JGH",V2)) %>% 
  mutate(anonymized_patient_id = paste0("CA_",as.numeric(str_split(V2, "_",simplify=TRUE)[,2])))
a <- a %>% inner_join(map,by=c("anonymized_patient_id"="anonymized_patient_id")) %>%
  mutate(tmp=paste0(V1,"_",V1))
a <- tmp0.5 %>% merge(a, by.y="tmp", by.x="s")

#italy_valenti
italy_valenti_geno <- readRDS("hgi_italy_valenti/all_chr.rds") %>% mutate(anonymized_patient_id = paste0("ITV_",IID)) 
b <- italy_valenti_geno %>% select(c(-FID,-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
valenti_pca <- tmp0.5 %>%
  mutate(anonymized_patient_id = paste0("ITV_",s)) 
b <- b %>% merge(valenti_pca, by="anonymized_patient_id")

#spain bujanda
spain_bujanda_key <- fread(paste0("/home/tomoko/01.data.QC/hgi_spain_bujanda/linking.file")) %>% mutate(`Tube.number` = as.character(`Tube.number`))
spain_bujanda_geno <- readRDS("hgi_spain_bujanda/all_chr.rds") %>% inner_join(spain_bujanda_key,by=c("FID"="Tube.number")) %>% 
  mutate(anonymized_patient_id = paste0("SBuj_",anonymized_patient_id)) 
c <- spain_bujanda_geno %>% select(c(-FID,-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
bujanda_pca <- tmp0.5 %>%
  mutate(s = as.character(s)) %>% inner_join(spain_bujanda_key,by=c("s"="Tube.number")) %>% mutate(anonymized_patient_id = paste0("SBuj_",anonymized_patient_id))
c <- c %>% merge(bujanda_pca, by="anonymized_patient_id", all.x=T)


#spain butti
spain_butti_geno <- readRDS("hgi_spain_butti/all_chr.rds")
d <- spain_butti_geno %>% select(c(-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
butti_pca <- tmp0.5 
butti_pca$s <- as.character(butti_pca$s)
d <- d %>% merge(butti_pca, by.x="FID", by.y="s",all.x=T)
#key <- readxl::read_excel("../../covid19-hgi-clinical-values/hgi_spain_butti/ID_codification_buti_spain.xlsx")
d <- d %>% mutate(anonymized_patient_id = paste0("SBut_",FID)) %>% rename(s = FID)


tmp <- readRDS("hgi_swe_ita_bel/all_chr.rds")
tmp <- tmp %>% inner_join(tmp0.5, by=c("FID"="s"))
# Select Swedish samples
e <- tmp %>% mutate(anonymized_patient_id=as.numeric(gsub("0_COV.COV-","",FID))) %>% filter(!is.na(anonymized_patient_id)) %>%
  mutate(anonymized_patient_id=paste0("SW_",anonymized_patient_id)) %>% rename(s=FID)
# Select Italian samples
ita_renieri_manifest <- fread("/home/tomoko/01.data.QC/hgi_italy/Italy.withCom_20201209TN.tsv")
f <- tmp %>% mutate(anonymized_patient_id=paste0(gsub("_","-",gsub("^0_COV.","",FID)))) %>% 
  filter(anonymized_patient_id %in% ita_renieri_manifest$anonymized_patient_id) %>%
  mutate(anonymized_patient_id=paste0("ITR_",anonymized_patient_id)) %>% rename(s=FID)
# Select Belgium samples
bel_manifest <- read.xlsx(paste0(path,"hgi_belgium/Covid19-EGA-ErasmeData-111020_TN.xls"),sheetName = "one_visit",stringsAsFactors=FALSE)
bel_manifest <- bel_manifest[-dim(bel_manifest)[1],]

g <- tmp %>% mutate(anonymized_patient_id=paste0(as.character(gsub("0_COV.","",FID)))) %>% 
  filter(anonymized_patient_id %in% bel_manifest$anonymized_patient_id) %>%
  mutate(anonymized_patient_id=paste0("BM_",anonymized_patient_id)) %>% rename(s=FID)

belgium_rahmouni_manifest <- readRDS(paste0("../../01.data.QC/hgi_belgium_rahmouni/belgium_rahmouni_clinical.rds"))
belgium_rahmouni_manifest <- belgium_rahmouni_manifest[!duplicated(belgium_rahmouni_manifest$anonymized_patient_id),]

belgium_rahmouni_key <- read.xlsx(paste0(path,"hgi_belgium_rahmouni/BelCovid_ULiege_01152021_TN_01192020.xlsx"),sheetName = "Sample sheet IDs",stringsAsFactors=FALSE,colClasses="character",startRow = 9)

belgium_rahmouni_manifest <- belgium_rahmouni_manifest %>% inner_join(belgium_rahmouni_key, by=c("anonymized_patient_id"="X53_example"))
o <- tmp %>% mutate(tmp=paste0(as.character(gsub("0_COV.","",FID)))) %>% 
  inner_join(belgium_rahmouni_manifest[,c("anonymized_patient_id", "X124_example")], by=c("tmp"="X124_example")) %>%
  mutate(anonymized_patient_id=paste0("BR_",anonymized_patient_id)) %>% rename(s=FID)

# Select germany shuflet
germany_schulte_manifest <- readxlsb::read_xlsb(paste0(path,"hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xlsb"),sheet = 2)
germany_schulte_manifest <- germany_schulte_manifest %>% mutate(FID = paste0(as.character(gsub("TUM_","",genotyping_id))))
h <- tmp %>% mutate(anonymized_patient_id=paste0(as.character(gsub("0_COV.","",FID)))) %>% 
  filter(anonymized_patient_id %in% germany_schulte_manifest$genotyping_id) %>%
  mutate(anonymized_patient_id=paste0("GS_",anonymized_patient_id)) %>% rename(s=FID)

#brasil
brasil_geno <- readRDS("hgi_brasil/all_chr.rds")
bra_manifest <- read.xlsx(paste0(path,"hgi_brasil/chr3_dataset_01132021.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character") %>%
  mutate(anonymized_patient_id=paste0("BZ_",anonymized_patient_id))

i <- brasil_geno %>% mutate(anonymized_patient_id=paste0("BZ_",str_split(FID, "_",simplify=TRUE)[,2])) %>% 
  mutate(anonymized_patient_id=ifelse(anonymized_patient_id == "BR_CORV", paste0("BZ_",str_split(FID, "_",simplify=TRUE)[,3]),anonymized_patient_id)) %>%
  filter(anonymized_patient_id %in% bra_manifest$anonymized_patient_id) 

brasil_pca <- tmp0.5 %>% filter(s %in% brasil_geno$FID)

i <- i %>% merge(brasil_pca, by.x="FID", by.y="s") %>% rename(s=FID)

#spain_alarcon
alarcon_geno <- fread("../../covid19-hgi-clinical-values/bgen/hgi_spain_alarcon/leadSNP.raw")
j <- alarcon_geno %>% mutate(anonymized_patient_id=paste0("SA_",IID)) %>% select(c(-FID,-IID,-PAT,-MAT,-SEX,-PHENOTYPE)) 
alarcon_pca <- tmp0.5 %>% 
  mutate(anonymized_patient_id=paste0("SA_",s)) 
j <- j %>% merge(alarcon_pca, by="anonymized_patient_id")
j <- j %>% mutate(`chr3:45823240:T:C_C` = `3:45864732:T:C_C`,
                  `chr6:31153455:T:C_C` = `6:31121232:T:C_C`,
                  `chr6:41534945:A:C_C` = `6:41502683:A:C_C`,
                  `chr9:133273813:T:C_C` = `9:136149229:T:C_C`,
                  `chr12:112919637:G:A_A` = 2 - `12:113357442:G:A_G`,
                  `chr17:46143984:A:G_G` = 2 - `17:44221350:G:A_A`,
                  `chr19:4719431:G:A_A` = `19:4719443:G:A_A`,
                  `chr19:10317045:T:A_A` = `19:10427721:T:A_A`,
                  `chr21:33242905:C:T_T` = `21:34615210:T:C_T`) 

#germany ludwi
ludwig_geno <- readRDS("hgi_germany_ludwig/all_chr.rds")
ludwig_pca <- tmp0.5
ludwig_geno <- ludwig_geno %>% inner_join(ludwig_pca, by=c("FID"="s"))
k <- ludwig_geno %>% mutate(anonymized_patient_id=paste0("GL_",FID)) %>% rename(s=FID)

#spain gomez
spain_gomez_geno <- readRDS("hgi_spain_gomez/all_chr.rds") 
l <- spain_gomez_geno %>% select(c(-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
spain_gomez_pca <- tmp0.5
spain_gomez_pca$s <- as.character(spain_gomez_pca$s)
l <- l %>% merge(spain_gomez_pca, by.x="FID", by.y="s",all.x=T)
#key <- readxl::read_excel("../../covid19-hgi-clinical-values/hgi_spain_butti/ID_codification_buti_spain.xlsx")
l <- l %>% mutate(anonymized_patient_id = paste0("SG_",FID)) %>% rename(s=FID)

#spain planas
planas_geno <- fread("../../covid19-hgi-clinical-values/bgen/hgi_spain_planas/leadSNP.raw")
m <- planas_geno %>% mutate(anonymized_patient_id=paste0("SP_",IID)) %>% select(c(-FID,-IID,-PAT,-MAT,-SEX,-PHENOTYPE)) 
m <- m %>% mutate(`chr3:45823240:T:C_C` = `3:45864732:T:C_C`,
                  `chr6:31153455:T:C_C` = `6:31121232:T:C_C`,
                  `chr9:133273813:T:C_C` = `9:136149229:T:C_C`,
                  `chr12:112919637:G:A_A` = 2 - `12:113357442:G:A_G`,
                  `chr17:46143984:A:G_G` = 2 - `17:44221350:G:A_A`,
                  `chr19:4719431:G:A_A` = `19:4719443:G:A_A`,
                  `chr19:10317045:T:A_A` = `19:10427721:T:A_A`,
                  `chr21:33242905:C:T_T` = `21:34615210:T:C_T`)
planas_pca <- tmp0.5 %>% 
  mutate(anonymized_patient_id=paste0("SP_",s)) 
m <- m %>% merge(planas_pca, by="anonymized_patient_id")

#norway
norway_geno <- readRDS("../../covid19-hgi-clinical-values/bgen/hgi_norway/all_chr.rds") 
n <- norway_geno %>% mutate(anonymized_patient_id = paste0("NO_",FID)) %>% select(c(-FID,-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
norway_pca <- tmp0.5 %>% 
  mutate(anonymized_patient_id=paste0("NO_",s)) 
n <- n %>% merge(norway_pca, by="anonymized_patient_id")

#italy_duga
italy_duga_geno <- readRDS("../../covid19-hgi-clinical-values/bgen/hgi_italy_duga/all_chr.rds") 
italy_duga_geno$FID <- as.character(italy_duga_geno$FID)
p <- italy_duga_geno %>% mutate(anonymized_patient_id = paste0("ITD_",FID)) %>% select(c(-FID,-IID,-PAT,-MAT,-SEX,-PHENOTYPE))
italy_duga_pca <- tmp0.5 %>% 
  mutate(anonymized_patient_id=paste0("ITD_",s)) 
p <- p %>% merge(italy_duga_pca, by="anonymized_patient_id")

all1 <- bind_rows(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p))
#all1 <- bind_rows(list(a,b,c,d,e,f,g,h,i,j,k,m,o))
all <- all1 %>% select(c("anonymized_patient_id","chr3:45823240:T:C_C","pop",                   
                        ,"PC1","PC2","PC3","PC4","PC5", "PC6" ,"PC7","PC8", "PC9", "PC10" ,"PC11", "PC12" ,"PC13","PC14","PC15",                   
                        ,"PC16","PC17", "PC18","PC19","PC20"))

saveRDS(all,file="all_ans_pca_raw.rds")

tmp <- readRDS("~/01.data.QC/curated_clinical/complication_dat_final.rds")
tmp <- tmp[,c(1:42,90)]
genoukb <- readRDS("~/01.data.QC/ukbb/all_chr.rds")
pcukb <- tmp0.5
genoukb <- genoukb %>% merge(pcukb, by.x="FID",by.y="s")
genoukb <- genoukb %>% mutate(anonymized_patient_id = FID,
                              study = "UKB") %>% rename(s=FID) %>% mutate(s=as.character(s))
genoukb <- genoukb %>% filter(anonymized_patient_id %in% ukb$anonymized_patient_id) %>% 
  mutate(anonymized_patient_id = paste0("UKB_",anonymized_patient_id))

ALL <- bind_rows(all1, genoukb)
TMP <- tmp %>% merge(ALL, by="anonymized_patient_id")

ibd <- ibd %>% filter(i %in% TMP$s | j %in% TMP$s)
ibd_list <- unique(c(ibd$i,ibd$j))
ibd_list <- ibd_list[ibd_list %in% TMP$s]

TMP$family <- 0
TMP$flag <- 0
for(a in seq(1,length(ibd_list))){
  ibd_tmp <- ibd %>% filter(i %in% ibd_list[a] | j %in% ibd_list[a])
  ibd_tmp_list <- c(ibd_tmp$i, ibd_tmp$j)
  ibd_tmp_list <- ibd_tmp_list[ ibd_tmp_list %in% TMP$s]
  TMP$flag[TMP$s %in% ibd_tmp_list & TMP$family > 0] <- 1
  TMP$family[TMP$s %in% ibd_tmp_list & TMP$family == 0] <- a
}

TMP %>% filter(flag > 0)

TMP1 <- TMP %>% filter(family == 0)
TMP2 <- TMP %>% filter(family > 0)
TMP2 <- TMP2[!(duplicated(TMP2$family)),]

TMP3 <- bind_rows(TMP1, TMP2)

write.table(TMP3$anonymized_patient_id, file="unrelated.list", sep="\t", quote=F, col.names = F, row.names = F)
