setwd("/home/tomoko/01.data.QC/")

library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lubridate)
library(readr)
library(data.table)

path <- "../covid19-hgi-clinical-values/"

## READ IN BELGIUM DATA
bel_manifest <- read.xlsx(paste0(path,"hgi_belgium/Covid19-EGA-ErasmeData-111020_TN.xls"),sheetName = "one_visit",stringsAsFactors=FALSE)
bel_manifest <- bel_manifest[-dim(bel_manifest)[1],]

belgium_rahmouni_manifest <- readRDS("hgi_belgium_rahmouni/belgium_rahmouni_clinical.rds")
belgium_rahmouni_manifest <- belgium_rahmouni_manifest[!duplicated(belgium_rahmouni_manifest$anonymized_patient_id),]

## READ IN SWEDISH DATA
swe_manifest1 <- read.xlsx(paste0(path,"hgi_swe/PronMed genetic fenotypes with WHO.xlsx"),sheetName = "Data",stringsAsFactors=FALSE)
swe_manifest2 <- read.xlsx(paste0(path,"hgi_swe/PronMed mortality etc.xlsx"),sheetName = "EXCEL_Mortality_20201216_2020-1",stringsAsFactors=FALSE)

swe_manifest2 <- swe_manifest2 %>% rename(anonymized_patient_id = Study.Subject.ID,
                                          highest_who_score = WHO.score)  
  
swe_manifest <- merge(swe_manifest1, swe_manifest2, by="anonymized_patient_id")

## READ IN SPANISH BUJANDA DATA
spain_bujanda_manifest1 <- readRDS("hgi_spain_bujanda/spain_bujanda_clinical.rds")
spain_bujanda_manifest2 <- readRDS("hgi_spain_bujanda/spain_bujanda_clinical_Basurto.rds") %>% mutate(date_of_death = as.Date(date_of_death, origin="1989-01-01"))
spain_bujanda_manifest3 <- readRDS("hgi_spain_bujanda/spain_bujanda_clinical_Cruces.rds") %>% mutate(date_of_death = as.Date(date_of_death, origin="1989-01-01"))
spain_bujanda_manifest <- bind_rows(spain_bujanda_manifest1, spain_bujanda_manifest2, spain_bujanda_manifest3)
## READ IN SPANISH BUTTI DATA
spain_butti_manifest <- read.xlsx(paste0(path,"hgi_spain_butti/Data_dictionary_HGI_Covid-19_V2_editedTN20201031.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")
spain_butti_manifest <- spain_butti_manifest %>% drop_na(anonymized_patient_id)
key <- read.xlsx(paste0(path,"hgi_spain_butti/ID_codification_buti_spain.xlsx"),sheetIndex = 1,stringsAsFactors=FALSE,colClasses="character")
key$digits3456_patient_id <- as.character(key$digits3456_patient_id)
spain_butti_manifest <- spain_butti_manifest %>% inner_join(key[,-3],by=c("anonymized_patient_id"="digits3456_patient_id")) 
spain_butti_manifest <- spain_butti_manifest %>% mutate(anonymized_patient_id = original_sample_ID) %>% dplyr::select(-original_sample_ID)

## READ IN ITALIAN DATA
ita_renieri_manifest <- fread(paste0("hgi_italy/Italy.withCom_20210119TN.tsv"))

## READ BRASIL
bra_manifest <- read.xlsx(paste0(path,"hgi_brasil/chr3_dataset_01132021_TN.xlsx"),sheetName = "Example_one_visit",stringsAsFactors=FALSE,colClasses="character")

## READ CANADA
can_manifest <- fread(paste0(path,"hgi_canada/bqc19_one_time_V2.tsv"))

## READ ITALY VALENTI
italy_valenti_manifest <- read.xlsx(paste0(path,"hgi_italy_valenti/Form_GWAS_COVID_DataDictConv_Ultima_mod_FM_RC_FM_15DIC2020.xlsx"),sheetName = "One_Time",stringsAsFactors=FALSE,colClasses="character")
italy_valenti_manifest <- italy_valenti_manifest[!duplicated(italy_valenti_manifest$anonymized_patient_id),]

## READ IN GERMAN_SCHULTE
germany_schulte_manifest <- readxlsb::read_xlsb(paste0(path,"hgi_germany_schulte/DataDictionary_Covid-19HGI_v2_COMRIMunich_Schulte_20201009.xlsb"),sheet = 2)
germany_schulte_manifest <- germany_schulte_manifest %>% mutate(anonymized_patient_id = paste0(as.character(genotyping_id)))

## READ IN GERMAN_LUDWIG
germany_ludwig_manifest <- readRDS("hgi_germany_ludwig/germany_ludwig_clinical.rds")

## READ IN SPAIN_ALARCON
spain_alarcon_manifest <- readRDS("hgi_spain_alarcon/spain_alarcon_clinical.rds")

## READ IN SPAIN_GOMEZ
spain_gomez_manifest <- readRDS("hgi_spain_gomez/spain_gomez_clinical.rds")

## READ IN SPAIN_PLANAS
spain_planas_manifest <- readRDS("hgi_spain_planas/spain_planas_clinical.rds")

## READ IN NORWAY
norway_manifest <- readRDS("hgi_norway/norway_clinical.rds")

## READ IN ITALY DUGA
italy_duga_manifest <- readRDS("hgi_italy_duga/italy_duga_clinical.rds")

## Variable name dictionary
identification <- c("anonymized_patient_id")

demographics <- c("age_at_diagnosis","sex","height","weight","smoking")

comorbidities <- 
  c("com_hiv",
    "com_immunocomp",
    "com_transplant",
    "com_autoimm_rheum",
    "com_type_i_diabetes",
    "com_type_ii_diabetes",
    "com_diabetes",
    "com_asthma",
    "com_chronic_pulm",
    "com_sleep_apnea",
    "com_liver",
    "com_gallbl",
    "com_pancreas",
    "com_chronic_kidney",
    "com_heart_failure",
    "com_hypertension",
    "com_infarction",
    "com_vascular",
    "com_stroke",
    "com_dementia",
    "com_neurological",
    "com_leukemia",
    "com_lymphoma",
    "com_malignant_solid",
    "com_dialysis",
    "com_afib",
    "com_dialysis")


hospital_info <- 
  c("hospitalization",
    "hospitalization_start",
    "hospitalization_end",
    "hospitalization_end_cause",
    "icu_admit",
    "icu_duration",
    "highest_respiratory_support",
    "days_ventilator",
    "hosp_dvt",
    "hosp_thrombo",
    "hosp_stroke",
    "hosp_infarction",
    "hosp_renal",
    "hosp_hepatic",
    "hosp_bleeding",
    "death",
    "cause_of_death",
    "date_of_death",
    "highest_who_score")


covid19 <-
  c("covid19_test",
    "covid19_test_date",
    "covid19_test_type",
    "covid19_first_symptoms_date")

medications <- 
  c("steroids",
    "biologics",
    "lmwh",
    "hcq",
    "remdesivir")

## START PROCESSING HOSPITAL INFORMATION

## Demographics ##
a <- bel_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics)))  %>% mutate_at(vars(-("anonymized_patient_id")),as.numeric) %>% mutate(study="belgium_migeotte") 
b <- swe_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(study="sweden") %>% mutate_at("anonymized_patient_id",as.character)
c <- spain_bujanda_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(weight=as.numeric(weight),study="spain_bujanda") %>% mutate_at("anonymized_patient_id",as.character)
d <- spain_butti_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(sex=as.numeric(sex),age_at_diagnosis=as.numeric(age_at_diagnosis),study="spain_butti") %>% mutate_at("anonymized_patient_id",as.character)
e <- ita_renieri_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(study="italy_renieri") 
f <- bra_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(study="brasil") %>% mutate_at("anonymized_patient_id",as.character)
g <- can_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(study="canada") %>% mutate_at("anonymized_patient_id",as.character)
h <- italy_valenti_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(study="italy_valenti") %>% mutate_at(c("sex"),as.numeric) %>%
  mutate(height = ifelse(height < 100, height*100, height))
i <- germany_schulte_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(study="germany_schulte") %>% mutate_at(c("sex"),as.numeric) %>% mutate_at("anonymized_patient_id",as.character)
j <- germany_ludwig_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(study="germany_ludwig") %>% mutate_at(c("sex"),as.numeric) %>% mutate_at("anonymized_patient_id",as.character)
k <- spain_alarcon_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,demographics))) %>% mutate(study="spain_alarcon") %>% mutate_at(c("sex"),as.numeric) %>% mutate_at("anonymized_patient_id",as.character)
l <- spain_gomez_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,demographics))) %>% mutate(study="spain_gomez") %>% mutate_at(c("sex"),as.numeric) %>% mutate_at("anonymized_patient_id",as.character)
m <- spain_planas_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,demographics))) %>% mutate(study="spain_planas") %>% mutate_at(c("sex"),as.numeric) %>% mutate_at("anonymized_patient_id",as.character)
n <- belgium_rahmouni_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,demographics))) %>% mutate(study="belgium_rahmouni") %>% mutate_at(c("sex"),as.numeric) %>% mutate_at("anonymized_patient_id",as.character)
o <- norway_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,demographics))) %>% mutate(study="norway") %>% mutate_at(c("sex"),as.numeric) %>% mutate_at("anonymized_patient_id",as.character)
p <- italy_duga_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,demographics))) %>% mutate(study="italy_duga") %>% mutate_at(c("sex"),as.numeric) %>% mutate_at("anonymized_patient_id",as.character)


# Combine
demog_data <- bind_rows(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p))

# Recode missing and fix small issues (height==1 in the swedish data and height or weight=0 in the spain_butti data)
demog_data <- demog_data %>% mutate_at(.vars = vars(-anonymized_patient_id), .funs = funs(ifelse(. == -1,NA,.))) %>% mutate(height=ifelse(height==1 | height==1.75 | height==0 | height==16 | height < 0,NA,height),weight=ifelse(weight==0 | weight < 0,NA,weight))

# Plotting - missing %
df <- demog_data %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name, ncol=3) + theme_bw() + ylab("% missing") + ylim(0,100) + theme(axis.text.x=element_text(angle=45, hjust=1)) 

# Plotting - continuos
df <- demog_data %>% select(anonymized_patient_id,age_at_diagnosis,height,weight,study) %>% pivot_longer(!c(anonymized_patient_id,study))
ggplot(df,aes(y=value,x=study)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + facet_wrap(~name,scales = "free_y") + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) 

# Plotting - dicotom
df <- demog_data %>% 
  select(anonymized_patient_id,sex,study) %>%
  pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(.==1,na.rm=T)/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name,scales = "free_y") + theme_bw() + ylab("%") + theme(axis.text.x=element_text(angle=45, hjust=1)) 

demog_data %>% mutate(temp=recode(study, belgium_migeotte="BM",sweden="SW",spain_bujanda="SBuj",spain_butti="SBut",italy_renieri="ITR",brasil="BZ",canada="CA",italy_valenti="ITV",germany_schulte="GS",germany_ludwig="GL",spain_alarcon="SA",
                                  spain_gomez="SG", spain_planas="SP", belgium_rahmouni="BR",norway="NO",italy_duga="ITD"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>% 
  dplyr::select(-temp) %>% saveRDS(paste0("curated_clinical/demog_data.rds"))
