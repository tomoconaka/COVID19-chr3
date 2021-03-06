setwd("/home/tomoko/01.data.QC/")

library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(lubridate)
library(readr)

parse_date_covid19 <- function(date,order_type)
{
  date_out <- as_date(rep(NA,length(date)))
  for (i in 1:length(date))
  {
    i <- as.character(date[i])
    if(i=="-1" | is.na(i) | i=="1899-12-29" | i=="1899-12-30") {
      o <- NA}
    else if (nchar(i)==5) {
      o <- as_date(as.numeric(i),origin="1899-12-30") 
    }
    else if (nchar(i)>5) {
      o <- as_date(parse_date_time(i,order_type))
    }
    else {o <- NA}
    date_out[which(i==date)] <- o
  }
  return(date_out)
}

path <- "../covid19-hgi-clinical-values/"

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
italy_valenti_manifest <- italy_valenti_manifest %>% mutate(hospitalization_end = as.Date(hospitalization_end, origin = "1899-12-30"))
italy_valenti_manifest <- italy_valenti_manifest[!duplicated(italy_valenti_manifest$anonymized_patient_id),]
italy_valenti_manifest <- italy_valenti_manifest %>% mutate(death = case_when(hospitalization_end_cause == 1 ~ 1,
                                                                              hospitalization_end_cause == 0 ~ 0,
                                                                              TRUE ~ -1),
                                                            cause_of_death = case_when(death == 1 ~ 1,
                                                                                       death == 0 ~ 0,
                                                                                       TRUE ~ -1),
                                                            date_of_death = hospitalization_end)
                                                                                      
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


## hospital_info ##
a <- bel_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="belgium_migeotte", hospitalization_start=parse_date_covid19(hospitalization_start,"dmy"),hospitalization_end=parse_date_covid19(hospitalization_end,"dmy"), date_of_death=parse_date_covid19(date_of_death,"dmy"), highest_who_score=as.numeric(highest_who_score)) %>%
  mutate_at(vars(-c("anonymized_patient_id","hospitalization_start", "hospitalization_end", "date_of_death", "study")), .funs = funs(as.numeric(.)))

b <- swe_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="sweden", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

c <- spain_bujanda_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="spain_bujanda", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd"),date_of_death=parse_date_covid19(date_of_death,"ymd")) %>% mutate_at("anonymized_patient_id",as.character) 

d <- spain_butti_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% mutate(study="spain_butti", hospitalization_start=parse_date_covid19(hospitalization_start,"dmy"),hospitalization_end=as.Date(hospitalization_end),
                                     date_of_death=parse_date_covid19(date_of_death,"dmy"),highest_respiratory_support=as.numeric(highest_respiratory_support)) %>% mutate_at("anonymized_patient_id",as.character) %>%
  dplyr::select(any_of(c(identification,hospital_info,"study")))

e <- ita_renieri_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="italy_renieri")

f <- bra_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::rename(death = death., com_afib = com_afib.) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="brasil", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)
f  <- f %>% mutate(highest_who_score = case_when(death == 1 ~ 10,
                                                 highest_respiratory_support == 2 ~ 7,
                                                 highest_respiratory_support == 1 ~ 6,
                                                 highest_respiratory_support == 0 ~ 5,
                                                 hospitalization == 1 ~ 4,
                                                 TRUE ~ 1))

g <- can_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="canada", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd"),date_of_death=parse_date_covid19(date_of_death,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

h <- italy_valenti_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% mutate(study="italy_valenti", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd"), date_of_death=parse_date_covid19(date_of_death,"ymd")) %>% mutate_at("anonymized_patient_id",as.character) %>% mutate_at("icu_admit",as.numeric)


i <- germany_schulte_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% mutate(study="germany_schulte", hospitalization_start=parse_date_covid19(hospitalization_start,"ymd"),hospitalization_end=parse_date_covid19(hospitalization_end,"ymd"),
                                                                                                    date_of_death=parse_date_covid19(date_of_death,"ymd"))
j <- germany_ludwig_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="germany_ludwig",date_of_death=parse_date_covid19(date_of_death,"ymd"))

k <- spain_alarcon_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,hospital_info))) %>% 
  mutate(study="spain_alarcon")

l <- spain_gomez_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,hospital_info))) %>% mutate(study="spain_gomez") 
m <- spain_planas_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,hospital_info))) %>% mutate(study="spain_planas")

n <- belgium_rahmouni_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,hospital_info))) %>% mutate(study="belgium_rahmouni") %>% mutate_at("anonymized_patient_id",as.character)
o <- norway_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,hospital_info))) %>% mutate(study="norway") %>% mutate_at("anonymized_patient_id",as.character)
p <- italy_duga_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,hospital_info))) %>% mutate(study="italy_duga") %>% mutate_at("anonymized_patient_id",as.character)


hosp_data <- bind_rows(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p))
hosp_data <- hosp_data %>% mutate(hospitalization = ifelse(anonymized_patient_id == 24170,1,hospitalization))
table(hosp_data[,c("highest_respiratory_support", "highest_who_score")])
hosp_data <- hosp_data %>% mutate(highest_who_score = case_when(death == 1 & cause_of_death == 1 ~ 10,
                                                highest_who_score < 7 & highest_who_score != -1 & highest_respiratory_support == 2 ~ 7,
                                                highest_who_score < 6 & highest_who_score != -1 & highest_respiratory_support == 1 ~ 6,
                                                highest_who_score < 5 & highest_who_score != -1 & highest_respiratory_support == 0 ~ 5,
                                                TRUE ~ highest_who_score
))
table(hosp_data[,c("highest_respiratory_support", "highest_who_score")])
table(hosp_data[,c("hospitalization", "highest_who_score")])
#hosp_data %>% filter(highest_who_score == -1 & highest_respiratory_support == -1 & study == "spain_planas")
#spain_planas, belgium_rahmouni
#hosp_data %>% filter(highest_who_score == 4 & highest_respiratory_support != -7)
#hosp_data <- hosp_data %>% mutate(highest_who_score = case_when(death == 1 ~ 10,
#                                                                highest_who_score > 7 ~ highest_who_score,
#                                                                highest_respiratory_support == 2 ~ 7,
#                                                                highest_respiratory_support == 1 ~ 6,
#                                                                highest_respiratory_support == 0 ~ 5,
#                                                                highest_respiratory_support == -1 ~ -1,
#                                                                hospitalization == 1 ~ 4,
#                                                                highest_who_score >= 1 ~ highest_who_score,
#                                                                TRUE ~ 1
#                                                                ))



hosp_data$hospitalization_start[as.Date(hosp_data$hospitalization_start) < as.Date("2000-01-01")] <- NA
hosp_data$hospitalization_end[as.Date(hosp_data$hospitalization_end) < as.Date("2000-01-01")] <- NA
hosp_data$date_of_death[as.Date(hosp_data$date_of_death) < as.Date("2000-01-01")] <- NA
hosp_data <- hosp_data %>% mutate(temp=recode(study, belgium_migeotte="BM",sweden="SW",spain_bujanda="SBuj",spain_butti="SBut",italy_renieri="ITR",brasil="BZ",canada="CA",italy_valenti="ITV",germany_schulte="GS",germany_ludwig="GL",spain_alarcon="SA",
                                              spain_gomez="SG", spain_planas="SP", belgium_rahmouni="BR",norway="NO",italy_duga="ITD"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>%
  dplyr::select(-temp) %>% 
  mutate_at(.vars = vars(-c(anonymized_patient_id,hospitalization_start,hospitalization_end,date_of_death)), .funs = funs(ifelse(. == -1,NA,.))) %>%
  mutate(hospitalization_duration=as.numeric(hospitalization_end-hospitalization_start))

hosp_data <- hosp_data %>% mutate(icu_duration = ifelse(icu_duration > 4000, -1, icu_duration))
hosp_data %>% saveRDS(paste0("curated_clinical/hosp_data.rds"))

# Plotting - missing %
df <- hosp_data  %>% pivot_longer(c("cause_of_death", "death")) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("% missing") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

df <- hosp_data %>% pivot_longer(c("icu_duration","days_ventilator","hospitalization_duration")) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("% missing") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Exclude variables with too many missing
hosp_data_data_miss_cleaned <- hosp_data %>% purrr::discard(~sum(is.na(.x))/length(.x)* 100 >=80)

# Plotting - dates
df <- hosp_data %>% select(anonymized_patient_id,hospitalization_start,hospitalization_end,date_of_death,study) %>% pivot_longer(!c(anonymized_patient_id,study))
ggplot(df,aes(y=value,x=study)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + facet_wrap(~name) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plotting - durations
df <- hosp_data %>% select(anonymized_patient_id,hospitalization_duration,icu_duration,days_ventilator,study) %>% pivot_longer(!c(anonymized_patient_id,study))
ggplot(df,aes(y=value,x=study)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + facet_wrap(~name,scales = "free_y") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plotting - dicotom
df <- hosp_data %>% select(anonymized_patient_id,matches("hosp_"),hospitalization,icu_admit,study) %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(.==1,na.rm=T)/length(.)*100)) %>% filter(value!=0)
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("%") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
