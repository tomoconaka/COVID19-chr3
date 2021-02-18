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
germany_schulte_manifest$covid19_first_symptoms_date[germany_schulte_manifest$anonymized_patient_id == "TUM_018"] <- as.Date("2020-03-07")

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

covid19 <-
  c("covid19_test",
    "covid19_test_date",
    "covid19_test_type",
    "covid19_first_symptoms_date")


## COVID19 info ##
a <- bel_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="belgium_migeotte",  covid19_test_date=parse_date_covid19(covid19_test_date,"dmy"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"dmy")) %>%
  mutate(covid19_test = as.numeric(covid19_test), covid19_test_type = as.numeric(covid19_test_type))

b <- swe_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="sweden", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

c <- spain_bujanda_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="spain_bujanda", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

d <- spain_butti_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="spain_butti", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

e <- ita_renieri_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="italy_renieri", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd"))

f <- bra_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id)%>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="brasil", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

g <- can_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="canada", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd")) %>% mutate_at("anonymized_patient_id",as.character)

h <- italy_valenti_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id)  %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="italy_valenti", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd"))

i <- germany_schulte_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="germany_schulte", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd"))

j <- germany_ludwig_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% 
  mutate(study="germany_ludwig", covid19_test_date=parse_date_covid19(covid19_test_date,"ymd"),covid19_first_symptoms_date=parse_date_covid19(covid19_first_symptoms_date,"ymd"))

k <- spain_alarcon_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% dplyr::select(any_of(c(identification,covid19))) %>% mutate(study="spain_alarcon")

l <- spain_gomez_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,covid19))) %>% mutate(study="spain_gomez") 

m <- spain_planas_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,covid19))) %>% mutate(study="spain_planas")

n <- belgium_rahmouni_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,covid19))) %>% mutate(study="belgium_rahmouni") %>% mutate_at("anonymized_patient_id",as.character)
o <- norway_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,covid19))) %>% mutate(study="norway") %>% mutate_at("anonymized_patient_id",as.character)
p <- italy_duga_manifest %>% filter(covid19_test == 1) %>% drop_na(anonymized_patient_id) %>% select(any_of(c(identification,covid19))) %>% mutate(study="italy_duga") %>% mutate_at("anonymized_patient_id",as.character)

covid_data <- bind_rows(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p))

covid_data <- covid_data %>% 
  mutate_at(.vars = vars(-c(anonymized_patient_id,covid19_first_symptoms_date,covid19_test_date)), .funs = funs(ifelse(. == -1,NA,.)))

covid_data$covid19_test_date[as.Date(covid_data$covid19_test_date) < as.Date("2000-01-01")] <- NA
covid_data$covid19_first_symptoms_date[as.Date(covid_data$covid19_first_symptoms_date) < as.Date("2000-01-01")] <- NA

covid_data %>% mutate(time=as.numeric(covid19_test_date-covid19_first_symptoms_date)) %>% filter(time < 0)


#%>% filter(covid19_first_symptoms_date < "2020-08-31" | is.na(covid19_first_symptoms_date))

covid_data <- covid_data %>% mutate(temp=recode(study, belgium_migeotte="BM",sweden="SW",spain_bujanda="SBuj",spain_butti="SBut",italy_renieri="ITR",brasil="BZ",canada="CA",italy_valenti="ITV",germany_schulte="GS",germany_ludwig="GL",spain_alarcon="SA",
                                                spain_gomez="SG", spain_planas="SP", belgium_rahmouni="BR",norway="NO",italy_duga="ITD"),anonymized_patient_id=paste0(temp,"_",anonymized_patient_id)) %>%
  dplyr::select(-temp) %>% saveRDS(paste0("curated_clinical/covid_data.rds"))

covid_data <- bind_rows(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o))

covid_data <- covid_data %>% 
  mutate_at(.vars = vars(-c(anonymized_patient_id,covid19_first_symptoms_date,covid19_test_date)), .funs = funs(ifelse(. == -1,NA,.)))

# Plotting - missing %
df <- covid_data %>% select(anonymized_patient_id,covid19_first_symptoms_date,covid19_test_date,study) %>% pivot_longer(!c(anonymized_patient_id,study)) 
df <- df %>% group_by(study,name) %>% summarise_all(funs(sum(is.na(.))/length(.)*100))
ggplot(df,aes(x=study,y=value)) + geom_bar(stat="identity") + facet_wrap(~name) + theme_bw() + ylab("% missing") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plotting - dates
df <- covid_data %>% select(anonymized_patient_id,covid19_first_symptoms_date,covid19_test_date,study) %>% pivot_longer(!c(anonymized_patient_id,study))
df$value[as.Date(df$value) < as.Date("2000-12-29")] <- NA
ggplot(df,aes(y=value,x=study)) + geom_violin(trim=FALSE, fill="gray") + geom_boxplot(width=0.1) + facet_wrap(~name) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

df <- covid_data %>% select(anonymized_patient_id,covid19_first_symptoms_date,covid19_test_date,study) 
df$covid19_first_symptoms_date[as.Date(df$covid19_first_symptoms_date) < as.Date("2000-12-29")] <- NA
