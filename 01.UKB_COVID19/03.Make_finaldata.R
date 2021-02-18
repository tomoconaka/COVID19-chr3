setwd("/project/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/06.clinical_values")

data <- readRDS("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/UKBall20210120.457941.rds")
test <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/covid19_result_20210120.txt.gz")
table(data$status)
#     0      1 
#444989  12952 
#hesin
hesin_diag <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/hesin_diag_20210122.txt.gz", sep="\t")
hesin <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/hesin_20210122.txt.gz", sep="\t")
hesin1 <- hesin[grepl("/2020$", hesin$epiend),]
#neg <- unique(data$ID[data$status == 0])
#tmp <- hesin1[hesin1$eid %in% neg,]
#hesin_new_neg <- merge(tmp, hesin_diag, by=c("eid", "ins_index"))
#hesin_new_neg <- hesin_new_neg %>% filter(diag_icd10 %in% c("U071","U072"))
#data <- data %>% mutate(status = ifelse(ID %in% unique(hesin_new_neg$eid), 1, status),
#                        COVIDdate = ifelse(ID %in% unique(hesin_new_neg$eid), NA, COVIDdate))
#table(data$status)
#0      1 
#450583   8597
data$COVIDdate <- as.Date(data$COVIDdate, origin="1970-01-01")
data <- data %>% filter(status == 1)
data <- data %>% mutate(covid19_test_type = case_when(!is.na(COVIDdate) ~ 1,
                                                      TRUE ~ -1))
#for(i in seq(1,length(unique(hesin_new_neg$eid)))){
#  tmp <- hesin_new_neg %>% filter(eid == unique(hesin_new_neg$eid)[i])
#  data$COVIDdate[data$ID == unique(hesin_new_neg$eid)[i]] <- min(as.Date(tmp$epistart, format="%d/%m/%Y"))
#}

identification <- c("anonymized_patient_id")
demographics <- c("age_at_diagnosis", "sex", "ancestry", "height", "weight", "smoking")
comorbidities <- c("com_hiv","com_diabetes",
                   "com_asthma","com_immunocomp","com_transplant","com_dementia",
                   "com_chronic_pulm",
                   "com_liver",
                   "com_chronic_kidney",
                   "com_dialysis",
                   "com_heart_failure",
                   "com_hypertension",
                   "com_stroke",
                   "com_afib",
                   "com_leukemia",
                   "com_lymphoma",
                   "com_malignant_solid")

hospital <- c("hospitalization",
              "hospitalization_start",
              "hospitalization_end",
              "death",
              "cause_of_death",
              "date_of_death",
              "icu_admit",
              "icu_duration",
              "highest_who_score",
              "highest_respiratory_support",
              "days_ventilator",
              "hosp_dvt",
              "hosp_thrombo",
              "hosp_stroke",
              "hosp_infarction",
              "hosp_renal",
              "hosp_bleeding",
              "hosp_hepatic")

covid19 <- c("covid19_test",
             "covid19_test_date",
             "covid19_test_type",
             "covid19_first_symptoms_date")

rx <- c("steroids",
        "biologics",
        "lmwh",
        "hcq",
        "remdesivir")

library(dplyr)
library(tidyr)

data <- data %>% rename(anonymized_patient_id = ID, age_at_diagnosis = AGE) %>% mutate(sex = ifelse(SEX == 0, 1, 0))
#demographics
#"height",50 "weight",21002 "smoking" 1239.0.0 1249.0.0	2644.0.0
infile <- "scratch/ukb27449_20688_fetch_demographics.tab.gz"
t <- fread(infile, sep="\t", na.strings = c("NA"))
t <- t %>% rename(ID = f.eid) %>% mutate(height = case_when(!is.na(f.50.0.0) ~ f.50.0.0,
                                                            !is.na(f.50.1.0) ~ f.50.1.0,
                                                            !is.na(f.50.2.0) ~ f.50.2.0,
                                                            TRUE ~ -1),
                                         weight = case_when(!is.na(f.21002.0.0) ~ f.21002.0.0,
                                                            !is.na(f.21002.1.0) ~ f.21002.1.0,
                                                            !is.na(f.21002.2.0) ~ f.21002.2.0,
                                                            TRUE ~ -1),
                                         smoking = case_when(f.1239.0.0 == 1 ~ 0,
                                                             f.1239.0.0 == 2 ~ 0, 
                                                             f.1239.1.0 == 1 ~ 0,
                                                             f.1239.1.0 == 2 ~ 0,
                                                             f.1239.2.0 == 1 ~ 0,
                                                             f.1239.2.0 == 2 ~ 0,
                                                             f.1239.0.0 %in% c(0,-3) & f.1249.0.0 == 1 ~ 1,
                                                             f.1239.0.0 == 0 & f.1249.0.0 %in% c(2,3) & f.2644.0.0 == 1 ~ 1,
                                                             f.1239.0.0 == -3 & f.1249.0.0 == 2 & f.2644.0.0 == 1 ~ 1,
                                                             f.1239.1.0 %in% c(0,-3) & f.1249.1.0 == 1 ~ 1,
                                                             f.1239.1.0 == 0 & f.1249.1.0 %in% c(2,3) & f.2644.1.0 == 1 ~ 1,
                                                             f.1239.1.0 == -3 & f.1249.1.0 == 2 & f.2644.1.0 == 1 ~ 1,
                                                             f.1239.2.0 %in% c(0,-3) & f.1249.2.0 == 1 ~ 1,
                                                             f.1239.2.0 == 0 & f.1249.2.0 %in% c(2,3) & f.2644.2.0 == 1 ~ 1,
                                                             f.1239.2.0 == -3 & f.1249.2.0 == 2 & f.2644.2.0 == 1 ~ 1,
                                                             f.1239.0.0 == 0 & f.1249.0.0 == 4 ~ 2,
                                                             f.1239.0.0 == 0 & f.1249.0.0 %in% c(2,3) & f.2644.0.0 == 0 ~ 2,
                                                             f.1239.1.0 == 0 & f.1249.1.0 == 4 ~ 2,
                                                             f.1239.1.0 == 0 & f.1249.1.0 %in% c(2,3) & f.2644.1.0 == 0 ~ 2,
                                                             f.1239.2.0 == 0 & f.1249.2.0 == 4 ~ 2,
                                                             f.1239.2.0 == 0 & f.1249.2.0 %in% c(2,3) & f.2644.2.0 == 0 ~ 2,
                                                             TRUE ~ -1
                                                             )
)
t <- t %>% select(c("ID", "height","weight","smoking"))
data <- data %>% merge(t, by.x="anonymized_patient_id", by.y="ID", all.x=T)

#comorbidities 6152 20002
infile <- "scratch/ukb27449_20688_fetch_comorbidities.tab.gz"
t <- fread(infile, sep="\t", na.strings = c("NA"))
colnames(t)[1] <- "ID"

#20002
f <- names(t)[names(t) %like% "f.20002."]
t1 <- t[, f.20002.x.x := do.call(paste0,list(c(.SD, "NA"), sep=",", collapse="")), by=ID, .SDcols=f]
f.20002 <- t1[,c("ID", "f.20002.x.x")]
data <- merge(data, f.20002, by.x="anonymized_patient_id", by.y="ID", all.x=T, sort=F)
data <- data %>% mutate(com_hiv = case_when(grepl("1439", f.20002.x.x) ~ 1,
                                            TRUE ~ 0),
                        com_dementia = case_when(grepl("1263", f.20002.x.x) ~ 1,
                                                 TRUE ~ 0),
                        com_diabetes = case_when(grepl("1220", f.20002.x.x) ~ 1,
                                                 grepl("1222", f.20002.x.x) ~ 1,
                                                 grepl("1223", f.20002.x.x) ~ 1,
                                                 TRUE ~ 0),
                        com_asthma = case_when(grepl("1111", f.20002.x.x) ~ 1,
                                                TRUE ~ 0),
                        com_chronic_pulm = case_when(grepl("1112", f.20002.x.x) ~ 1,
                                                     grepl("1113", f.20002.x.x) ~ 1,
                                                   TRUE ~ 0),
                        com_liver = case_when(grepl("1158", f.20002.x.x) ~ 1,
                                              grepl("1604", f.20002.x.x) ~ 1,
                                              TRUE ~ 0),
                        com_chronic_kidney = case_when(grepl("1192", f.20002.x.x) ~ 1,
                                                       grepl("1193", f.20002.x.x) ~ 1,
                                                       grepl("1194", f.20002.x.x) ~ 1,
                                                       TRUE ~ 0),
                        com_heart_failure = case_when(grepl("1076", f.20002.x.x) ~ 1,
                                                      TRUE ~ 0),
                        com_hypertension = case_when(grepl("1065", f.20002.x.x) ~ 1,
                                                     grepl("1072", f.20002.x.x) ~ 1,
                                                     TRUE ~ 0),
                        com_stroke = case_when(grepl("1583", f.20002.x.x) ~ 1,
                                               grepl("1081", f.20002.x.x) ~ 1,
                                               TRUE ~ 0),
                        com_afib = case_when(grepl("1471", f.20002.x.x) ~ 1,
                                               TRUE ~ 0)
                        ) %>% select(c(-f.20002.x.x))

##6152
f <- names(t)[names(t) %like% "f.6152."]
t1 <- t[, f.6152.x.x := do.call(paste0,list(c(.SD, "NA"), sep=",", collapse="")), by=ID, .SDcols=f]
f.6152 <- t1[,c("ID", "f.6152.x.x")]
data <- merge(data, f.6152, by.x="anonymized_patient_id", by.y="ID", all.x=T, sort=F)
data <- data %>% mutate(com_asthma = case_when(grepl("8", f.6152.x.x) ~ 1,
                                               TRUE ~ com_asthma),
                        com_chronic_pulm = case_when(grepl("6", f.6152.x.x) ~ 1,
                                                     TRUE ~ com_chronic_pulm)
) %>% select(c(-f.6152.x.x))
## Docter diagnosed 
#asthma 22127
#DM 2443
#cancer 2453
#COPD 22130 22129 22128

f.22127 <- t[,c("ID", "f.22127.0.0")]
data <- merge(data, f.22127, by.x="anonymized_patient_id", by.y="ID", all.x=T, sort=F)
data <- data %>% mutate(com_asthma = case_when(grepl(1, f.22127.0.0) ~ 1,
                                               TRUE ~ com_asthma)
) %>% select(c(-f.22127.0.0))

f.2443 <- t[,c("ID", "f.2443.0.0")]
data <- merge(data, f.2443, by.x="anonymized_patient_id", by.y="ID", all.x=T, sort=F)
data <- data %>% mutate(com_diabetes = case_when(grepl("1", f.2443.0.0) ~ 1,
                                               TRUE ~ com_diabetes)
) %>% select(c(-f.2443.0.0))

f.2453 <- t[,c("ID", "f.2453.0.0")]
data <- merge(data, f.2453, by.x="anonymized_patient_id", by.y="ID", all.x=T, sort=F)
data <- data %>% mutate(com_cancer = case_when(grepl("1", f.2453.0.0) ~ 1,
                                                 TRUE ~ 0)
) %>% select(c(-f.2453.0.0))

f.22130 <- t[,c("ID", "f.22130.0.0")]
data <- merge(data, f.22130, by.x="anonymized_patient_id", by.y="ID", all.x=T, sort=F)
data <- data %>% mutate(com_chronic_pulm = case_when(grepl("1", f.22130.0.0) ~ 1,
                                               TRUE ~ com_chronic_pulm)
) %>% select(c(-f.22130.0.0))

f.22129 <- t[,c("ID", "f.22129.0.0")]
data <- merge(data, f.22129, by.x="anonymized_patient_id", by.y="ID", all.x=T, sort=F)
data <- data %>% mutate(com_chronic_pulm = case_when(grepl("1", f.22129.0.0) ~ 1,
                                                     TRUE ~ com_chronic_pulm)
) %>% select(c(-f.22129.0.0))
f.22128 <- t[,c("ID", "f.22128.0.0")]
data <- merge(data, f.22128, by.x="anonymized_patient_id", by.y="ID", all.x=T, sort=F)
data <- data %>% mutate(com_chronic_pulm = case_when(grepl("1", f.22128.0.0) ~ 1,
                                                     TRUE ~ com_chronic_pulm)
) %>% select(c(-f.22128.0.0))

# Diagnoses - main ICD9
f <- names(t)[names(t) %like% "f.41203."]
t1 <- t[, f.41203.x.x := do.call(paste0,list(c(.SD, "NA"), sep=",", collapse="")), by=ID, .SDcols=f]
f.41203 <- t1[,c("ID", "f.41203.x.x")]
data <- merge(data, f.41203, by.x="anonymized_patient_id", by.y="ID", all.x=T)
data <- data %>% mutate(com_transplant = 0,
                        com_immunocomp = 0,
                        com_dementia = 0)
data <- data %>% mutate(com_hiv = case_when(grepl("27917", f.41203.x.x) ~ 1,
                                            TRUE ~ com_hiv),
                        com_transplant = case_when(grepl("^9968", f.41203.x.x) ~ 1,
                                                   grepl("^V42", f.41203.x.x) ~ 1,
                                                   TRUE ~ com_transplant),
                        com_immunocomp = case_when(grepl("^9631", f.41203.x.x) ~ 1,
                                                   grepl("^E9331", f.41203.x.x) ~ 1,
                                                   grepl("^279", f.41203.x.x) ~ 1,
                                                   TRUE ~ com_immunocomp),
                        com_dementia = case_when(grepl("^290", f.41203.x.x) ~ 1,
                                                 grepl("^2912", f.41203.x.x) ~ 1,
                                                 grepl("^2941", f.41203.x.x) ~ 1,
                                                 TRUE ~ com_dementia),
                        com_diabetes = case_when(grepl("^250", f.41203.x.x) ~ 1,
                                                 grepl("^3572", f.41203.x.x) ~ 1,
                                                 grepl("^5881", f.41203.x.x) ~ 1,
                                                 TRUE ~ com_diabetes),
                        com_asthma = case_when(grepl("^493", f.41203.x.x) ~ 1,
                                               TRUE ~ com_asthma),
                        com_chronic_pulm = case_when(grepl("^4910", f.41203.x.x) ~ 1,
                                                     grepl("^4911", f.41203.x.x) ~ 1,
                                                     grepl("^4912", f.41203.x.x) ~ 1,
                                                     grepl("^4918", f.41203.x.x) ~ 1,
                                                     grepl("^4919", f.41203.x.x) ~ 1,
                                                     grepl("^4929", f.41203.x.x) ~ 1,
                                                     grepl("^496", f.41203.x.x) ~ 1,
                                                     grepl("^492", f.41203.x.x) ~ 1,
                                                     TRUE ~ com_chronic_pulm),
                        com_liver = case_when(grepl("^5710", f.41203.x.x) ~ 1,
                                              grepl("^5712", f.41203.x.x) ~ 1,
                                              grepl("^5713", f.41203.x.x) ~ 1,
                                              grepl("^5714", f.41203.x.x) ~ 1,
                                              grepl("^5715", f.41203.x.x) ~ 1,
                                              grepl("^5716", f.41203.x.x) ~ 1,
                                              grepl("^5718", f.41203.x.x) ~ 1,
                                              grepl("^5719", f.41203.x.x) ~ 1,
                                              grepl("^5721", f.41203.x.x) ~ 1,
                                              grepl("^5722", f.41203.x.x) ~ 1,
                                              grepl("^5723", f.41203.x.x) ~ 1,
                                              grepl("^5724", f.41203.x.x) ~ 1,
                                              grepl("^5728", f.41203.x.x) ~ 1,
                                              grepl("^573", f.41203.x.x) ~ 1,
                                              TRUE ~ com_liver),
                        com_chronic_kidney = case_when(grepl("^585", f.41203.x.x) ~ 1,
                                                       grepl("^586", f.41203.x.x) ~ 1,
                                                       TRUE ~ com_chronic_kidney),
                        com_heart_failure = case_when(grepl("^428", f.41203.x.x) ~ 1,
                                                      TRUE ~ com_heart_failure),
                        com_hypertension = case_when(grepl("^401", f.41203.x.x) ~ 1,
                                                     grepl("^405", f.41203.x.x) ~ 1,
                                                     TRUE ~ com_hypertension),
                        com_afib = case_when(grepl("4273", f.41203.x.x) ~ 1,
                                             TRUE ~ com_afib),
                        com_cancer = case_when(grepl("^14", f.41203.x.x) ~ 1,
                                               grepl("^15", f.41203.x.x) ~ 1,
                                               grepl("^16", f.41203.x.x) ~ 1,
                                               grepl("^17", f.41203.x.x) ~ 1,
                                               grepl("^18", f.41203.x.x) ~ 1,
                                               grepl("^19", f.41203.x.x) ~ 1,
                                               grepl("^20", f.41203.x.x) ~ 1,
                                               grepl("^23", f.41203.x.x) ~ 1,
                                               TRUE ~ com_cancer),
) %>% select(c(-f.41203.x.x))

# Diagnoses - secondary ICD9
f <- names(t)[names(t) %like% "f.41205."]
t1 <- t[, f.41205.x.x := do.call(paste0,list(c(.SD, "NA"), sep=",", collapse="")), by=ID, .SDcols=f]
f.41205 <- t1[,c("ID", "f.41205.x.x")]
data <- merge(data, f.41205, by.x="anonymized_patient_id", by.y="ID", all.x=T)
data <- data %>% mutate(com_hiv = case_when(grepl("27917", f.41205.x.x) ~ 1,
                                            TRUE ~ com_hiv),
                        com_transplant = case_when(grepl("^9968", f.41205.x.x) ~ 1,
                                                   grepl("^V42", f.41205.x.x) ~ 1,
                                                   TRUE ~ com_transplant),
                        com_immunocomp = case_when(grepl("^9631", f.41205.x.x) ~ 1,
                                                   grepl("^E9331", f.41205.x.x) ~ 1,
                                                   grepl("^279", f.41205.x.x) ~ 1,
                                                   TRUE ~ com_immunocomp),
                        com_dementia = case_when(grepl("^290", f.41205.x.x) ~ 1,
                                                 grepl("^2912", f.41205.x.x) ~ 1,
                                                 grepl("^2941", f.41205.x.x) ~ 1,
                                                 TRUE ~ com_dementia),
                        com_diabetes = case_when(grepl("^250", f.41205.x.x) ~ 1,
                                                 grepl("^3572", f.41205.x.x) ~ 1,
                                                 grepl("^5881", f.41205.x.x) ~ 1,
                                                 TRUE ~ com_diabetes),
                        com_asthma = case_when(grepl("^493", f.41205.x.x) ~ 1,
                                               TRUE ~ com_asthma),
                        com_chronic_pulm = case_when(grepl("^4910", f.41205.x.x) ~ 1,
                                                     grepl("^4911", f.41205.x.x) ~ 1,
                                                     grepl("^4912", f.41205.x.x) ~ 1,
                                                     grepl("^4918", f.41205.x.x) ~ 1,
                                                     grepl("^4919", f.41205.x.x) ~ 1,
                                                     grepl("^4929", f.41205.x.x) ~ 1,
                                                     grepl("^496", f.41205.x.x) ~ 1,
                                                     grepl("^492", f.41205.x.x) ~ 1,
                                                     TRUE ~ com_chronic_pulm),
                        com_liver = case_when(grepl("^5710", f.41205.x.x) ~ 1,
                                              grepl("^5712", f.41205.x.x) ~ 1,
                                              grepl("^5713", f.41205.x.x) ~ 1,
                                              grepl("^5714", f.41205.x.x) ~ 1,
                                              grepl("^5715", f.41205.x.x) ~ 1,
                                              grepl("^5716", f.41205.x.x) ~ 1,
                                              grepl("^5718", f.41205.x.x) ~ 1,
                                              grepl("^5719", f.41205.x.x) ~ 1,
                                              grepl("^5721", f.41205.x.x) ~ 1,
                                              grepl("^5722", f.41205.x.x) ~ 1,
                                              grepl("^5723", f.41205.x.x) ~ 1,
                                              grepl("^5724", f.41205.x.x) ~ 1,
                                              grepl("^5728", f.41205.x.x) ~ 1,
                                              grepl("^573", f.41205.x.x) ~ 1,
                                              TRUE ~ com_liver),
                        com_chronic_kidney = case_when(grepl("^585", f.41205.x.x) ~ 1,
                                                       grepl("^586", f.41205.x.x) ~ 1,
                                                       TRUE ~ com_chronic_kidney),
                        com_heart_failure = case_when(grepl("^428", f.41205.x.x) ~ 1,
                                                      TRUE ~ com_heart_failure),
                        com_hypertension = case_when(grepl("^401", f.41205.x.x) ~ 1,
                                                     grepl("^405", f.41205.x.x) ~ 1,
                                                     TRUE ~ com_hypertension),
                        com_afib = case_when(grepl("4273", f.41205.x.x) ~ 1,
                                             TRUE ~ com_afib),
                        com_cancer = case_when(grepl("^14", f.41205.x.x) ~ 1,
                                               grepl("^15", f.41205.x.x) ~ 1,
                                               grepl("^16", f.41205.x.x) ~ 1,
                                               grepl("^17", f.41205.x.x) ~ 1,
                                               grepl("^18", f.41205.x.x) ~ 1,
                                               grepl("^19", f.41205.x.x) ~ 1,
                                               grepl("^20", f.41205.x.x) ~ 1,
                                               TRUE ~ com_cancer),
) %>% select(c(-f.41205.x.x))

##Hesin

for(i in seq(1,dim(data)[1])){
  tmp1 <- hesin %>% filter(eid == data$anonymized_patient_id[i]) %>%
    filter(as.Date(epistart, "%d/%m/%Y") < as.Date(data$COVIDdate[i]))
  tmp <- hesin_diag %>% filter(eid == data$anonymized_patient_id[i]) %>% filter(ins_index %in% tmp1$ins_index)
  if(dim(tmp)[1]>0){
  icdlist <- unique(tmp$diag_icd10)
  data$com_hiv[i][any(grepl("^B20", icdlist))] <- 1
  data$com_hiv[i][any(grepl("^B21", icdlist))] <- 1
  data$com_hiv[i][any(grepl("^B22", icdlist))] <- 1
  data$com_hiv[i][any(grepl("^B23", icdlist))] <- 1
  data$com_hiv[i][any(grepl("^B24", icdlist))] <- 1
  data$com_transplant[i][any(grepl("^T86", icdlist))] <- 1
  data$com_transplant[i][any(grepl("^N165", icdlist))] <- 1
  data$com_transplant[i][any(grepl("^Y480", icdlist))] <- 1
  data$com_transplant[i][any(grepl("^Z94", icdlist))] <- 1
  data$com_immunocomp[i][any(grepl("^D80", icdlist))] <- 1
  data$com_immunocomp[i][any(grepl("^D81", icdlist))] <- 1
  data$com_immunocomp[i][any(grepl("^D82", icdlist))] <- 1
  data$com_immunocomp[i][any(grepl("^D83", icdlist))] <- 1
  data$com_immunocomp[i][any(grepl("^D84", icdlist))] <- 1
  data$com_immunocomp[i][any(grepl("^D89", icdlist))] <- 1
  data$com_immunocomp[i][any(grepl("^T451", icdlist))] <- 1
  data$com_immunocomp[i][any(grepl("^Y434", icdlist))] <- 1
  data$com_dementia[i][any(grepl("^F00", icdlist))] <- 1
  data$com_dementia[i][any(grepl("^F01", icdlist))] <- 1
  data$com_dementia[i][any(grepl("^F02", icdlist))] <- 1
  data$com_dementia[i][any(grepl("^F03", icdlist))] <- 1
  data$com_dementia[i][any(grepl("^F05", icdlist))] <- 1
  data$com_diabetes[i][any(grepl("^E10", icdlist))] <- 1
  data$com_diabetes[i][any(grepl("^E11", icdlist))] <- 1
  data$com_diabetes[i][any(grepl("^E12", icdlist))] <- 1
  data$com_diabetes[i][any(grepl("^E13", icdlist))] <- 1
  data$com_diabetes[i][any(grepl("^E14", icdlist))] <- 1
  data$com_asthma[i][any(grepl("^J45", icdlist))] <- 1
  data$com_asthma[i][any(grepl("^J46", icdlist))] <- 1
  data$com_chronic_pulm[i][any(grepl("^J41", icdlist))] <- 1
  data$com_chronic_pulm[i][any(grepl("^J42", icdlist))] <- 1
  data$com_chronic_pulm[i][any(grepl("^J43", icdlist))] <- 1
  data$com_chronic_pulm[i][any(grepl("^J44", icdlist))] <- 1
  data$com_liver[i][any(grepl("^K7", icdlist))] <- 1
  data$com_chronic_kidney[i][any(grepl("^N18", icdlist))] <- 1
  data$com_chronic_kidney[i][any(grepl("^N19", icdlist))] <- 1
  data$com_chronic_kidney[i][any(grepl("^I120", icdlist))] <- 1
  data$com_chronic_kidney[i][any(grepl("^I131", icdlist))] <- 1
  data$com_heart_failure[i][any(grepl("^I110", icdlist))] <- 1
  data$com_heart_failure[i][any(grepl("^I119", icdlist))] <- 1
  data$com_heart_failure[i][any(grepl("^I130", icdlist))] <- 1
  data$com_heart_failure[i][any(grepl("^I132", icdlist))] <- 1
  data$com_heart_failure[i][any(grepl("^I150", icdlist))] <- 1
  data$com_hypertension[i][any(grepl("^I1", icdlist))] <- 1
  data$com_stroke[i][any(grepl("^I60", icdlist))] <- 1
  data$com_stroke[i][any(grepl("^I61", icdlist))] <- 1
  data$com_stroke[i][any(grepl("^I62", icdlist))] <- 1
  data$com_stroke[i][any(grepl("^I63", icdlist))] <- 1
  data$com_stroke[i][any(grepl("^I64", icdlist))] <- 1
  data$com_stroke[i][any(grepl("^I65", icdlist))] <- 1
  data$com_stroke[i][any(grepl("^I66", icdlist))] <- 1
  data$com_afib[i][any(grepl("^I48", icdlist))] <- 1
  data$com_cancer[i][any(grepl("^C", icdlist))] <- 1    
  }
}

saveRDS(data, file="tmp.rds")

data <- readRDS("tmp.rds")
#hospitalization
hesin_diag <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/hesin_diag_20210122.txt.gz", sep="\t")
hesin <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/hesin_20210122.txt.gz", sep="\t")
hesin1 <- hesin[grepl("/2020$", hesin$epiend),]

tmp <- hesin1[hesin1$eid == unique(data$anonymized_patient_id)[1],]
TMP <- tmp
for(i in seq(2,length(unique(data$anonymized_patient_id)))){
  tmp <- hesin1[hesin1$eid == unique(data$anonymized_patient_id)[i],]
  TMP <- rbind(TMP,tmp)
}
hesin_new <- merge(TMP, hesin_diag, by=c("eid", "ins_index"))
hesin_oper <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/hesin_oper_20210122.txt.gz", sep="\t")
hesin_new1 <- merge(TMP, hesin_oper, by=c("eid", "ins_index"))

data$hospitalization <- data$hospital
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  tmp1 <- tmp %>% filter(as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$COVIDdate) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$COVIDdate)  + 30 )
  tmp1 <- tmp1 %>%  filter(as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$COVIDdate) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$COVIDdate))
  if(dim(tmp1)[1] > 0){
    data$hospitalization[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
    data$hospitalization_start[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- as.Date(tmp1$epistart[1], "%d/%m/%Y")
    data$hospitalization_end[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- as.Date(tmp1$epiend[1], "%d/%m/%Y")
  }
  }
table(data$hospitalization)
#   0    1 
# 7056 1541 
hesin_covid <- hesin_new %>% filter(diag_icd10 %in% c("U071", "U072"))
for(i in seq(1,length(unique(hesin_covid$eid)))){
  tmp <- hesin_covid[hesin_covid$eid == unique(hesin_covid$eid)[i],]
  if(dim(tmp)[1] > 0){
    data$hospitalization[data$anonymized_patient_id == unique(hesin_covid$eid)[i]] <- 1
    data$hospitalization_start[data$anonymized_patient_id == unique(hesin_covid$eid)[i]] <- as.Date(tmp$epistart[1], "%d/%m/%Y")
    data$hospitalization_end[data$anonymized_patient_id == unique(hesin_covid$eid)[i]] <- as.Date(tmp$epiend[1], "%d/%m/%Y")
    }  
  }
}
table(data$hospitalization)
as.Date(tmp$epistart[1], "%d/%m/%Y")
data$hospitalization_start <- as.Date(data$hospitalization_start, origin="1970-01-01")
data$hospitalization_end <- as.Date(data$hospitalization_end, origin="1970-01-01")

data <- data %>% mutate(COVIDdate = case_when(!is.na(hospitalization_start) & as.Date(COVIDdate) > as.Date(hospitalization_start) ~ as.Date(hospitalization_start),
                                              TRUE ~ as.Date(COVIDdate)))
                        
data %>% filter(anonymized_patient_id == 1783711) %>% str()

table(data$hospitalization)

#   0    1 
#6899 1698
data <- data %>% rename(covid19_test_date = COVIDdate, covid19_test = status)
data <- data %>% select(-date_of_death)

d <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/death_20210121.txt.gz", sep="\t")
unique(d$date_of_death[grepl("11\\/2020$", d$date_of_death)])
d1 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/death_cause_20210121.txt.gz", sep="\t")
D <- merge(d, d1, by=c("eid", "ins_index"))
colnames(D)[1] <- "ID"
D1 <- D[as.Date(D$date_of_death, "%d/%m/%Y") >= "2020/03/16",]
D1 <- D1 %>% rename(anonymized_patient_id = ID)

data1 <- data %>% merge(D1[,c("anonymized_patient_id", "date_of_death","cause_icd10")], by="anonymized_patient_id", all.x=T)

data2 <- data1 %>% group_by(anonymized_patient_id) %>% mutate(death = case_when(!is.na(date_of_death) ~ 1,
                                                                                as.Date(covid19_test_date) >= max(as.Date(data1$date_of_death, format="%d/%m/%Y"), na.rm=T) ~ -1,
                                                                                TRUE ~ 0),
                                                              cause_of_death = case_when(any(cause_icd10 %in% c("U071", "U072")) ~ 1,
                                                                                         death == 1 ~ 2,
                                                                                         TRUE ~ 0)) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>% select(-c("cause_icd10"))

data <- data2
table(data[,c("death","cause_of_death")])

#critical care
c <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/hesin_critical_20210122.txt.gz", sep="\t")
c <- c[c$eid %in% data$anonymized_patient_id]
c <- c %>% rename(anonymized_patient_id = eid)
data1 <- data %>% merge(c, by="anonymized_patient_id", all.x=T)
data1 <- data1 %>% filter(as.Date(ccstartdate, "%d/%m/%Y") <= 30 + as.Date(covid19_test_date))
data1 <- data1 %>% filter(as.Date(ccdisdate, "%d/%m/%Y") >= as.Date(covid19_test_date))
data1 <- data1 %>% select(c("ccstartdate","ccdisdate","anonymized_patient_id", "bressupdays", "aressupdays", "acardsupdays", "rensupdays")) %>%
  mutate(icu_duration = as.numeric(as.Date(ccdisdate, "%d/%m/%Y") - as.Date(ccstartdate, "%d/%m/%Y")))

crit_id <- unique(data1$anonymized_patient_id)
data2 <- data %>% merge(data1, by="anonymized_patient_id", all.x=T) %>%
  group_by(anonymized_patient_id) %>%
  mutate_at(vars(c("bressupdays", "aressupdays", "acardsupdays", "rensupdays", "icu_duration")), max, na.rm=T) %>% 
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE)

data2 <- data2 %>% mutate(icu_admit = case_when(anonymized_patient_id %in% crit_id ~ 1,
                                                as.Date(covid19_test_date) >= max(as.Date(data1$ccstartdate, "%d/%m/%Y")) ~ -1,
                                                TRUE ~ 0),
                          hospitalization_start = case_when(is.na(hospitalization_start) & icu_admit == 1 ~ as.Date(ccstartdate, "%d/%m/%Y"),
                                                            TRUE ~ hospitalization_start),
                          hospitalization_end = case_when(is.na(hospitalization_end) & icu_admit == 1 ~ as.Date(ccdisdate, "%d/%m/%Y"),
                                                            TRUE ~ hospitalization_end),
                          icu_duration = ifelse(is.infinite(icu_duration), NA, icu_duration),
                          days_ventilator = ifelse(is.infinite(aressupdays), NA, aressupdays))

data2 <- data2 %>% mutate(highest_who_score = case_when(cause_of_death == 1 ~ 10,
                                                        acardsupdays > 0 & aressupdays > 0 ~ 9,
                                                        acardsupdays > 0 ~ 8,
                                                        aressupdays > 0 ~ 7,
                                                        bressupdays > 0 ~ 5,
                                                        hospitalization == 1 ~ 4,
                                                        hospitalization == 0 ~ 1,
                                                        TRUE ~ -1
                                                        ))

data2 <- data2 %>% mutate(highest_respiratory_support = case_when(aressupdays > 0 ~ 2,
                                                        bressupdays > 0 ~ 0,
                                                        hospitalization == 1 ~ -7,
                                                        hospitalization == 0 ~ -7,
                                                        TRUE ~ -1
))

data <- data2 %>% select(-c("ccstartdate","ccdisdate","bressupdays","aressupdays","acardsupdays"))
data <- data %>% mutate(date_of_death = as.Date(date_of_death, "%d/%m/%Y"))
data <- data %>% filter(covid19_test_date <= max(as.Date(data$date_of_death), na.rm = T))

##resp_severe
data$resp_severe <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  if(dim(TMP)>0){
    tmp1 <- tmp %>% filter(diag_icd10 %in% c("U071", "U072") | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
    tmp1 <- tmp1 %>%  filter(diag_icd10 %in% c("U071", "U072") | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
    data$resp_severe[data$anonymized_patient_id == unique(hesin_new$eid)[i] & any(tmp1$diag_icd10 %in% c("J80","J9600","J9609","Z991"))] <- 1
  }
}
for(i in seq(1,length(unique(hesin_new1$eid)))){
  tmprev <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],] %>% filter(diag_icd10 %in% c("U071", "U072")) %>% select(ins_index)
  tmp <- hesin_new1[hesin_new1$eid == unique(hesin_new1$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new1$eid)[i],]
  if(dim(TMP)>0){
  tmp1 <- tmp %>% filter(ins_index %in% tmprev$ins_index | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
  tmp1 <- tmp1 %>%  filter(ins_index %in% tmprev$ins_index | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
  data$resp_severe[data$anonymized_patient_id == unique(hesin_new1$eid)[i] & any(tmp1$oper4 %in% c("E851","E852"))] <- 1
  }
}
table(data$resp_severe)#209

data <- data %>% mutate(resp_severe = case_when(highest_respiratory_support > 0 ~ 1,
                                                resp_severe == 1 ~ 1,
                                                as.Date(covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T) ~ -1,
                                                death == 1 ~ -1,
                                                TRUE ~ 0))
table(data$resp_severe)#215

##resp_mild
data$resp_mild <- 0
data <- data %>% mutate(resp_mild = case_when(highest_respiratory_support == 0 ~ 1,
                                              resp_severe == 1 ~ 1,
                                              as.Date(covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T) ~ -1,
                                              death == 1 ~ -1,
                                              TRUE ~ 0))

data <- data %>% mutate(highest_respiratory_support = case_when(highest_respiratory_support == 2 ~ 2,
                                                                resp_severe == 1 ~ 1,
                                                                resp_mild == 1 ~ 0,
                                                                TRUE ~ highest_respiratory_support))

table(data[c("resp_severe","highest_respiratory_support")])
table(data[c("resp_mild","highest_respiratory_support")])

data$resp_severe[data$resp_mild == 1 & data$highest_respiratory_support == 0] <- -1

##hosp_dvt I81 I82* 
data$hosp_dvt <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  if(dim(TMP)>0){
  tmp1 <- tmp %>% filter(diag_icd10 %in% c("U071", "U072") | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
  tmp1 <- tmp1 %>%  filter(diag_icd10 %in% c("U071", "U072") | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
  data$hosp_dvt[data$anonymized_patient_id == unique(hesin_new$eid)[i] & "I81" %in% tmp1$diag_icd10] <- 1
  data$hosp_dvt[data$anonymized_patient_id == unique(hesin_new$eid)[i] & any(grepl("^I82", tmp1$diag_icd10))] <- 1
  }
}
sum(data$hosp_dvt)#7
data$hosp_dvt[data$hosp_dvt == 0 & as.Date(data$covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T)] <- -1
table(data$hosp_dvt)
##hosp_thrombo I26*
data$hosp_thrombo <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  if(dim(TMP)>0){
  tmp1 <- tmp %>% filter(diag_icd10 %in% c("U071", "U072") | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
  tmp1 <- tmp1 %>%  filter(diag_icd10 %in% c("U071", "U072") | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
  if(any(grepl("^I26", tmp1$diag_icd10))){
    data$hosp_thrombo[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
  }
}
sum(data$hosp_thrombo)#25
data$hosp_thrombo[data$hosp_thrombo == 0 & as.Date(data$covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T)] <- -1
table(data$hosp_thrombo)
##hosp_stroke I61,I62, I63, I64,I65,I66*, 
data$hosp_stroke <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  if(dim(TMP)>0){
  tmp1 <- tmp %>% filter(diag_icd10 %in% c("U071", "U072") | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
  tmp1 <- tmp1 %>%  filter(diag_icd10 %in% c("U071", "U072") | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
  if(any(grepl("^I60", tmp1$diag_icd10), grepl("^I61", tmp1$diag_icd10), grepl("^I62", tmp1$diag_icd10),
         grepl("^I63", tmp1$diag_icd10), grepl("^I64", tmp1$diag_icd10), grepl("^I65", tmp1$diag_icd10),
         grepl("^I66", tmp1$diag_icd10))){
    data$hosp_stroke[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
  }
}
sum(data$hosp_stroke)#38
data$hosp_stroke[data$hosp_stroke == 0 & as.Date(data$covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T)] <- -1
table(data$hosp_stroke)
#hosp_infarction I21*
data$hosp_infarction <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  if(dim(TMP)>0){
  tmp1 <- tmp %>% filter(diag_icd10 %in% c("U071", "U072") | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
  tmp1 <- tmp1 %>%  filter(diag_icd10 %in% c("U071", "U072") | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
  if(any(grepl("^I21", tmp1$diag_icd10))){
    data$hosp_infarction[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
  }
}
sum(data$hosp_infarction)#15
data$hosp_infarction[data$hosp_infarction == 0 & as.Date(data$covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T)] <- -1
#hosp_renal N17*
data$hosp_renal <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  if(dim(TMP)>0){
  tmp1 <- tmp %>% filter(diag_icd10 %in% c("U071", "U072") | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
  tmp1 <- tmp1 %>%  filter(diag_icd10 %in% c("U071", "U072") | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
  if(any(grepl("^N17", tmp1$diag_icd10))){
    data$hosp_renal[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
  }
}
sum(data$hosp_renal)#209

data <- data %>% mutate(hosp_renal = case_when(rensupdays > 0 ~ 1,
                                                TRUE ~ hosp_renal))
sum(data$hosp_renal)#216
data$hosp_renal[data$hosp_renal== 0 & as.Date(data$covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T)] <- -1

#hosp_bleeding I60*, I61*, I62*, K250 K252 K260 K262 K270 K272 K280 K282 K625 K922 I850 
data$hosp_bleeding <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  if(dim(TMP)>0){
  tmp1 <- tmp %>% filter(diag_icd10 %in% c("U071", "U072") | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
  tmp1 <- tmp1 %>%  filter(diag_icd10 %in% c("U071", "U072") | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
  if(any(grepl("^I60", tmp1$diag_icd10) | grepl("^I61", tmp1$diag_icd10) | grepl("^I62", tmp1$diag_icd10) |
         grepl("K250", tmp1$diag_icd10) | grepl("K252", tmp1$diag_icd10) | grepl("K260", tmp1$diag_icd10) |
         grepl("K262", tmp1$diag_icd10) | grepl("K272", tmp1$diag_icd10) | grepl("K270", tmp1$diag_icd10) |
         grepl("K280", tmp1$diag_icd10) | grepl("K282", tmp1$diag_icd10) | grepl("K625", tmp1$diag_icd10) | 
         grepl("K922", tmp1$diag_icd10) | grepl("I850", tmp1$diag_icd10))){
    data$hosp_bleeding[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
  }
}
sum(data$hosp_bleeding)#22
data$hosp_bleeding[data$hosp_bleeding== 0 & as.Date(data$covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T)] <- -1

#hosp_hepatic K720
data$hosp_hepatic <- 0
for(i in seq(1,length(unique(hesin_new$eid)))){
  tmp <- hesin_new[hesin_new$eid == unique(hesin_new$eid)[i],]
  TMP <- data[data$anonymized_patient_id == unique(hesin_new$eid)[i],]
  if(dim(TMP)>0){
  tmp1 <- tmp %>% filter(diag_icd10 %in% c("U071", "U072") | as.Date(epistart, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date) + 30 | as.Date(admidate, "%d/%m/%Y") <= as.Date(TMP$covid19_test_date)  + 30 )
  tmp1 <- tmp1 %>%  filter(diag_icd10 %in% c("U071", "U072") | as.Date(epiend, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date) | as.Date(disdate, "%d/%m/%Y") >= as.Date(TMP$covid19_test_date))
  if(any(grepl("K720", tmp1$diag_icd10))){
    data$hosp_hepatic[data$anonymized_patient_id == unique(hesin_new$eid)[i]] <- 1
  }
  }
}
sum(data$hosp_hepatic)#1
data$hosp_hepatic[data$hosp_hepatic == 0 & as.Date(data$covid19_test_date) > max(as.Date(data$date_of_death), na.rm = T)] <- -1

data <- data %>% select(-c(SEX,CENTRE,rensupdays))

hosps <- c("resp_severe","resp_mild","hosp_dvt","hosp_thrombo","hosp_stroke","hosp_infarction","hosp_renal","hosp_bleeding","hosp_hepatic")

data1 <- data %>% mutate_at(.vars = c(hosps), .funs = funs(ifelse(as.Date(covid19_test_date) > as.Date("2020-12-09"), -1, .)))
data1 <- data1 %>% filter(death != -1)

saveRDS(data1, file="scratch/clinical_value_20210217.rds")
