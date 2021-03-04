setwd("/home/tomoko/02.associations")
demog_data <-  readRDS("../01.data.QC/curated_clinical/demog_data.rds")
hosp_data <- readRDS("../01.data.QC/curated_clinical/hosp_data.rds")
comorb_data <- readRDS("../01.data.QC/curated_clinical/comorb_data.rds")
covid_data <- readRDS("../01.data.QC/curated_clinical/covid_data.rds")

tmp <- covid_data %>% mutate(time =  covid19_test_date - covid19_first_symptoms_date)
mean(tmp$time, na.rm=T)#6.602385
lab_data <- readRDS("../01.data.QC/curated_clinical/all_lab_data.rds")

identification <- c("anonymized_patient_id")

hospital_info <- 
  c("icu_admit","icu_duration", "death","date_of_death","cause_of_death","hospitalization_duration",
    "highest_respiratory_support","days_ventilator","hosp_hepatic",
    "highest_who_score", "hosp_dvt", "hosp_thrombo", "hosp_stroke",
    "hosp_infarction", "hosp_bleeding", "hospitalization","hospitalization_start","hosp_renal","hospitalization_end")

covid_info <- c("covid19_test_date","covid19_test","covid19_first_symptoms_date")

blood_values <- c("lab_result_date","lab_lymphocytes","lab_neutrophils","lab_monocytes","lab_eosinophils",
                  "lab_basophils","lab_alt","lab_trop_t","lab_trop_i","lab_creatinine")

data <- demog_data %>% inner_join(hosp_data[c(identification,hospital_info)], by="anonymized_patient_id")
data <- data %>% merge(covid_data[c(identification,covid_info)], by="anonymized_patient_id")

data <- data %>% merge(comorb_data[,-15], by="anonymized_patient_id")

data <- data %>% filter(covid19_test == 1)
df <- lab_data %>% select(c(identification,blood_values)) %>% merge(data, all = T) %>% 
  mutate(time = case_when(!is.na(covid19_test_date) & !is.na(lab_result_date) ~ as.numeric(as.Date(lab_result_date) - as.Date(covid19_test_date)),
                          TRUE ~ 0)) %>% filter(time >= -2 & time <= 30) %>% group_by(anonymized_patient_id) %>%
  mutate_at(vars(c(any_of(blood_values[c(-1,-2)]))),max,na.rm=T) %>% 
  mutate_at(vars(c(any_of(blood_values[2]))),min,na.rm=T) %>%
  ungroup()  %>% 
  distinct(anonymized_patient_id, .keep_all=TRUE) %>% select(c(identification,blood_values))
data <- data %>% merge(df, by="anonymized_patient_id", all.x=T)
#hosp_cardiac
data <- data %>% mutate(hosp_cardiac = case_when(lab_trop_i > 16 ~ 1,
                                                 study == "spain_bujanda" & lab_trop_t >= 15 ~ 1,
                                                 study == "spain_butti" & lab_trop_t >= 45 ~ 1,
                                                 study == "italy_renieri" & hosp_infarction == 1 ~ 1,
                                                 study == "italy_valenti" & lab_trop_t > 0.40 ~ 1,
                                                 study == "belgium_migeotte" & lab_trop_t > 15 ~ 1,
                                                 study == "brasil" & lab_trop_t > 15 ~ 1,
                                                 study == "canada" & lab_trop_t > 15 ~ 1,
                                                 study == "germany_ludwig" & lab_trop_t > 14 ~ 1,
                                                 study == "germany_schulte" & lab_trop_t > 14 ~ 1,
                                                 study == "belgium_rahmouni" & lab_trop_i > 34.20 ~ 1,
                                                 study == "norway" & lab_trop_t > 30 ~ 1,
                                                 study == "spain_gomez" & lab_trop_t > 14 & grepl("HUVR_", anonymized_patient_id) ~ 1,
                                                 study == "spain_gomez" & lab_trop_t > 19.8 & !grepl("HUVR_", anonymized_patient_id) ~ 1,
                                                 study == "italy_duga" & lab_trop_t  > 19.8 ~ 1,
                                                 is.na(lab_trop_i) & is.na(lab_trop_t) ~ -1,
                                                 TRUE ~ 0))

#cre > 3 x ULN def
data <- data %>% mutate(hosp_aki3 = case_when(com_chronic_kidney == 1 ~ -1,
                                              study == "spain_bujanda" & lab_creatinine > 1.2*3 ~ 1,
                                              study == "spain_butti" & lab_creatinine >= 0.95*3 ~ 1,
                                              study == "italy_renieri" & lab_creatinine > 1.2*3 ~ 1,
                                              study == "italy_valenti" & lab_creatinine > 1.2*3 ~ 1,
                                              study == "belgium_migeotte" & lab_creatinine > 1.2*3 ~ 1,
                                              study == "brasil" & lab_creatinine > 1.2*3 ~ 1,
                                              study == "canada" & lab_creatinine > 1.2*3 ~ 1,
                                              study == "sweden" & lab_creatinine > 105*0.0113*3 ~ 1,
                                              study == "germany_ludwig" & lab_creatinine > 1.2*3 ~ 1,
                                              study == "germany_schulte" & lab_creatinine > 1.2*3 ~ 1,
                                              study == "belgium_rahmouni" & lab_trop_i > 1.18*3 ~ 1,
                                              study == "norway" & lab_trop_t > 1.2*3 ~ 1,
                                              study == "spain_gomez" & lab_creatinine > 1.5*3 & grepl("HUVR_", anonymized_patient_id) ~ 1,
                                              study == "spain_gomez" & lab_creatinine > 1.2*3 & !grepl("HUVR_", anonymized_patient_id) ~ 1,
                                              study == "italy_duga" & lab_creatinine > 1.17*3 ~ 1,
                                              is.na(lab_creatinine) ~ -1,
                                              TRUE ~ 0))

#cre > 1.5 x ULN def
data <- data %>% mutate(hosp_aki2 = case_when(com_chronic_kidney == 1 ~ -1,
                                              study == "spain_bujanda" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "spain_butti" & lab_creatinine >= 0.95*1.5 ~ 1,
                                              study == "italy_renieri" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "italy_valenti" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "belgium_migeotte" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "brasil" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "canada" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "sweden" & lab_creatinine > 105*0.0113*1.5 ~ 1,
                                              study == "germany_ludwig" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "germany_schulte" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "belgium_rahmouni" & lab_creatinine > 1.18*1.5 ~ 1,
                                              study == "norway" & lab_creatinine > 1.2*1.5 ~ 1,
                                              study == "spain_gomez" & lab_creatinine > 1.5*1.5 & grepl("HUVR_", anonymized_patient_id) ~ 1,
                                              study == "spain_gomez" & lab_creatinine > 1.2*1.5 & !grepl("HUVR_", anonymized_patient_id) ~ 1,
                                              study == "italy_duga" & lab_creatinine > 1.17*1.5 ~ 1,
                                              is.na(lab_creatinine) ~ -1,
                                              TRUE ~ 0))
table(data$hosp_aki3)
table(data$hosp_aki2)

#liver ALT > 5*ULN
data <- data %>% mutate(hosp_liver3 = case_when(com_liver == 1 ~ -1,
                                                study == "spain_bujanda" & lab_alt > 41*5 ~ 1,
                                                study == "spain_butti" & lab_alt >= 35*5 ~ 1,
                                                study == "italy_renieri" & lab_alt > 40*5 ~ 1,
                                                study == "italy_valenti" & lab_alt > 41*5 ~ 1,
                                                study == "belgium_migeotte" & lab_alt > 41*5 ~ 1,
                                                study == "brasil" & lab_alt > 41*5 ~ 1,
                                                study == "canada" & lab_alt > 40*5 ~ 1,
                                                study == "sweden" & lab_alt > 1.2/0.0166*5 ~ 1,
                                                study == "germany_ludwig" & lab_alt > 56*5 ~ 1,
                                                study == "germany_schulte" & lab_alt > 50*5 ~ 1,
                                                study == "belgium_rahmouni" & lab_alt > 55*5 ~ 1,
                                                study == "norway" & lab_alt > 70*5 ~ 1,
                                                study == "spain_gomez" & lab_alt > 45*5 & grepl("HUVR_", anonymized_patient_id) ~ 1,
                                                study == "spain_gomez" & lab_alt > 40*5 & !grepl("HUVR_", anonymized_patient_id) ~ 1,
                                                study == "italy_duga" & lab_alt > 51*5 ~ 1,
                                                is.na(lab_alt) ~ -1,
                                                TRUE ~ 0))

table(data$hosp_liver3)

#liver ALT > 3*ULN
data <- data %>% mutate(hosp_liver2 = case_when(com_liver == 1 ~ -1,
                                                study == "spain_bujanda" & lab_alt > 41*3 ~ 1,
                                                study == "spain_butti" & lab_alt >= 35*3 ~ 1,
                                                study == "italy_renieri" & lab_alt > 40*3 ~ 1,
                                                study == "italy_valenti" & lab_alt > 41*3 ~ 1,
                                                study == "belgium_migeotte" & lab_alt > 41*3 ~ 1,
                                                study == "brasil" & lab_alt > 41*3 ~ 1,
                                                study == "canada" & lab_alt > 40*3 ~ 1,
                                                study == "sweden" & lab_alt > 1.2/0.0166*3 ~ 1,
                                                study == "germany_ludwig" & lab_alt > 56*3 ~ 1,
                                                study == "germany_schulte" & lab_alt > 50*3 ~ 1,
                                                study == "belgium_rahmouni" & lab_alt > 55*3 ~ 1,
                                                study == "norway" & lab_alt > 70*3 ~ 1,
                                                study == "spain_gomez" & lab_alt > 45*3 & grepl("HUVR_", anonymized_patient_id) ~ 1,
                                                study == "spain_gomez" & lab_alt > 40*3 & !grepl("HUVR_", anonymized_patient_id) ~ 1,
                                                study == "italy_duga" & lab_alt > 51*3 ~ 1,
                                                is.na(lab_alt) ~ -1,
                                                TRUE ~ 0))
table(data$hosp_liver2)
data <- data %>% mutate(resp_severe = case_when(highest_respiratory_support > 0 ~ 1,
                                                highest_respiratory_support == 0 ~ -1,
                                                highest_respiratory_support == -7 & death == 0 ~ 0,
                                                hospitalization == 0 & death == 0 ~ 0,
                                                TRUE ~ -1),
                        resp_mild = case_when(highest_respiratory_support >= 0 ~ 1,
                                              highest_respiratory_support == -7 & death == 0 ~ 0,
                                              hospitalization == 0 & death == 0 ~ 0,
                                              TRUE ~ -1))

##lymphocytopenia 0.5
data <- data %>% mutate(lymphocytopenia3 = case_when(lab_lymphocytes < 0.5 ~ 1,
                                                     lab_lymphocytes >= 0.5 ~ 0,
                                                     TRUE ~ -1))
data <- data %>% mutate(lymphocytopenia2 = case_when(lab_lymphocytes < 0.8 ~ 1,
                                                     lab_lymphocytes >= 0.8 ~ 0,
                                                     TRUE ~ -1))

data <- data %>% mutate(vte = case_when(hosp_thrombo == 1 ~ 1,
                                        hosp_dvt == 1 ~ 1,
                                        hosp_thrombo == 0 & death == 0 ~ 0,
                                        hosp_dvt == 0 & death == 0 ~ 0,
                                        hospitalization == 0 & death == 0 ~ 0,
                                        TRUE ~ -1),
                        cvd = case_when(hosp_infarction == 1 ~ 1,
                                        hosp_stroke == 1 ~ 1,
                                        hosp_cardiac == 1 ~ 1,
                                        hosp_infarction == 0 & death == 0 ~ 0,
                                        hosp_stroke == 0 & death == 0 ~ 0,
                                        hosp_cardiac == 0 & death == 0 ~ 0,
                                        hospitalization == 0 & death == 0 ~ 0,
                                        TRUE ~ -1),
                        aki = case_when(hosp_renal == 1 ~ 1,
                                        hosp_aki3 == 1 ~ 1,
                                        hosp_aki2 == 1 ~ 1,
                                        hosp_renal == 0 & death == 0 ~ 0,
                                        hosp_aki2 == 0 & death == 0 ~ 0,
                                        hosp_aki3 == 0 & death == 0 ~ 0,
                                        hospitalization == 0 & death == 0 ~ 0,
                                        TRUE ~ -1),
                        hepatic = case_when(hosp_hepatic == 1 & study != "italy_renieri" ~ 1,
                                            hosp_liver3 == 1 ~ 1,
                                            hosp_liver2 == 1 ~ 1,
                                            hosp_hepatic == 0 & death == 0 ~ 0,
                                            hosp_liver3 == 0 & death == 0 ~ 0,
                                            hosp_liver2 == 0 & death == 0 ~ 0,
                                            hospitalization == 0 & death == 0 ~ 0,
                                            TRUE ~ -1))

data %>% saveRDS("../01.data.QC/curated_clinical/complication_dat.rds")
