setwd("/home/tomoko/02.associations")
data <- readRDS("../01.data.QC/curated_clinical/complication_dat.rds")
geno <- readRDS("../covid19-hgi-clinical-values/bgen/all_ans_pca_raw.rds")
data1 <- data %>% inner_join(geno, by="anonymized_patient_id")
data1 %>% saveRDS("../01.data.QC/curated_clinical/complication_dat_geno.rds")
data1 <- data1 %>% mutate_at(.vars = c(colnames(data1)[grepl("com_",colnames(data1))]), .funs = funs(as.numeric(.)))
ukb <- readRDS("../01.data.QC/ukbb/clinical_value_20210217.rds")
genoukb <- readRDS("../01.data.QC/ukbb/all_chr.rds")
pcukb <- fread("~/11.pca/cases_controls_pca_sup_pops_0.5_probs.txt")

genoukb <- genoukb %>% merge(pcukb, by.x="FID",by.y="s")

ukb <- ukb %>% inner_join(genoukb, by=c("anonymized_patient_id"="FID")) %>% mutate(anonymized_patient_id = as.character(anonymized_patient_id))
ukb <- ukb  %>% mutate(anonymized_patient_id = paste0("UKB_",anonymized_patient_id))
ukb <- ukb %>% mutate(age2 = age_at_diagnosis^2)

ukb <- ukb %>% rename(hepatic = hosp_hepatic, aki = hosp_renal) %>% mutate(a1 = case_when(death == 1 ~ 1,
                                                                                          resp_severe == 1 ~ 1,
                                                                                          highest_who_score >= 6 ~ 1,
                                                                                          resp_severe == 0 & death == 0 ~ 0,
                                                                                          TRUE ~ -1),
                                                                           vte = case_when(hosp_thrombo == 1 ~ 1,
                                                                                           hosp_dvt == 1 ~ 1,
                                                                                           hosp_thrombo == 0 ~ 0,
                                                                                           hosp_dvt == 0 ~ 0,
                                                                                           TRUE ~ -1),
                                                                           cvd = case_when(hosp_infarction == 1 ~ 1,
                                                                                           hosp_stroke == 1 ~ 1,
                                                                                           hosp_infarction == 0 ~ 0,
                                                                                           hosp_stroke == 0 ~ 0,
                                                                                           TRUE ~ -1))

ukb <- ukb %>% select(any_of(colnames(data1))) %>% mutate(anonymized_patient_id = as.character(anonymized_patient_id))
ukb <- ukb %>% mutate(study = "UKB")
final <- bind_rows(ukb, data1)
final <- final %>% mutate(a1 = case_when(death == 1 ~ 1,
                                         resp_severe == 1 ~ 1,
                                         resp_severe == 0 & death == 0 ~ 0,
                                         TRUE ~ -1))
final$hospitalization[final$a1 == 1] <- 1
final$hospitalization[final$death == 1] <- 1
final$hospitalization[final$icu_admit == 1] <- 1
final <- final %>% mutate(cardiac = case_when(cvd == 1 ~ 1,
                                              hosp_cardiac == 1 ~ 1,
                                              hosp_cardiac == 0 & death == 0 ~ 0,
                                              cvd == 0 & death == 0 ~ 0,
                                              TRUE ~ -1),
                          icu_admit_rev = case_when(icu_admit == 1 ~ 1,
                                                    is.na(hospitalization) ~ -1,
                                                    icu_admit == 0 & hospitalization == 1 ~ -1,
                                                icu_admit == -1 & hospitalization == 0 ~ 0,
                                                is.na(icu_admit) & hospitalization == 0 ~ 0,
                                                TRUE ~ icu_admit))

final <- final %>% mutate(vte = case_when(vte == 1 ~ 1,
                                          vte == 0 & death == 0 ~ 0,
                                          hospitalization == 0 & death == 0 ~ 0,
                                          TRUE ~ -1),
                          cvd = case_when(cvd == 1 ~ 1,
                                          cvd == 0 & death == 0 ~ 0,
                                          hospitalization == 0 & death == 0 ~ 0,
                                          TRUE ~ -1),
                          aki = case_when(aki == 1 ~ 1,
                                          aki == 0 & death == 0 ~ 0,
                                          hospitalization == 0 & death == 0 ~ 0,
                                          TRUE ~ -1),
                          hepatic = case_when(hepatic == 1 ~ 1,
                                              hepatic == 0 & death == 0 ~ 0,
                                              hospitalization == 0 & death == 0 ~ 0,
                                              TRUE ~ -1),
                          pe = case_when(hosp_thrombo == 1 ~ 1,
                                         vte == 0 & death == 0 ~ 0,
                                         hospitalization == 0 & death == 0 ~ 0,
                                         TRUE ~ -1))
final <- final %>% mutate(highest_who_score = case_when(death == 1 ~ 10,
                                                        resp_severe == 1 & highest_who_score >= 6 ~ highest_who_score,
                                                        resp_severe == 1 & highest_who_score < 6 ~ 6,
                                                        resp_mild == 1 & highest_who_score > 5 ~ highest_who_score,
                                                        resp_mild == 1 & highest_who_score <= 5 ~ 5,
                                                        hospitalization == 1 & highest_who_score > 4 ~ highest_who_score,
                                                        hospitalization == 1 & highest_who_score <= 4 ~ 4,
                                                        TRUE ~ highest_who_score
))
final <- final %>% mutate(death = case_when(highest_who_score == 10 ~ 1,
                                            TRUE ~ death))

table(final[,c("highest_who_score","highest_respiratory_support")])
table(final[,c("highest_who_score","resp_severe")])
table(final[,c("highest_who_score","death")])

final <- final %>% mutate(hospitalization = case_when(resp_mild == 1 ~ 1,
                                                      TRUE ~ hospitalization))
final <- final %>% mutate(nation = case_when(study == "UKB" ~ "UK",
                                             grepl("belgium", study) ~ "Belgium",
                                             grepl("brasil", study) ~ "Union",
                                             grepl("canada", study) ~ "Union",
                                             grepl("germany", study) ~ "Union",
                                             grepl("italy", study) ~ "Italy",
                                             grepl("norway", study) ~ "Union",
                                             grepl("spain", study) ~ "Spain",
                                             grepl("sweden", study) ~ "Union"
                                             ))

data_EUR <- final %>% filter(pop == "EUR")
saveRDS(final, "../01.data.QC/curated_clinical/complication_dat_final.rds")

