setwd("/home/tomoko/02.associations")
final <- readRDS("../01.data.QC/curated_clinical/complication_dat_final.rds")

final <- final %>% mutate(icu_admit_rev = case_when(is.na(hospitalization) ~ -1,
                                                    icu_admit == 0 & hospitalization == 1 ~ 0,
                                                    hospitalization == 0 ~ 0,
                                                    TRUE ~ icu_admit))

snps <- colnames(final)[grepl("chr", colnames(final))][3:12]
outcome <- c("resp_severe","resp_mild", "icu_admit" ,"vte", "cardiac",
             "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")

library(tidyr)
library(dplyr)
library(table1)
final <- final %>% mutate(height = ifelse(height <0 ,NA, height),
                          weight = ifelse(weight <0, NA, weight))
final <- final %>% mutate(BMI = weight/((height/100)^2))


my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

rndr <- function(x, name, ...) {
  if (!is.numeric(x)) return(c(render.categorical.default(x)))
  what <- switch(name,
                 age_at_diagnosis = c("Mean (SD)", render.missing.default(x)),
                 BMI  = c("Mean (SD)",render.missing.default(x)))
  parse.abbrev.render.code(c("", what))(x)
}

comorb <- c("sex","smoking",colnames(final)[grepl("com_", colnames(final))])

final <- final %>% mutate_at(.vars=c(comorb), .funs=funs(ifelse(is.na(.), -1, .)))

final$sex <- 
  factor(final$sex, levels=c(1, -1),
         labels=c("Female", "Missing"))

final$smoking <- 
  factor(final$smoking, levels=c(0,1,2,-1),
         labels=c("Ever", 
                  "Ever",
                  "Never",
                  "Missing"))

label(final$sex) <- "Sex"
label(final$age_at_diagnosis) <- "Age (years)"
label(final$smoking) <- "Smoking status"

final$com_diabetes <- 
  factor(final$com_diabetes, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_diabetes) <- "Diabetes Mellitus"

final$com_asthma <- 
  factor(final$com_asthma, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_asthma) <- "Asthma"

final$com_chronic_pulm <- 
  factor(final$com_chronic_pulm, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_chronic_pulm) <- "COPD"

final$com_liver <- 
  factor(final$com_liver, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_liver) <- "Liver disease"

final$com_chronic_kidney <- 
  factor(final$com_chronic_kidney, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_chronic_kidney) <- "Chronic kidney disease"

final$com_heart_failure <- 
  factor(final$com_heart_failure, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_heart_failure) <- "Chronic heart failure"

final$com_stroke <- 
  factor(final$com_stroke, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_stroke) <- "Stroke"

final$com_hypertension <- 
  factor(final$com_hypertension, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_hypertension) <- "Hypertenstion"

final$com_cancer <- 
  factor(final$com_cancer, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_cancer) <- "Cancer"

final$com_dementia <- 
  factor(final$com_dementia, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_dementia) <- "Dementia"

final$com_immunocomp <- 
  factor(final$com_immunocomp, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_immunocomp) <- "Immunocompromised state"


final$com_transplant <- 
  factor(final$com_transplant, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$com_transplant) <- "Transplantation"

final$pop <- 
  factor(final$pop, levels=unique(final$pop),
         labels=c("European", "South Asian", "African", "others", "East Asian","Admixed American"))
label(final$pop) <- "Ancestry"

#final$hospitalization[is.na(final$hospitalization)] <- -1
final$hospitalization <- 
  factor(final$hospitalization, levels=c(1),
       labels=c("Hospitalized"))
label(final$hospitalization) <- "Hospitalization"

final$icu_admit_rev <- 
  factor(final$icu_admit_rev, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$icu_admit_rev) <- "ICU admission"

final$resp_severe[is.na(final$resp_severe)] <- -1
final <- final %>% mutate(resp = case_when(resp_severe == 1 ~ 2,
                                           resp_mild == 1 ~ 1,
                                           resp_severe == 0 ~ 0,
                                           resp_mild == 0 ~ 0,
                                           TRUE ~ -1))
  
final$resp <- 
  factor(final$resp, levels=c(2, 1, -1),
         labels=c("Severe respiratory failure", "Oxygen supplementation", "Missing"))
label(final$resp) <- "Respiratory failure"

final$death[is.na(final$death)] <- -1

final$death <- 
  factor(final$death, levels=c(0, 1, -1),
         labels=c("Survived", "Deceased","Missing"))
label(final$death) <- "Death Status"

final$hepatic[is.na(final$hepatic)] <- -1
final$hepatic <- 
  factor(final$hepatic, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$hepatic) <- "Hepatic injury"

final$cardiac[is.na(final$cardiac)] <- -1
final$cardiac <- 
  factor(final$cardiac, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$cardiac) <- "Cardiovascular complications"

final$aki[is.na(final$aki)] <- -1
final$aki <- 
  factor(final$aki, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$aki) <- "Kidney injury"

final$vte[is.na(final$vte)] <- -1
final$vte <- 
  factor(final$vte, levels=c(1, -1),
         labels=c("Present", "Missing"))
label(final$vte) <- "Venous thromboembolism"

label(final$highest_who_score) <- "WHO severity score"

table1(~ sex + age_at_diagnosis + pop + icu_admit_rev + death + resp + 
         hepatic + cardiac + aki + vte | hospitalization, data=final, overall = "Total", render = rndr, render.missing.default = TRUE)

final <- final %>% mutate(STUDY = case_when(study == "UKB" ~ "UKB",
                                            study == "belgium_migeotte" ~ "BelCovid_1",
                                            study == "belgium_rahmouni" ~ "BelCovid_2",
                                            study == "brasil" ~ "BRACOVID",
                                            study == "canada" ~ "BQC19",
                                            study == "germany_ludwig" ~ "BosCo",
                                            study == "germany_schulte" ~ "COMRI",
                                            study == "italy_duga" ~ "COVID19-Host(a)ge_4",
                                            study == "italy_renieri" ~ "GEN-COVID",
                                            study == "italy_valenti" ~ "FoGS",
                                            study == "norway" ~ "NorCoV2",
                                            study == "spain_alarcon" ~ "SPGRX",
                                            study == "spain_bujanda" ~ "COVID19-Host(a)ge_1",
                                            study == "spain_butti" ~ "COVID19-Host(a)ge_2",
                                            study == "spain_gomez" ~ "COVID19-Host(a)ge_3",
                                            study == "spain_planas" ~ "INMUNGEN-CoV2",
                                            study == "sweden" ~ "SweCovid"
                                            ))

table1(~ sex + age_at_diagnosis + pop + BMI + smoking + 
         com_cancer + com_chronic_kidney + com_chronic_pulm + com_heart_failure + com_transplant +
         com_diabetes + hospitalization + icu_admit_rev + death + resp + 
         hepatic + cardiac + aki + vte | STUDY, data=final, overall = "Total", render = rndr, render.missing.default = TRUE)
