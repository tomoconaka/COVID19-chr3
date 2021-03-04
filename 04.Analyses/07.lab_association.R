setwd("/home/tomoko/02.associations")
final <- readRDS("../01.data.QC/curated_clinical/complication_dat_final.rds")
snps <- colnames(final)[grepl("chr", colnames(final))][3:12]
outcome <- c("resp_severe","resp_mild", "icu_admit" ,"vte", "cardiac",
             "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")
final <- final %>% mutate_at(.vars = c(outcome), .funs = funs(ifelse(. == -1, NA, .)))
final <- final %>% mutate(age2 = age_at_diagnosis*age_at_diagnosis)

identification <- c("anonymized_patient_id")

lab_data <- readRDS("../01.data.QC/curated_clinical/all_lab_data.rds")

hospital_info <- 
  c("a1")

covid19 <-
  c("covid19_test_date")

blood_values <- c("lab_result_date",
                  "lab_wbc",
                  "lab_lymphocytes",
                  "lab_neutrophils",
                  "lab_monocytes",
                  "lab_platelets",
                  "lab_crp",
                  "lab_trop_t",
                  "lab_ast",
                  "lab_alt",
                  "lab_bilirubin",
                  "lab_ldh",
                  "lab_ggt",
                  "lab_alp",
                  "lab_d_dimer",
                  "lab_il_6",
                  "lab_serum_ferritin",
                  "lab_procalcitonin",
                  "lab_ck",
                  "lab_fibrinogen",
                  "lab_creatinine")

data <- final %>% select(-c(colnames(final)[grepl("lab_",colnames(final))]))
data <- lab_data %>% merge(data, by="anonymized_patient_id", all.x=T) 

data_EUR <- data %>% filter(pop == "EUR")

data_EUR <- data_EUR %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))

blood_values_max <- c("lab_wbc","lab_neutrophils","lab_monocytes","lab_platelets","lab_crp","lab_trop_t",      
                      "lab_ast","lab_alt","lab_bilirubin","lab_ldh","lab_ggt","lab_alp",
                      "lab_d_dimer","lab_il_6","lab_serum_ferritin", "lab_procalcitonin", 
                      "lab_ck" ,"lab_creatinine","lab_fibrinogen")

blood_values_min <- c("lab_lymphocytes")

for(j in c(1)){
  out <- data.frame(matrix(0,19,6))
  colnames(out) <- c("biomarker", "beta", "se", "pvalue","N","SNP")
  for(i in c(1:19)){
    df <- data_EUR %>% mutate(time = case_when(!is.na(covid19_test_date) & !is.na(lab_result_date) ~ as.numeric(as.Date(lab_result_date) - as.Date(covid19_test_date)),
                                               TRUE ~ 0)) %>% filter(time >= -2 & time <= 30) %>% group_by(anonymized_patient_id) %>%
      mutate_at(vars(c(any_of(blood_values_max[i]))),max,na.rm=T) %>% 
      ungroup()  %>% 
      distinct(anonymized_patient_id, .keep_all=TRUE) %>%
      pivot_longer(c(any_of(blood_values_max[i])))
    df <- df %>% filter(!is.na(value) & !is.infinite(value))
    df$value <- scale(log(df$value))
    LM <- glmer(as.formula(paste0("value ~ snp + age_at_diagnosis + sex + (1 | study) + PC1 + PC2 + PC3 + PC4 + PC5")), dat=df,family = "gaussian")
    out[i,1] <- blood_values_max[i]
    out[i,2] <- summary(LM)$coefficient["snp",1]
    out[i,3] <- summary(LM)$coefficient["snp",2]
    out[i,4] <- parameters::p_value(LM)[2,2]
    out[i,5] <- dim(df[!is.na(df$value),])[1]
    out[i,6] <- snps[j]
    }
  print(out[order(out$pvalue),])
  write.table(out, file="results/worst_EUR_biomarker.results", append=T,quote = F, col.names = F, row.names = F, sep="\t")
}


for(j in c(1)){
  out <- data.frame(matrix(0,1,6))
  colnames(out) <- c("biomarker", "beta", "se", "pvalue","N","SNP")
  for(i in c(1:1)){
    df <- data_EUR %>% mutate(time = case_when(!is.na(covid19_test_date) & !is.na(lab_result_date) ~ as.numeric(as.Date(lab_result_date) - as.Date(covid19_test_date)),
                                               TRUE ~ 0)) %>% filter(time >= -2 & time <= 30) %>% group_by(anonymized_patient_id) %>%
      mutate_at(vars(c(any_of(blood_values_min[i]))),min,na.rm=T) %>% 
      ungroup()  %>% 
      distinct(anonymized_patient_id, .keep_all=TRUE) %>%
      pivot_longer(c(any_of(blood_values_min[i])))
    df <- df %>% filter(!is.na(value) & !is.infinite(value))
    df$value <- scale(log(df$value))
    LM <- glmer(as.formula(paste0("value ~ snp + age_at_diagnosis + sex + (1 | study) + PC1 + PC2 + PC3 + PC4 + PC5")), dat=df,family = "gaussian")
    out[i,1] <- blood_values_min[i]
    out[i,2] <- summary(LM)$coefficient["snp",1]
    out[i,3] <- summary(LM)$coefficient["snp",2]
    out[i,4] <- parameters::p_value(LM)[2,2]
    out[i,5] <- dim(df[!is.na(df$value),])[1]
    out[i,6] <- snps[j]
  }
  print(out[order(out$pvalue),])
  write.table(out, file="results/worst_EUR_biomarker.results", append=T,quote = F, col.names = F, row.names = F, sep="\t")
}

out <- fread("results/worst_EUR_biomarker.results")
colnames(out) <- c("biomarker", "beta", "se", "pvalue","N","SNP")
write.table(out, file="results/worst_EUR_biomarker.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")

out <- data.frame(matrix(0,19,6))
colnames(out) <- c("biomarker", "beta", "se", "pvalue","N","SNP")
for(i in c(1:19)){
    df <- data_EUR %>% mutate(time = case_when(!is.na(covid19_test_date) & !is.na(lab_result_date) ~ as.numeric(as.Date(lab_result_date) - as.Date(covid19_test_date)),
                                               TRUE ~ 0)) %>% filter(time >= -2 & time <= 30) %>% group_by(anonymized_patient_id) %>%
      mutate_at(vars(c(any_of(blood_values_max[i]))),max,na.rm=T) %>% 
      ungroup()  %>% 
      distinct(anonymized_patient_id, .keep_all=TRUE) %>%
      pivot_longer(c(any_of(blood_values_max[i])))
    df <- df %>% filter(!is.na(value) & !is.infinite(value))
    df$value <- log(df$value)
    df <- df %>% drop_na(value)
    LM <- glm(as.formula(paste0("value ~ study")), dat=df,family = "gaussian")
    df$value <- scale(residuals(LM))
    LM <- glmer(as.formula(paste0("a1 ~ value + age_at_diagnosis + sex + (1 | nation)")), dat=df,family = "binomial")
    out[i,1] <- blood_values_max[i]
    out[i,2] <- summary(LM)$coefficient[2,1]
    out[i,3] <- summary(LM)$coefficient[2,2]
    out[i,4] <- summary(LM)$coefficient[2,4]
    out[i,5] <- dim(df[!is.na(df$value),])[1]
    out[i,6] <- "a1"
  }
print(out[order(out$pvalue),])
write.table(out, file="results/worst_EUR_biomarker_clinical.results", append=T,quote = F, col.names = F, row.names = F, sep="\t")

out <- fread("results/worst_EUR_biomarker_clinical.results")
colnames(out) <- c("biomarker", "beta", "se", "pvalue","N","SNP")
write.table(out, file="results/worst_EUR_biomarker_clinical.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")
