setwd("/home/tomoko/02.associations")
snps <- colnames(final)[grepl("chr", colnames(final))][3:12]
outcome <- c("resp_severe","resp_mild", "icu_admit_rev" ,"vte", "cardiac",
             "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")

final <- readRDS("../01.data.QC/curated_clinical/complication_dat_final.rds")
final <- final %>% mutate(height = ifelse(height == -1, NA, height))
final <- final %>% mutate(BMI = weight/((height/100)^2))
final <- final %>% mutate_at(.vars = c(outcome), .funs = funs(ifelse(. == -1, NA, .)))
final <- final %>% mutate(age2 = age_at_diagnosis*age_at_diagnosis)

comorb <- c("com_cancer","com_chronic_kidney","com_chronic_pulm",
            "com_transplant","com_heart_failure","com_diabetes")

final <- final %>% mutate(BMIhigh = case_when(BMI >= 30 ~ 1,
                                              BMI < 30 ~ 0,
                                              TRUE ~ -1),
                          smoke = case_when(smoking %in% c(0,1) ~ 1,
                                            smoking == 2 ~ 0,
                                            TRUE ~ -1))
final$BMIhigh[final$BMIhigh == -1] <- NA
final$smoke[final$smoke == -1] <- NA
final <- final %>% mutate(com = case_when(com_cancer == 1 ~ 1,
                                          com_chronic_kidney == 1 ~ 1,
                                          com_chronic_pulm == 1 ~ 1,
                                          com_transplant == 1 ~ 1,
                                          com_heart_failure == 1 ~ 1,
                                          com_diabetes == 1 ~ 1,
                                          BMI >= 30 ~ 1,
                                          smoking %in% c(0,1) ~ 1,
                                          com_cancer == 0 & com_chronic_kidney == 0 & com_chronic_pulm == 0 & com_transplant == 0 &
                                          com_heart_failure == 0 & com_diabetes == 0 & smoking == 2 & BMI < 30 ~ 0,
                                          TRUE ~ -1))

final$com[final$com == -1] <- NA

final <- final %>% mutate(com1 = com_cancer + com_chronic_kidney + com_chronic_pulm + 
                            com_transplant + com_heart_failure + com_diabetes + BMIhigh + smoke)

final <- final %>% mutate(com2 = case_when(com1 >= 2 ~ 2,
                                           com1 == 1 ~ 1,
                                           com1 == 0 ~ 0))
final <- final %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))

j <- 1
for(l in c(0:2)){
  out <- data.frame(matrix(0, 11, 10))
  colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC", "N_comorb","Ncomorb_interact")
  for(i in seq(1,11)){
    data_EUR <- final %>% filter(pop == "EUR") %>% filter(com2 == l)
    data_EUR <- data_EUR %>% filter(pop == "EUR") %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))
    LM <- glmer(as.formula(paste0(outcome[i]," ~  snp + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
    out[i,1] <- outcome[i]
    out[i,2] <- summary(LM)$coefficient[2,1]
    out[i,3] <- summary(LM)$coefficient[2,2]
    out[i,4] <- summary(LM)$coefficient[2,4]
    out[i,5] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% sum(na.rm=T)
    out[i,6] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% summarise(., sum(!is.na(.))) - out[i,5]
    out[i,7] <- snps[j]
    out[i,8] <- min(round(sum(data_EUR[snps[j]], na.rm=T)),  2*dim(data_EUR)[1] - round(sum(data_EUR[snps[j]], na.rm=T)))
    out[i,9] <- l
    data_EUR <- final %>% filter(pop == "EUR") 
    data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))
    LM <- glmer(as.formula(paste0(outcome[i]," ~  snp*com1 + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
    out[i,10] <- summary(LM)$coefficient["snp:com1",4]
  }
  print(out)
  write.table(out, file="results/clinical_EUR_comorbidities.results", append=T,quote = F, col.names = F, row.names = F, sep="\t")
}

out <- fread("results/clinical_EUR_comorbidities.results")
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC", "N_comorb","Ncomorb_interact")
out <- out %>% filter(Outcome != 0)
write.table(out, file="results/clinical_EUR_comorbidities.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")


outcome <- c("resp_severe","resp_mild", "icu_admit" ,"vte", "cardiac",
             "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")
final <- final %>% mutate_at(.vars = c(outcome), .funs = funs(ifelse(. == -1, NA, .)))
data_EUR <- final %>% filter(pop == "EUR") %>% drop_na(PC1) %>% filter(hospitalization == 1)

for(l in c(0:2)){
  out <- data.frame(matrix(0, 10, 10))
  colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC", "N_comorb","Ncomorb_interact")
  for(i in c(1:6,8:11)){
    data_EUR <- final %>% filter(pop == "EUR") %>% filter(hospitalization == 1) %>% drop_na(PC1) %>% filter(com2 == l)
    data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))
    LM <- glmer(as.formula(paste0(outcome[i]," ~  snp + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
    out[i,1] <- outcome[i]
    out[i,2] <- summary(LM)$coefficient[2,1]
    out[i,3] <- summary(LM)$coefficient[2,2]
    out[i,4] <- summary(LM)$coefficient[2,4]
    out[i,5] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% sum(na.rm=T)
    out[i,6] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% summarise(., sum(!is.na(.))) - out[i,5]
    out[i,7] <- snps[j]
    out[i,8] <- min(round(sum(data_EUR[snps[j]], na.rm=T)),  2*dim(data_EUR)[1] - round(sum(data_EUR[snps[j]], na.rm=T)))
    out[i,9] <- l
    data_EUR <- final %>% filter(pop == "EUR") %>% filter(hospitalization == 1)
    data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))
    LM <- glmer(as.formula(paste0(outcome[i]," ~  snp*com1 + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
    out[i,10] <- summary(LM)$coefficient["snp:com1",4]
  }
  print(out)
  write.table(out, file="results/clinical_EUR_comorbidities_hosp.results", append=T,quote = F, col.names = F, row.names = F, sep="\t")
}

out <- fread("results/clinical_EUR_comorbidities_hosp.results")
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC", "N_comorb","Ncomorb_interact")
out <- out %>% filter(Outcome != 0)
write.table(out, file="results/clinical_EUR_comorbidities_hosp.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")
