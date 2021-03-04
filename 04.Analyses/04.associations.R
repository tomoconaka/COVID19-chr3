setwd("/home/tomoko/02.associations")
final <- readRDS("../01.data.QC/curated_clinical/complication_dat_final.rds")
snps <- colnames(final)[grepl("chr", colnames(final))][3:12]
outcome <- c("resp_severe","resp_mild", "icu_admit_rev" ,"vte", "cardiac",
            "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")

final <- final %>% mutate_at(.vars = c(outcome), .funs = funs(ifelse(. == -1, NA, .)))

final <- final %>% mutate_at(.vars = c("smoke","BMI",outcome), .funs = funs(ifelse(. == "NA", NA, .)))

final <- final %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))
data_EUR <- final %>% filter(pop == "EUR")
#data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))

j <- 1
out <- data.frame(matrix(0, 11, 12))
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC","age_interact","beta_age_snp","se_age_snp","pvalue_age_snp")
for(i in c(1:11)){
  LM <- glmer(as.formula(paste0(outcome[i]," ~  snp  + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
  summary(LM)
  out[i,1] <- outcome[i]
  out[i,2] <- summary(LM)$coefficient[2,1]
  out[i,3] <- summary(LM)$coefficient[2,2]
  out[i,4] <- summary(LM)$coefficient[2,4]
  out[i,5] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% sum(na.rm=T)
  out[i,6] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% summarise(., sum(!is.na(.))) - out[i,5]
  out[i,7] <- snps[j]
  out[i,8] <- min(round(sum(data_EUR[snps[j]], na.rm=T)),  2*dim(data_EUR)[1] - round(sum(data_EUR[snps[j]], na.rm=T)))
  LM <- glmer(as.formula(paste0(outcome[i]," ~  snp*age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
  out[i,9] <- summary(LM)$coefficient["snp:age_at_diagnosis",4]
  tmp <- data_EUR %>% drop_na(outcome[i])
  tmp <- tmp[unlist(tmp[,outcome[i]]) == 1,]
  LM <- glmer(as.formula(paste0("age_at_diagnosis ~  snp + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=tmp, family ="gaussian")
  out[i,10] <- summary(LM)$coefficient["snp",1]
  out[i,11] <- summary(LM)$coefficient["snp",2]
  out[i,12] <- parameters::p_value(LM)[2,2]
}
print(out)
write.table(out, file="results/clinical_EUR.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")


data_EUR <- final %>% filter(pop == "AFR")
outcome <- c("resp_severe","resp_mild", "icu_admit_rev" ,"vte", "cardiac",
             "hospitalization", "a1","aki","hepatic","hosp_thrombo")
data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))

out <- data.frame(matrix(0, 10, 12))
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC","beta_age_snp","se_age_snp","pvalue_age_snp","caseMAC")
for(i in c(1:10)){
  LM <- glmer(as.formula(paste0(outcome[i]," ~  snp  + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
  summary(LM)
  out[i,1] <- outcome[i]
  out[i,2] <- summary(LM)$coefficient[2,1]
  out[i,3] <- summary(LM)$coefficient[2,2]
  out[i,4] <- summary(LM)$coefficient[2,4]
  out[i,5] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% sum(na.rm=T)
  out[i,6] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% summarise(., sum(!is.na(.))) - out[i,5]
  out[i,7] <- snps[j]
  out[i,8] <- min(round(sum(data_EUR[snps[j]], na.rm=T)),  2*dim(data_EUR)[1] - round(sum(data_EUR[snps[j]], na.rm=T)))
  tmp <- data_EUR %>% drop_na(outcome[i])
  tmp <- tmp[unlist(tmp[,outcome[i]]) == 1,]
  out[i,12] <- min(round(sum(tmp[snps[j]], na.rm=T)),  2*dim(tmp)[1] - round(sum(tmp[snps[j]], na.rm=T)))
  LM <- glmer(as.formula(paste0("age_at_diagnosis ~  snp + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=tmp, family ="gaussian")
  out[i,9] <- summary(LM)$coefficient["snp",1]
  out[i,10] <- summary(LM)$coefficient["snp",2]
  out[i,11] <- parameters::p_value(LM)[2,2]
  out[i,12] <- min(round(sum(tmp[snps[j]], na.rm=T)),  2*dim(tmp)[1] - round(sum(tmp[snps[j]], na.rm=T)))
}
print(out)
write.table(out, file="results/clinical_AFR.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")

data_EUR <- final %>% filter(pop == "SAS")
data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))

out <- data.frame(matrix(0, 10, 12))
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC","beta_age_snp","se_age_snp","pvalue_age_snp","caseMAC")
for(i in c(1:10)){
  LM <- glmer(as.formula(paste0(outcome[i]," ~  snp  + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
  summary(LM)
  out[i,1] <- outcome[i]
  out[i,2] <- summary(LM)$coefficient[2,1]
  out[i,3] <- summary(LM)$coefficient[2,2]
  out[i,4] <- summary(LM)$coefficient[2,4]
  out[i,5] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% sum(na.rm=T)
  out[i,6] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% summarise(., sum(!is.na(.))) - out[i,5]
  out[i,7] <- snps[j]
  out[i,8] <- min(round(sum(data_EUR[snps[j]], na.rm=T)),  2*dim(data_EUR)[1] - round(sum(data_EUR[snps[j]], na.rm=T)))
  tmp <- data_EUR %>% drop_na(outcome[i])
  tmp <- tmp[unlist(tmp[,outcome[i]]) == 1,]
  out[i,12] <- min(round(sum(tmp[snps[j]], na.rm=T)),  2*dim(tmp)[1] - round(sum(tmp[snps[j]], na.rm=T)))
  LM <- glmer(as.formula(paste0("age_at_diagnosis ~  snp + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=tmp, family ="gaussian")
  out[i,9] <- summary(LM)$coefficient["snp",1]
  out[i,10] <- summary(LM)$coefficient["snp",2]
  out[i,11] <- parameters::p_value(LM)[2,2]
}
print(out)
write.table(out, file="results/clinical_SAS.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")

data_EUR <- final %>% filter(pop == "EAS")
data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))

out <- data.frame(matrix(0, 10, 12))
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC","beta_age_snp","se_age_snp","pvalue_age_snp","caseMAC")
for(i in c(1:10)){
  out[i,1] <- outcome[i]
  out[i,5] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% sum(na.rm=T)
  out[i,6] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% summarise(., sum(!is.na(.))) - out[i,5]
  out[i,7] <- snps[j]
  out[i,8] <- min(round(sum(data_EUR[snps[j]], na.rm=T)),  2*dim(data_EUR)[1] - round(sum(data_EUR[snps[j]], na.rm=T)))
  LM <- glmer(as.formula(paste0(outcome[i]," ~  snp  + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
  summary(LM)
  out[i,2] <- summary(LM)$coefficient[2,1]
  out[i,3] <- summary(LM)$coefficient[2,2]
  out[i,4] <- summary(LM)$coefficient[2,4]
  tmp <- data_EUR %>% drop_na(outcome[i])
  tmp <- tmp[unlist(tmp[,outcome[i]]) == 1,]
  out[i,12] <- min(round(sum(tmp[snps[j]], na.rm=T)),  2*dim(tmp)[1] - round(sum(tmp[snps[j]], na.rm=T)))
  LM <- glmer(as.formula(paste0("age_at_diagnosis ~  snp + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=tmp, family ="gaussian")
  out[i,9] <- summary(LM)$coefficient["snp",1]
  out[i,10] <- summary(LM)$coefficient["snp",2]
  out[i,11] <- parameters::p_value(LM)[2,2]
}
print(out)
write.table(out, file="results/clinical_EAS.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")

data_EUR <- final %>% filter(pop == "AMR")

data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))

out <- data.frame(matrix(0, 10, 12))
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC","beta_age_snp","se_age_snp","pvalue_age_snp","caseMAC")
for(i in c(1:10)){
  out[i,1] <- outcome[i]
  out[i,5] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% sum(na.rm=T)
  out[i,6] <- data_EUR %>% drop_na(snps[j]) %>% drop_na(PC1) %>% dplyr::select(c(outcome[i])) %>% summarise(., sum(!is.na(.))) - out[i,5]
  out[i,7] <- snps[j]
  out[i,8] <- min(round(sum(data_EUR[snps[j]], na.rm=T)),  2*dim(data_EUR)[1] - round(sum(data_EUR[snps[j]], na.rm=T)))
  LM <- glmer(as.formula(paste0(outcome[i]," ~  snp  + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=data_EUR, family ="binomial")
  summary(LM)
  out[i,2] <- summary(LM)$coefficient[2,1]
  out[i,3] <- summary(LM)$coefficient[2,2]
  out[i,4] <- summary(LM)$coefficient[2,4]
  tmp <- data_EUR %>% drop_na(outcome[i])
  tmp <- tmp[unlist(tmp[,outcome[i]]) == 1,]
  out[i,12] <- min(round(sum(tmp[snps[j]], na.rm=T)),  2*dim(tmp)[1] - round(sum(tmp[snps[j]], na.rm=T)))
  LM <- glmer(as.formula(paste0("age_at_diagnosis ~  snp + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")), data=tmp, family ="gaussian")
  out[i,9] <- summary(LM)$coefficient["snp",1]
  out[i,10] <- summary(LM)$coefficient["snp",2]
  out[i,11] <- parameters::p_value(LM)[2,2]
}
print(out)
write.table(out, file="results/clinical_AMR.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")
