setwd("/home/tomoko/02.associations")
final <- readRDS("../01.data.QC/curated_clinical/complication_dat_final.rds")
snps <- colnames(final)[grepl("chr", colnames(final))][3:12]
outcome <- c("resp_severe","resp_mild", "icu_admit_rev" ,"vte", "cardiac",
             "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")

final <- final %>% mutate(height = ifelse(height == -1, NA, height))
final <- final %>% mutate(BMI = weight/((height/100)^2))
final <- final %>% mutate_at(.vars = c(outcome), .funs = funs(ifelse(. == -1, NA, .)))
final <- final %>% mutate(age2 = age_at_diagnosis*age_at_diagnosis)
final <- final %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))
data_EUR <- final %>% filter(age_at_diagnosis <= 60) %>% filter(pop == "EUR")# %>% drop_na(PC1) #%>% filter(study == "UKB")
data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))

j <- 1
out <- data.frame(matrix(0, 11, 8))
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC")
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
}
print(out)
write.table(out, file="results/clinical_EUR_age_young60.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")

data_EUR <- final %>% filter(age_at_diagnosis > 60) %>% filter(pop == "EUR")# %>% drop_na(PC1)
data_EUR <- data_EUR %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5"), .funs = funs(scale(.)))
out <- data.frame(matrix(0, 11, 8))
colnames(out) <- c("Outcome", "beta","se","pvalue", "N_case", "N_control", "snp", "MAC")
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
}
print(out)

write.table(out, file="results/clinical_EUR_age_old60.results", append=F,quote = F, col.names = T, row.names = F, sep="\t")
