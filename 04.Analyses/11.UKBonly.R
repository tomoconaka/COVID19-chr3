setwd("/home/tomoko/02.associations")
outcome <- c("resp_severe","resp_mild", "icu_admit_rev" ,"vte", "cardiac",
             "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")
final <- readRDS("../01.data.QC/curated_clinical/complication_dat_final.rds")
final <- final %>% mutate(height = ifelse(height == -1, NA, height))
final <- final %>% mutate(BMI = weight/((height/100)^2))
final <- final %>% mutate_at(.vars = c(outcome), .funs = funs(ifelse(. == -1, NA, .)))
final <- final %>% mutate(age2 = age_at_diagnosis*age_at_diagnosis)
#data_EUR <- final %>% filter(age_at_diagnosis <= 60) %>% filter(pop == "EUR") %>% drop_na(PC1) #%>% filter(study == "UKB")
snps <- colnames(final)[grepl("chr", colnames(final))][3:12]

final <- final %>% mutate(smoke = case_when(smoking == 2 ~ 0,
                                            smoking == 1 ~ 1,
                                            smoking == 0 ~ 1,
                                            smoking == -1 ~ -1))
final$smoke[final$smoke == -1] <- NA
library(jtools)
library(dplyr)
library(tidyr)
final <- final %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))
final <- final %>% mutate(Male = ifelse(sex == 1, 0, 1))
final <- final %>% mutate(BMIhigh = ifelse(BMI >= 30, 1, 0))
final <- final %>% filter(study == "UKB")
outcome <- c("a1" ,"resp_severe" , "icu_admit_rev" , "vte" , "cardiac", "aki", "hepatic","hospitalization")
for(i in seq(1,8)){
  age60 <- final
  age60 <- age60 %>% filter(pop == "EUR")
  age60 <- age60[,c("age_at_diagnosis", "Male", "smoke","BMIhigh","com_cancer","com_chronic_kidney","com_chronic_pulm",
                    "com_transplant","com_heart_failure","com_diabetes", outcome[i], "snp")]
  age60 <- age60[complete.cases(age60),]
  LM1 <- glm(paste0(outcome[i]," ~ snp + age_at_diagnosis + Male + smoke + BMIhigh + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes + com_cancer"), data=age60, family=binomial(link="logit"))
  tmp <- data.frame(summ(LM1)$coeftable)
  tmp <- tmp %>% mutate(riskfactor=rownames(summ(LM1)$coeftable))
  tmp <- tmp[-1,]
  tmp <- tmp %>% mutate(OR = exp(Est.),
                        LL = exp(Est. + qnorm(0.025)*S.E.),
                        UL = exp(Est. + qnorm(0.975)*S.E.),
                        outcome = outcome[i],
                        agerange = "ALL")
  tmp <- tmp %>% mutate(freq=case_when(riskfactor == "snp" ~ sum(age60$snp == 1)/dim(age60)[1],
                                       riskfactor == "Male" ~ sum(age60$Male == 1)/dim(age60)[1],
                                       riskfactor == "smoke" ~ sum(age60$smoke == 1)/dim(age60)[1],
                                       riskfactor == "BMIhigh" ~ sum(age60$BMIhigh == 1)/dim(age60)[1],
                                       riskfactor == "com_cancer" ~ sum(age60$com_cancer == 1)/dim(age60)[1],
                                       riskfactor == "com_chronic_kidney" ~ sum(age60$com_chronic_kidney == 1)/dim(age60)[1],
                                       riskfactor == "com_chronic_pulm" ~ sum(age60$com_chronic_pulm == 1)/dim(age60)[1],
                                       riskfactor == "com_transplant" ~ sum(age60$com_transplant == 1)/dim(age60)[1],
                                       riskfactor == "com_heart_failure" ~ sum(age60$com_heart_failure == 1)/dim(age60)[1],
                                       riskfactor == "com_diabetes" ~ sum(age60$com_diabetes == 1)/dim(age60)[1]
  ))
  tmp1 <- tmp %>% select(c("outcome","agerange","riskfactor","OR","LL","UL","p","freq"))
  tmp1$Ncase <- sum(age60[,outcome[i]])
  tmp1$Ncontrol <- dim(age60)[1] - sum(age60[,outcome[i]])
  tmp1 %>% write.table("results/UKB.multivariate.table.tsv",row.names = F, col.names = F, append=T,sep="\t")
}

for(i in seq(1,8)){
  age60 <- final %>% filter(age_at_diagnosis <= 60)
  age60 <- age60 %>% filter(pop == "EUR")
  age60 <- age60[,c("age_at_diagnosis", "Male", "smoke","BMIhigh","com_cancer","com_chronic_kidney","com_chronic_pulm",
                    "com_transplant","com_heart_failure","com_diabetes", outcome[i], "snp")]
  age60 <- age60[complete.cases(age60),]
  LM1 <- glm(paste0(outcome[i]," ~ snp + age_at_diagnosis + Male + smoke + BMIhigh + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes + com_cancer"), data=age60, family=binomial(link="logit"))
  tmp <- data.frame(summ(LM1)$coeftable)
  tmp <- tmp %>% mutate(riskfactor=rownames(summ(LM1)$coeftable))
  tmp <- tmp[-1,]
  tmp <- tmp %>% mutate(OR = exp(Est.),
                        LL = exp(Est. + qnorm(0.025)*S.E.),
                        UL = exp(Est. + qnorm(0.975)*S.E.),
                        outcome = outcome[i],
                        agerange = "Age≤60")
  tmp <- tmp %>% mutate(freq=case_when(riskfactor == "snp" ~ sum(age60$snp == 1)/dim(age60)[1],
                                       riskfactor == "Male" ~ sum(age60$Male == 1)/dim(age60)[1],
                                       riskfactor == "smoke" ~ sum(age60$smoke == 1)/dim(age60)[1],
                                       riskfactor == "BMIhigh" ~ sum(age60$BMIhigh == 1)/dim(age60)[1],
                                       riskfactor == "com_cancer" ~ sum(age60$com_cancer == 1)/dim(age60)[1],
                                       riskfactor == "com_chronic_kidney" ~ sum(age60$com_chronic_kidney == 1)/dim(age60)[1],
                                       riskfactor == "com_chronic_pulm" ~ sum(age60$com_chronic_pulm == 1)/dim(age60)[1],
                                       riskfactor == "com_transplant" ~ sum(age60$com_transplant == 1)/dim(age60)[1],
                                       riskfactor == "com_heart_failure" ~ sum(age60$com_heart_failure == 1)/dim(age60)[1],
                                       riskfactor == "com_diabetes" ~ sum(age60$com_diabetes == 1)/dim(age60)[1]
  ))
  tmp1 <- tmp %>% select(c("outcome","agerange","riskfactor","OR","LL","UL","p","freq"))
  tmp1$Ncase <- sum(age60[,outcome[i]])
  tmp1$Ncontrol <- dim(age60)[1] - sum(age60[,outcome[i]])
  tmp1 %>% write.table("results/UKB.multivariate.table.tsv",row.names = F, col.names = F, append=T,sep="\t")
}

for(i in seq(1,8)){
  age60 <- final %>% filter(age_at_diagnosis > 60)
  age60 <- age60 %>% filter(pop == "EUR")
  age60 <- age60[,c("age_at_diagnosis", "Male", "smoke","BMIhigh","com_cancer","com_chronic_kidney","com_chronic_pulm",
                    "com_transplant","com_heart_failure","com_diabetes", outcome[i], "snp")]
  age60 <- age60[complete.cases(age60),]
  LM1 <- glm(paste0(outcome[i]," ~ snp + age_at_diagnosis + Male + smoke + BMIhigh + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes + com_cancer"), data=age60, family=binomial(link="logit"))
  tmp <- data.frame(summ(LM1)$coeftable)
  tmp <- tmp %>% mutate(riskfactor=rownames(summ(LM1)$coeftable))
  tmp <- tmp[-1,]
  tmp <- tmp %>% mutate(OR = exp(Est.),
                        LL = exp(Est. + qnorm(0.025)*S.E.),
                        UL = exp(Est. + qnorm(0.975)*S.E.),
                        outcome = outcome[i],
                        agerange = "Age>60")
  tmp <- tmp %>% mutate(freq=case_when(riskfactor == "snp" ~ sum(age60$snp == 1)/dim(age60)[1],
                                       riskfactor == "Male" ~ sum(age60$Male == 1)/dim(age60)[1],
                                       riskfactor == "smoke" ~ sum(age60$smoke == 1)/dim(age60)[1],
                                       riskfactor == "BMIhigh" ~ sum(age60$BMIhigh == 1)/dim(age60)[1],
                                       riskfactor == "com_cancer" ~ sum(age60$com_cancer == 1)/dim(age60)[1],
                                       riskfactor == "com_chronic_kidney" ~ sum(age60$com_chronic_kidney == 1)/dim(age60)[1],
                                       riskfactor == "com_chronic_pulm" ~ sum(age60$com_chronic_pulm == 1)/dim(age60)[1],
                                       riskfactor == "com_transplant" ~ sum(age60$com_transplant == 1)/dim(age60)[1],
                                       riskfactor == "com_heart_failure" ~ sum(age60$com_heart_failure == 1)/dim(age60)[1],
                                       riskfactor == "com_diabetes" ~ sum(age60$com_diabetes == 1)/dim(age60)[1]
  ))
  tmp1 <- tmp %>% select(c("outcome","agerange","riskfactor","OR","LL","UL","p","freq"))
  tmp1$Ncase <- sum(age60[,outcome[i]])
  tmp1$Ncontrol <- dim(age60)[1] - sum(age60[,outcome[i]])
  tmp1 %>% write.table("results/UKB.multivariate.table.tsv",row.names = F, col.names = F, append=T,sep="\t")
}

out <- fread("results/UKB.multivariate.table.tsv")
colnames(out) <- c("outcome","agerange","riskfactor","OR","LL","UL","p","freq","Ncase","Ncontrol")
out %>% write.xlsx("results/UKB.multivariate.table.xlsx",row.names = F)
out %>% write.table("results/UKB.multivariate.table.tsv",row.names = F, col.names = F, append=F,sep="\t")
