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
  tmp1 %>% write.table("results/multivariate.table.tsv",row.names = F, col.names = F, append=T,sep="\t")
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
  tmp1 %>% write.table("results/multivariate.table.tsv",row.names = F, col.names = F, append=T,sep="\t")
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
  tmp1 %>% write.table("results/multivariate.table.tsv",row.names = F, col.names = F, append=T,sep="\t")
}

out <- fread("results/multivariate.table.tsv")
colnames(out) <- c("outcome","agerange","riskfactor","OR","LL","UL","p","freq","Ncase","Ncontrol")
out %>% write.xlsx("results/multivariate.table.xlsx",row.names = F)
library(Hmisc)
reclassification_rev <- function (data, cOutcome, predrisk1, predrisk2, cutoff) 
{
  c1 <- cut(predrisk1, breaks = cutoff, include.lowest = TRUE, 
            right = FALSE)
  c2 <- cut(predrisk2, breaks = cutoff, include.lowest = TRUE, 
            right = FALSE)
  tabReclas <- table(`Initial Model` = c1, `Updated Model` = c2)
  cat(" _________________________________________\n")
  cat(" \n     Reclassification table    \n")
  cat(" _________________________________________\n")
  ta <- table(c1, c2, data[, cOutcome])
  cat("\n Outcome: absent \n  \n")
  TabAbs <- ta[, , 1]
  tab1 <- cbind(TabAbs, ` % reclassified` = round((rowSums(TabAbs) - 
                                                     diag(TabAbs))/rowSums(TabAbs), 2) * 100)
  names(dimnames(tab1)) <- c("Initial Model", "Updated Model")
  print(tab1)
  cat("\n \n Outcome: present \n  \n")
  TabPre <- ta[, , 2]
  tab2 <- cbind(TabPre, ` % reclassified` = round((rowSums(TabPre) - 
                                                     diag(TabPre))/rowSums(TabPre), 2) * 100)
  names(dimnames(tab2)) <- c("Initial Model", "Updated Model")
  print(tab2)
  cat("\n \n Combined Data \n  \n")
  Tab <- tabReclas
  tab <- cbind(Tab, ` % reclassified` = round((rowSums(Tab) - 
                                                 diag(Tab))/rowSums(Tab), 2) * 100)
  names(dimnames(tab)) <- c("Initial Model", "Updated Model")
  print(tab)
  cat(" _________________________________________\n")
  c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
  c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
  y <- improveProb(x1 = predrisk1, x2 = predrisk2, y = data[, 
                                                            cOutcome])
  return(y)}

out <- data.frame(matrix(0, 12, 15))
colnames(out) <- c("outcome","agerange","model","AUC","AUC.LL","AUC.UL","AUC.p-value","NRI","NRI.LL", "NRI.UL","NRI.p-value","IDI","IDI.LL","IDI.UL","IDI.p-value")
for(i in seq(1,8)){
  age60 <- final
  age60 <- age60 %>% filter(pop == "EUR")
  age60 <- age60[,c("age_at_diagnosis", "sex", "smoke","BMI","com_cancer","com_chronic_kidney","com_chronic_pulm",
                    "com_transplant","com_heart_failure","com_diabetes", outcome[i], "chr3:45823240:T:C_C")]
  age60 <- age60[complete.cases(age60),]
  LM1 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex"), data=age60, family=binomial(link="logit"))
  LM2 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex + smoke"), data=age60, family=binomial(link="logit"))
  LM3 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex + BMI "), data=age60, family=binomial(link="logit"))
  LM4 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_cancer"), data=age60, family=binomial(link="logit"))
  LM5 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_chronic_pulm"), data=age60, family=binomial(link="logit"))
  LM6 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_chronic_kidney"), data=age60, family=binomial(link="logit"))
  LM7 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_transplant"), data=age60, family=binomial(link="logit"))
  LM8 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_heart_failure"), data=age60, family=binomial(link="logit"))
  LM9 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_diabetes"), data=age60, family=binomial(link="logit"))
  LM10 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + `chr3:45823240:T:C_C`"), data=age60, family=binomial(link="logit"))
  LM11 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + smoke + BMI + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes"), data=age60, family=binomial(link="logit"))
  LM12 <- glm(paste0(outcome[i]," ~ `chr3:45823240:T:C_C` + age_at_diagnosis + sex + smoke + BMI + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes"), data=age60, family=binomial(link="logit"))
  age60$pred1 <- predict(LM1, age60)
  age60$pred2 <- predict(LM2, age60)
  age60$pred3 <- predict(LM3, age60)
  age60$pred4 <- predict(LM4, age60)
  age60$pred5 <- predict(LM5, age60)
  age60$pred6 <- predict(LM6, age60)
  age60$pred7 <- predict(LM7, age60)
  age60$pred8 <- predict(LM8, age60)
  age60$pred9 <- predict(LM9, age60)
  age60$pred10 <- predict(LM10, age60)
  age60$pred11 <- predict(LM11, age60)
  age60$pred12 <- predict(LM12, age60)
  age60 <- data.frame(age60)
  riskmodel1 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2), cGenPreds=0, cGenPredsCat=0)
  riskmodel2 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,3), cGenPreds=0, cGenPredsCat=0)
  riskmodel3 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2), cGenPreds=0, cGenPredsCat=0)
  riskmodel4 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,5), cGenPreds=0, cGenPredsCat=0)
  riskmodel5 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,6), cGenPreds=0, cGenPredsCat=0)
  riskmodel6 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,7), cGenPreds=0, cGenPredsCat=0)
  riskmodel7 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,8), cGenPreds=0, cGenPredsCat=0)
  riskmodel8 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,9), cGenPreds=0, cGenPredsCat=0)
  riskmodel9 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,10), cGenPreds=0, cGenPredsCat=0)
  riskmodel10 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2), cGenPreds=12, cGenPredsCat=0)
  riskmodel11 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2,3,5:10), cGenPreds=0, cGenPredsCat=0)
  riskmodel12 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2,3,5:10), cGenPreds=12, cGenPredsCat=0)
  predRisk1 <- predRisk(riskModel=riskmodel1)
  predRisk2 <- predRisk(riskmodel2)
  predRisk3 <- predRisk(riskmodel3)
  predRisk4 <- predRisk(riskmodel4)
  predRisk5 <- predRisk(riskmodel5)
  predRisk6 <- predRisk(riskmodel6)
  predRisk7 <- predRisk(riskmodel7)
  predRisk8 <- predRisk(riskmodel8)
  predRisk9 <- predRisk(riskmodel9)
  predRisk10 <- predRisk(riskmodel10)
  predRisk11 <- predRisk(riskmodel11)
  predRisk12 <- predRisk(riskmodel12)
  value<-as.vector(quantile(predRisk2))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk2, cutoff)
  response <- as.numeric(unlist(age60[,outcome[i]]))
  ROC1 <- roc(response=response, predictor=age60$pred1) 
  ROC2 <- roc(response=response, predictor=age60$pred2) 
  ROC3 <- roc(response=response, predictor=age60$pred3) 
  ROC4 <- roc(response=response, predictor=age60$pred4)
  ROC5 <- roc(response=response, predictor=age60$pred5)
  ROC6 <- roc(response=response, predictor=age60$pred6)
  ROC7 <- roc(response=response, predictor=age60$pred7)
  ROC8 <- roc(response=response, predictor=age60$pred8)
  ROC9 <- roc(response=response, predictor=age60$pred9)
  ROC10 <- roc(response=response, predictor=age60$pred10)
  ROC11 <- roc(response=response, predictor=age60$pred11)
  ROC12 <- roc(response=response, predictor=age60$pred12)
  out[1,1] <- outcome[i]
  out[1,2] <- "ALL"
  out[1,3] <- "age+sex"
  out[1,4] <- ROC1$auc
  out[1,5] <- ci.auc(ROC1, conf.level=0.95)[1]
  out[1,6] <- ci.auc(ROC1, conf.level=0.95)[3]
  out[2,1] <- outcome[i]
  out[2,2] <- "ALL"
  out[2,3] <- "age+sex+smoking"
  out[2,4] <- ROC2$auc
  out[2,5] <- ci.auc(ROC2, conf.level=0.95)[1]
  out[2,6] <- ci.auc(ROC2, conf.level=0.95)[3]
  out[2,7] <- roc.test(ROC1, ROC2, method="delong", alternative="two.sided")$p
  out[2,8] <- y$nri
  out[2,9] <- y$nri -  1.96 * y$se.nri
  out[2,10] <- y$nri + 1.96 * y$se.nri
  out[2,11] <- pnorm(-abs(y$z.nri))
  out[2,12] <- y$idi
  out[2,13] <- y$idi -  1.96 * y$se.idi
  out[2,14] <- y$idi + 1.96 * y$se.idi
  out[2,15] <- pnorm(-abs(y$z.idi))
  out[3,1] <- outcome[i]
  out[3,2] <- "ALL"
  out[3,3] <- "age+sex+BMI"
  out[3,4] <- ROC3$auc
  out[3,5] <- ci.auc(ROC3, conf.level=0.95)[1]
  out[3,6] <- ci.auc(ROC3, conf.level=0.95)[3]
  out[3,7] <- roc.test(ROC1, ROC3, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk3))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk3, cutoff)
  out[3,8] <- y$nri
  out[3,9] <- y$nri -  1.96 * y$se.nri
  out[3,10] <- y$nri + 1.96 * y$se.nri
  out[3,11] <- pnorm(-abs(y$z.nri))
  out[3,12] <- y$idi
  out[3,13] <- y$idi -  1.96 * y$se.idi
  out[3,14] <- y$idi + 1.96 * y$se.idi
  out[3,15] <- pnorm(-abs(y$z.idi))
  out[4,1] <- outcome[i]
  out[4,2] <- "ALL"
  out[4,3] <- "age+sex+cancer"
  out[4,4] <- ROC4$auc
  out[4,5] <- ci.auc(ROC4, conf.level=0.95)[1]
  out[4,6] <- ci.auc(ROC4, conf.level=0.95)[3]
  out[4,7] <- roc.test(ROC1, ROC4, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk4))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk4, cutoff)
  out[4,8] <- y$nri
  out[4,9] <- y$nri -  1.96 * y$se.nri
  out[4,10] <- y$nri + 1.96 * y$se.nri
  out[4,11] <- pnorm(-abs(y$z.nri))
  out[4,12] <- y$idi
  out[4,13] <- y$idi -  1.96 * y$se.idi
  out[4,14] <- y$idi + 1.96 * y$se.idi
  out[4,15] <- pnorm(-abs(y$z.idi))
  out[5,1] <- outcome[i]
  out[5,2] <- "ALL"
  out[5,3] <- "age+sex+CKD"
  out[5,4] <- ROC5$auc
  out[5,5] <- ci.auc(ROC5, conf.level=0.95)[1]
  out[5,6] <- ci.auc(ROC5, conf.level=0.95)[3]
  out[5,7] <- roc.test(ROC1, ROC5, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk5))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk5, cutoff)
  out[5,8] <- y$nri
  out[5,9] <- y$nri -  1.96 * y$se.nri
  out[5,10] <- y$nri + 1.96 * y$se.nri
  out[5,11] <- pnorm(-abs(y$z.nri))
  out[5,12] <- y$idi
  out[5,13] <- y$idi -  1.96 * y$se.idi
  out[5,14] <- y$idi + 1.96 * y$se.idi
  out[5,15] <- pnorm(-abs(y$z.idi))
  out[6,1] <- outcome[i]
  out[6,2] <- "ALL"
  out[6,3] <- "age+sex+COPD"
  out[6,4] <- ROC6$auc
  out[6,5] <- ci.auc(ROC6, conf.level=0.95)[1]
  out[6,6] <- ci.auc(ROC6, conf.level=0.95)[3]
  out[6,7] <- roc.test(ROC1, ROC6, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk6))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk6, cutoff)
  out[6,8] <- y$nri
  out[6,9] <- y$nri -  1.96 * y$se.nri
  out[6,10] <- y$nri + 1.96 * y$se.nri
  out[6,11] <- pnorm(-abs(y$z.nri))
  out[6,12] <- y$idi
  out[6,13] <- y$idi -  1.96 * y$se.idi
  out[6,14] <- y$idi + 1.96 * y$se.idi
  out[6,15] <- pnorm(-abs(y$z.idi))
  out[7,1] <- outcome[i]
  out[7,2] <- "ALL"
  out[7,3] <- "age+sex+transplantation"
  out[7,4] <- ROC7$auc
  out[7,5] <- ci.auc(ROC7, conf.level=0.95)[1]
  out[7,6] <- ci.auc(ROC7, conf.level=0.95)[3]
  out[7,7] <- roc.test(ROC1, ROC7, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk7))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk7, cutoff)
  out[7,8] <- y$nri
  out[7,9] <- y$nri -  1.96 * y$se.nri
  out[7,10] <- y$nri + 1.96 * y$se.nri
  out[7,11] <- pnorm(-abs(y$z.nri))
  out[7,12] <- y$idi
  out[7,13] <- y$idi -  1.96 * y$se.idi
  out[7,14] <- y$idi + 1.96 * y$se.idi
  out[7,15] <- pnorm(-abs(y$z.idi))
  out[8,1] <- outcome[i]
  out[8,2] <- "ALL"
  out[8,3] <- "age+sex+CHF"
  out[8,4] <- ROC8$auc
  out[8,5] <- ci.auc(ROC8, conf.level=0.95)[1]
  out[8,6] <- ci.auc(ROC8, conf.level=0.95)[3]
  out[8,7] <- roc.test(ROC1, ROC8, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk8))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk8, cutoff)
  out[8,8] <- y$nri
  out[8,9] <- y$nri -  1.96 * y$se.nri
  out[8,10] <- y$nri + 1.96 * y$se.nri
  out[8,11] <- pnorm(-abs(y$z.nri))
  out[8,12] <- y$idi
  out[8,13] <- y$idi -  1.96 * y$se.idi
  out[8,14] <- y$idi + 1.96 * y$se.idi
  out[8,15] <- pnorm(-abs(y$z.idi))
  out[9,1] <- outcome[i]
  out[9,2] <- "ALL"
  out[9,3] <- "age+sex+diabetes"
  out[9,4] <- ROC9$auc
  out[9,5] <- ci.auc(ROC9, conf.level=0.95)[1]
  out[9,6] <- ci.auc(ROC9, conf.level=0.95)[3]
  out[9,7] <- roc.test(ROC1, ROC9, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk9))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk9, cutoff)
  out[9,8] <- y$nri
  out[9,9] <- y$nri -  1.96 * y$se.nri
  out[9,10] <- y$nri + 1.96 * y$se.nri
  out[9,11] <- pnorm(-abs(y$z.nri))
  out[9,12] <- y$idi
  out[9,13] <- y$idi -  1.96 * y$se.idi
  out[9,14] <- y$idi + 1.96 * y$se.idi
  out[9,15] <- pnorm(-abs(y$z.idi))
  out[10,1] <- outcome[i]
  out[10,2] <- "ALL"
  out[10,3] <- "age+sex+rs10490770"
  out[10,4] <- ROC10$auc
  out[10,5] <- ci.auc(ROC10, conf.level=0.95)[1]
  out[10,6] <- ci.auc(ROC10, conf.level=0.95)[3]
  out[10,7] <- roc.test(ROC1, ROC10, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk10))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk10, cutoff)
  out[10,8] <- y$nri
  out[10,9] <- y$nri -  1.96 * y$se.nri
  out[10,10] <- y$nri + 1.96 * y$se.nri
  out[10,11] <- pnorm(-abs(y$z.nri))
  out[10,12] <- y$idi
  out[10,13] <- y$idi -  1.96 * y$se.idi
  out[10,14] <- y$idi + 1.96 * y$se.idi
  out[10,15] <- pnorm(-abs(y$z.idi))
  out[11,1] <- outcome[i]
  out[11,2] <- "ALL"
  out[11,3] <- "age+sex+allriskfactors"
  out[11,4] <- ROC11$auc
  out[11,5] <- ci.auc(ROC11, conf.level=0.95)[1]
  out[11,6] <- ci.auc(ROC11, conf.level=0.95)[3]
  out[11,7] <- roc.test(ROC1, ROC11, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk11))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk11, cutoff)
  out[11,8] <- y$nri
  out[11,9] <- y$nri -  1.96 * y$se.nri
  out[11,10] <- y$nri + 1.96 * y$se.nri
  out[11,11] <- pnorm(-abs(y$z.nri))
  out[11,12] <- y$idi
  out[11,13] <- y$idi -  1.96 * y$se.idi
  out[11,14] <- y$idi + 1.96 * y$se.idi
  out[11,15] <- pnorm(-abs(y$z.idi))
  out[12,1] <- outcome[i]
  out[12,2] <- "ALL"
  out[12,3] <- "age+sex+allriskfactors+rs10490770"
  out[12,4] <- ROC12$auc
  out[12,5] <- ci.auc(ROC12, conf.level=0.95)[1]
  out[12,6] <- ci.auc(ROC12, conf.level=0.95)[3]
  out[12,7] <- roc.test(ROC11, ROC12, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk12))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk11, predrisk2=predRisk12, cutoff)
  out[12,8] <- y$nri
  out[12,9] <- y$nri -  1.96 * y$se.nri
  out[12,10] <- y$nri + 1.96 * y$se.nri
  out[12,11] <- pnorm(-abs(y$z.nri))
  out[12,12] <- y$idi
  out[12,13] <- y$idi -  1.96 * y$se.idi
  out[12,14] <- y$idi + 1.96 * y$se.idi
  out[12,15] <- pnorm(-abs(y$z.idi))
  out %>% write.table("results/AUC.tsv",row.names = F, col.names = F, append=T,sep="\t")
}

out <- data.frame(matrix(0, 12, 15))
colnames(out) <- c("outcome","agerange","model","AUC","AUC.LL","AUC.UL","AUC.p-value","NRI","NRI.LL", "NRI.UL","NRI.p-value","IDI","IDI.LL","IDI.UL","IDI.p-value")
for(i in seq(1,8)){
  age60 <- final
  age60 <- age60 %>% filter(pop == "EUR") %>% filter(age_at_diagnosis <= 60)
  age60 <- age60[,c("age_at_diagnosis", "sex", "smoke","BMI","com_cancer","com_chronic_kidney","com_chronic_pulm",
                    "com_transplant","com_heart_failure","com_diabetes", outcome[i], "chr3:45823240:T:C_C")]
  age60 <- age60[complete.cases(age60),]
  LM1 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex"), data=age60, family=binomial(link="logit"))
  LM2 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex + smoke"), data=age60, family=binomial(link="logit"))
  LM3 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex + BMI "), data=age60, family=binomial(link="logit"))
  LM4 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_cancer"), data=age60, family=binomial(link="logit"))
  LM5 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_chronic_pulm"), data=age60, family=binomial(link="logit"))
  LM6 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_chronic_kidney"), data=age60, family=binomial(link="logit"))
  LM7 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_transplant"), data=age60, family=binomial(link="logit"))
  LM8 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_heart_failure"), data=age60, family=binomial(link="logit"))
  LM9 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_diabetes"), data=age60, family=binomial(link="logit"))
  LM10 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + `chr3:45823240:T:C_C`"), data=age60, family=binomial(link="logit"))
  LM11 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + smoke + BMI + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes"), data=age60, family=binomial(link="logit"))
  LM12 <- glm(paste0(outcome[i]," ~ `chr3:45823240:T:C_C` + age_at_diagnosis + sex + smoke + BMI + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes"), data=age60, family=binomial(link="logit"))
  age60$pred1 <- predict(LM1, age60)
  age60$pred2 <- predict(LM2, age60)
  age60$pred3 <- predict(LM3, age60)
  age60$pred4 <- predict(LM4, age60)
  age60$pred5 <- predict(LM5, age60)
  age60$pred6 <- predict(LM6, age60)
  age60$pred7 <- predict(LM7, age60)
  age60$pred8 <- predict(LM8, age60)
  age60$pred9 <- predict(LM9, age60)
  age60$pred10 <- predict(LM10, age60)
  age60$pred11 <- predict(LM11, age60)
  age60$pred12 <- predict(LM12, age60)
  age60 <- data.frame(age60)
  riskmodel1 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2), cGenPreds=0, cGenPredsCat=0)
  riskmodel2 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,3), cGenPreds=0, cGenPredsCat=0)
  riskmodel3 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2), cGenPreds=0, cGenPredsCat=0)
  riskmodel4 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,5), cGenPreds=0, cGenPredsCat=0)
  riskmodel5 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,6), cGenPreds=0, cGenPredsCat=0)
  riskmodel6 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,7), cGenPreds=0, cGenPredsCat=0)
  riskmodel7 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,8), cGenPreds=0, cGenPredsCat=0)
  riskmodel8 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,9), cGenPreds=0, cGenPredsCat=0)
  riskmodel9 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,10), cGenPreds=0, cGenPredsCat=0)
  riskmodel10 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2), cGenPreds=12, cGenPredsCat=0)
  riskmodel11 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2,3,5:10), cGenPreds=0, cGenPredsCat=0)
  riskmodel12 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2,3,5:10), cGenPreds=12, cGenPredsCat=0)
  predRisk1 <- predRisk(riskModel=riskmodel1)
  predRisk2 <- predRisk(riskmodel2)
  predRisk3 <- predRisk(riskmodel3)
  predRisk4 <- predRisk(riskmodel4)
  predRisk5 <- predRisk(riskmodel5)
  predRisk6 <- predRisk(riskmodel6)
  predRisk7 <- predRisk(riskmodel7)
  predRisk8 <- predRisk(riskmodel8)
  predRisk9 <- predRisk(riskmodel9)
  predRisk10 <- predRisk(riskmodel10)
  predRisk11 <- predRisk(riskmodel11)
  predRisk12 <- predRisk(riskmodel12)
  value<-as.vector(quantile(predRisk2))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk2, cutoff)
  response <- as.numeric(unlist(age60[,outcome[i]]))
  ROC1 <- roc(response=response, predictor=age60$pred1) 
  ROC2 <- roc(response=response, predictor=age60$pred2) 
  ROC3 <- roc(response=response, predictor=age60$pred3) 
  ROC4 <- roc(response=response, predictor=age60$pred4)
  ROC5 <- roc(response=response, predictor=age60$pred5)
  ROC6 <- roc(response=response, predictor=age60$pred6)
  ROC7 <- roc(response=response, predictor=age60$pred7)
  ROC8 <- roc(response=response, predictor=age60$pred8)
  ROC9 <- roc(response=response, predictor=age60$pred9)
  ROC10 <- roc(response=response, predictor=age60$pred10)
  ROC11 <- roc(response=response, predictor=age60$pred11)
  ROC12 <- roc(response=response, predictor=age60$pred12)
  out[1,1] <- outcome[i]
  out[1,2] <- "Age≤60"
  out[1,3] <- "age+sex"
  out[1,4] <- ROC1$auc
  out[1,5] <- ci.auc(ROC1, conf.level=0.95)[1]
  out[1,6] <- ci.auc(ROC1, conf.level=0.95)[3]
  out[2,1] <- outcome[i]
  out[2,2] <- "Age≤60"
  out[2,3] <- "age+sex+smoking"
  out[2,4] <- ROC2$auc
  out[2,5] <- ci.auc(ROC2, conf.level=0.95)[1]
  out[2,6] <- ci.auc(ROC2, conf.level=0.95)[3]
  out[2,7] <- roc.test(ROC1, ROC2, method="delong", alternative="two.sided")$p
  out[2,8] <- y$nri
  out[2,9] <- y$nri -  1.96 * y$se.nri
  out[2,10] <- y$nri + 1.96 * y$se.nri
  out[2,11] <- pnorm(-abs(y$z.nri))
  out[2,12] <- y$idi
  out[2,13] <- y$idi -  1.96 * y$se.idi
  out[2,14] <- y$idi + 1.96 * y$se.idi
  out[2,15] <- pnorm(-abs(y$z.idi))
  out[3,1] <- outcome[i]
  out[3,2] <- "Age≤60"
  out[3,3] <- "age+sex+BMI"
  out[3,4] <- ROC3$auc
  out[3,5] <- ci.auc(ROC3, conf.level=0.95)[1]
  out[3,6] <- ci.auc(ROC3, conf.level=0.95)[3]
  out[3,7] <- roc.test(ROC1, ROC3, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk3))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk3, cutoff)
  out[3,8] <- y$nri
  out[3,9] <- y$nri -  1.96 * y$se.nri
  out[3,10] <- y$nri + 1.96 * y$se.nri
  out[3,11] <- pnorm(-abs(y$z.nri))
  out[3,12] <- y$idi
  out[3,13] <- y$idi -  1.96 * y$se.idi
  out[3,14] <- y$idi + 1.96 * y$se.idi
  out[3,15] <- pnorm(-abs(y$z.idi))
  out[4,1] <- outcome[i]
  out[4,2] <- "Age≤60"
  out[4,3] <- "age+sex+cancer"
  out[4,4] <- ROC4$auc
  out[4,5] <- ci.auc(ROC4, conf.level=0.95)[1]
  out[4,6] <- ci.auc(ROC4, conf.level=0.95)[3]
  out[4,7] <- roc.test(ROC1, ROC4, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk4))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk4, cutoff)
  out[4,8] <- y$nri
  out[4,9] <- y$nri -  1.96 * y$se.nri
  out[4,10] <- y$nri + 1.96 * y$se.nri
  out[4,11] <- pnorm(-abs(y$z.nri))
  out[4,12] <- y$idi
  out[4,13] <- y$idi -  1.96 * y$se.idi
  out[4,14] <- y$idi + 1.96 * y$se.idi
  out[4,15] <- pnorm(-abs(y$z.idi))
  out[5,1] <- outcome[i]
  out[5,2] <- "Age≤60"
  out[5,3] <- "age+sex+CKD"
  out[5,4] <- ROC5$auc
  out[5,5] <- ci.auc(ROC5, conf.level=0.95)[1]
  out[5,6] <- ci.auc(ROC5, conf.level=0.95)[3]
  out[5,7] <- roc.test(ROC1, ROC5, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk5))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk5, cutoff)
  out[5,8] <- y$nri
  out[5,9] <- y$nri -  1.96 * y$se.nri
  out[5,10] <- y$nri + 1.96 * y$se.nri
  out[5,11] <- pnorm(-abs(y$z.nri))
  out[5,12] <- y$idi
  out[5,13] <- y$idi -  1.96 * y$se.idi
  out[5,14] <- y$idi + 1.96 * y$se.idi
  out[5,15] <- pnorm(-abs(y$z.idi))
  out[6,1] <- outcome[i]
  out[6,2] <- "Age≤60"
  out[6,3] <- "age+sex+COPD"
  out[6,4] <- ROC6$auc
  out[6,5] <- ci.auc(ROC6, conf.level=0.95)[1]
  out[6,6] <- ci.auc(ROC6, conf.level=0.95)[3]
  out[6,7] <- roc.test(ROC1, ROC6, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk6))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk6, cutoff)
  out[6,8] <- y$nri
  out[6,9] <- y$nri -  1.96 * y$se.nri
  out[6,10] <- y$nri + 1.96 * y$se.nri
  out[6,11] <- pnorm(-abs(y$z.nri))
  out[6,12] <- y$idi
  out[6,13] <- y$idi -  1.96 * y$se.idi
  out[6,14] <- y$idi + 1.96 * y$se.idi
  out[6,15] <- pnorm(-abs(y$z.idi))
  out[7,1] <- outcome[i]
  out[7,2] <- "Age≤60"
  out[7,3] <- "age+sex+transplantation"
  out[7,4] <- ROC7$auc
  out[7,5] <- ci.auc(ROC7, conf.level=0.95)[1]
  out[7,6] <- ci.auc(ROC7, conf.level=0.95)[3]
  out[7,7] <- roc.test(ROC1, ROC7, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk7))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk7, cutoff)
  out[7,8] <- y$nri
  out[7,9] <- y$nri -  1.96 * y$se.nri
  out[7,10] <- y$nri + 1.96 * y$se.nri
  out[7,11] <- pnorm(-abs(y$z.nri))
  out[7,12] <- y$idi
  out[7,13] <- y$idi -  1.96 * y$se.idi
  out[7,14] <- y$idi + 1.96 * y$se.idi
  out[7,15] <- pnorm(-abs(y$z.idi))
  out[8,1] <- outcome[i]
  out[8,2] <- "Age≤60"
  out[8,3] <- "age+sex+CHF"
  out[8,4] <- ROC8$auc
  out[8,5] <- ci.auc(ROC8, conf.level=0.95)[1]
  out[8,6] <- ci.auc(ROC8, conf.level=0.95)[3]
  out[8,7] <- roc.test(ROC1, ROC8, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk8))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk8, cutoff)
  out[8,8] <- y$nri
  out[8,9] <- y$nri -  1.96 * y$se.nri
  out[8,10] <- y$nri + 1.96 * y$se.nri
  out[8,11] <- pnorm(-abs(y$z.nri))
  out[8,12] <- y$idi
  out[8,13] <- y$idi -  1.96 * y$se.idi
  out[8,14] <- y$idi + 1.96 * y$se.idi
  out[8,15] <- pnorm(-abs(y$z.idi))
  out[9,1] <- outcome[i]
  out[9,2] <- "Age≤60"
  out[9,3] <- "age+sex+diabetes"
  out[9,4] <- ROC9$auc
  out[9,5] <- ci.auc(ROC9, conf.level=0.95)[1]
  out[9,6] <- ci.auc(ROC9, conf.level=0.95)[3]
  out[9,7] <- roc.test(ROC1, ROC9, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk9))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk9, cutoff)
  out[9,8] <- y$nri
  out[9,9] <- y$nri -  1.96 * y$se.nri
  out[9,10] <- y$nri + 1.96 * y$se.nri
  out[9,11] <- pnorm(-abs(y$z.nri))
  out[9,12] <- y$idi
  out[9,13] <- y$idi -  1.96 * y$se.idi
  out[9,14] <- y$idi + 1.96 * y$se.idi
  out[9,15] <- pnorm(-abs(y$z.idi))
  out[10,1] <- outcome[i]
  out[10,2] <- "Age≤60"
  out[10,3] <- "age+sex+rs10490770"
  out[10,4] <- ROC10$auc
  out[10,5] <- ci.auc(ROC10, conf.level=0.95)[1]
  out[10,6] <- ci.auc(ROC10, conf.level=0.95)[3]
  out[10,7] <- roc.test(ROC1, ROC10, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk10))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk10, cutoff)
  out[10,8] <- y$nri
  out[10,9] <- y$nri -  1.96 * y$se.nri
  out[10,10] <- y$nri + 1.96 * y$se.nri
  out[10,11] <- pnorm(-abs(y$z.nri))
  out[10,12] <- y$idi
  out[10,13] <- y$idi -  1.96 * y$se.idi
  out[10,14] <- y$idi + 1.96 * y$se.idi
  out[10,15] <- pnorm(-abs(y$z.idi))
  out[11,1] <- outcome[i]
  out[11,2] <- "Age≤60"
  out[11,3] <- "age+sex+allriskfactors"
  out[11,4] <- ROC11$auc
  out[11,5] <- ci.auc(ROC11, conf.level=0.95)[1]
  out[11,6] <- ci.auc(ROC11, conf.level=0.95)[3]
  out[11,7] <- roc.test(ROC1, ROC11, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk11))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk11, cutoff)
  out[11,8] <- y$nri
  out[11,9] <- y$nri -  1.96 * y$se.nri
  out[11,10] <- y$nri + 1.96 * y$se.nri
  out[11,11] <- pnorm(-abs(y$z.nri))
  out[11,12] <- y$idi
  out[11,13] <- y$idi -  1.96 * y$se.idi
  out[11,14] <- y$idi + 1.96 * y$se.idi
  out[11,15] <- pnorm(-abs(y$z.idi))
  out[12,1] <- outcome[i]
  out[12,2] <- "Age≤60"
  out[12,3] <- "age+sex+allriskfactors+rs10490770"
  out[12,4] <- ROC12$auc
  out[12,5] <- ci.auc(ROC12, conf.level=0.95)[1]
  out[12,6] <- ci.auc(ROC12, conf.level=0.95)[3]
  out[12,7] <- roc.test(ROC11, ROC12, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk12))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk11, predrisk2=predRisk12, cutoff)
  out[12,8] <- y$nri
  out[12,9] <- y$nri -  1.96 * y$se.nri
  out[12,10] <- y$nri + 1.96 * y$se.nri
  out[12,11] <- pnorm(-abs(y$z.nri))
  out[12,12] <- y$idi
  out[12,13] <- y$idi -  1.96 * y$se.idi
  out[12,14] <- y$idi + 1.96 * y$se.idi
  out[12,15] <- pnorm(-abs(y$z.idi))
  out %>% write.table("results/AUC.tsv",row.names = F, col.names = F, append=T,sep="\t")
}

out <- data.frame(matrix(0, 12, 15))
colnames(out) <- c("outcome","agerange","model","AUC","AUC.LL","AUC.UL","AUC.p-value","NRI","NRI.LL", "NRI.UL","NRI.p-value","IDI","IDI.LL","IDI.UL","IDI.p-value")
for(i in seq(1,8)){
  age60 <- final
  age60 <- age60 %>% filter(pop == "EUR") %>% filter(age_at_diagnosis > 60)
  age60 <- age60[,c("age_at_diagnosis", "sex", "smoke","BMI","com_cancer","com_chronic_kidney","com_chronic_pulm",
                    "com_transplant","com_heart_failure","com_diabetes", outcome[i], "chr3:45823240:T:C_C")]
  age60 <- age60[complete.cases(age60),]
  LM1 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex"), data=age60, family=binomial(link="logit"))
  LM2 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex + smoke"), data=age60, family=binomial(link="logit"))
  LM3 <- glm(paste0(outcome[i]," ~ age_at_diagnosis + sex + BMI "), data=age60, family=binomial(link="logit"))
  LM4 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_cancer"), data=age60, family=binomial(link="logit"))
  LM5 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_chronic_pulm"), data=age60, family=binomial(link="logit"))
  LM6 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_chronic_kidney"), data=age60, family=binomial(link="logit"))
  LM7 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_transplant"), data=age60, family=binomial(link="logit"))
  LM8 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_heart_failure"), data=age60, family=binomial(link="logit"))
  LM9 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + com_diabetes"), data=age60, family=binomial(link="logit"))
  LM10 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + `chr3:45823240:T:C_C`"), data=age60, family=binomial(link="logit"))
  LM11 <- glm(paste0(outcome[i]," ~  age_at_diagnosis + sex + smoke + BMI + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes"), data=age60, family=binomial(link="logit"))
  LM12 <- glm(paste0(outcome[i]," ~ `chr3:45823240:T:C_C` + age_at_diagnosis + sex + smoke + BMI + com_cancer + com_chronic_kidney + 
                    com_chronic_pulm + com_transplant + com_heart_failure + com_diabetes"), data=age60, family=binomial(link="logit"))
  age60$pred1 <- predict(LM1, age60)
  age60$pred2 <- predict(LM2, age60)
  age60$pred3 <- predict(LM3, age60)
  age60$pred4 <- predict(LM4, age60)
  age60$pred5 <- predict(LM5, age60)
  age60$pred6 <- predict(LM6, age60)
  age60$pred7 <- predict(LM7, age60)
  age60$pred8 <- predict(LM8, age60)
  age60$pred9 <- predict(LM9, age60)
  age60$pred10 <- predict(LM10, age60)
  age60$pred11 <- predict(LM11, age60)
  age60$pred12 <- predict(LM12, age60)
  age60 <- data.frame(age60)
  riskmodel1 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2), cGenPreds=0, cGenPredsCat=0)
  riskmodel2 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,3), cGenPreds=0, cGenPredsCat=0)
  riskmodel3 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2), cGenPreds=0, cGenPredsCat=0)
  riskmodel4 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,5), cGenPreds=0, cGenPredsCat=0)
  riskmodel5 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,6), cGenPreds=0, cGenPredsCat=0)
  riskmodel6 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,7), cGenPreds=0, cGenPredsCat=0)
  riskmodel7 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,8), cGenPreds=0, cGenPredsCat=0)
  riskmodel8 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,9), cGenPreds=0, cGenPredsCat=0)
  riskmodel9 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2,10), cGenPreds=0, cGenPredsCat=0)
  riskmodel10 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1), cNonGenPredsCat=c(2), cGenPreds=12, cGenPredsCat=0)
  riskmodel11 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2,3,5:10), cGenPreds=0, cGenPredsCat=0)
  riskmodel12 <- fitLogRegModel(data=age60, cOutcome=11, cNonGenPreds=c(1,4), cNonGenPredsCat=c(2,3,5:10), cGenPreds=12, cGenPredsCat=0)
  predRisk1 <- predRisk(riskModel=riskmodel1)
  predRisk2 <- predRisk(riskmodel2)
  predRisk3 <- predRisk(riskmodel3)
  predRisk4 <- predRisk(riskmodel4)
  predRisk5 <- predRisk(riskmodel5)
  predRisk6 <- predRisk(riskmodel6)
  predRisk7 <- predRisk(riskmodel7)
  predRisk8 <- predRisk(riskmodel8)
  predRisk9 <- predRisk(riskmodel9)
  predRisk10 <- predRisk(riskmodel10)
  predRisk11 <- predRisk(riskmodel11)
  predRisk12 <- predRisk(riskmodel12)
  value<-as.vector(quantile(predRisk2))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk2, cutoff)
  response <- as.numeric(unlist(age60[,outcome[i]]))
  ROC1 <- roc(response=response, predictor=age60$pred1) 
  ROC2 <- roc(response=response, predictor=age60$pred2) 
  ROC3 <- roc(response=response, predictor=age60$pred3) 
  ROC4 <- roc(response=response, predictor=age60$pred4)
  ROC5 <- roc(response=response, predictor=age60$pred5)
  ROC6 <- roc(response=response, predictor=age60$pred6)
  ROC7 <- roc(response=response, predictor=age60$pred7)
  ROC8 <- roc(response=response, predictor=age60$pred8)
  ROC9 <- roc(response=response, predictor=age60$pred9)
  ROC10 <- roc(response=response, predictor=age60$pred10)
  ROC11 <- roc(response=response, predictor=age60$pred11)
  ROC12 <- roc(response=response, predictor=age60$pred12)
  out[1,1] <- outcome[i]
  out[1,2] <- "Age>60"
  out[1,3] <- "age+sex"
  out[1,4] <- ROC1$auc
  out[1,5] <- ci.auc(ROC1, conf.level=0.95)[1]
  out[1,6] <- ci.auc(ROC1, conf.level=0.95)[3]
  out[2,1] <- outcome[i]
  out[2,2] <- "Age>60"
  out[2,3] <- "age+sex+smoking"
  out[2,4] <- ROC2$auc
  out[2,5] <- ci.auc(ROC2, conf.level=0.95)[1]
  out[2,6] <- ci.auc(ROC2, conf.level=0.95)[3]
  out[2,7] <- roc.test(ROC1, ROC2, method="delong", alternative="two.sided")$p
  out[2,8] <- y$nri
  out[2,9] <- y$nri -  1.96 * y$se.nri
  out[2,10] <- y$nri + 1.96 * y$se.nri
  out[2,11] <- pnorm(-abs(y$z.nri))
  out[2,12] <- y$idi
  out[2,13] <- y$idi -  1.96 * y$se.idi
  out[2,14] <- y$idi + 1.96 * y$se.idi
  out[2,15] <- pnorm(-abs(y$z.idi))
  out[3,1] <- outcome[i]
  out[3,2] <- "Age>60"
  out[3,3] <- "age+sex+BMI"
  out[3,4] <- ROC3$auc
  out[3,5] <- ci.auc(ROC3, conf.level=0.95)[1]
  out[3,6] <- ci.auc(ROC3, conf.level=0.95)[3]
  out[3,7] <- roc.test(ROC1, ROC3, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk3))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk3, cutoff)
  out[3,8] <- y$nri
  out[3,9] <- y$nri -  1.96 * y$se.nri
  out[3,10] <- y$nri + 1.96 * y$se.nri
  out[3,11] <- pnorm(-abs(y$z.nri))
  out[3,12] <- y$idi
  out[3,13] <- y$idi -  1.96 * y$se.idi
  out[3,14] <- y$idi + 1.96 * y$se.idi
  out[3,15] <- pnorm(-abs(y$z.idi))
  out[4,1] <- outcome[i]
  out[4,2] <- "Age>60"
  out[4,3] <- "age+sex+cancer"
  out[4,4] <- ROC4$auc
  out[4,5] <- ci.auc(ROC4, conf.level=0.95)[1]
  out[4,6] <- ci.auc(ROC4, conf.level=0.95)[3]
  out[4,7] <- roc.test(ROC1, ROC4, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk4))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk4, cutoff)
  out[4,8] <- y$nri
  out[4,9] <- y$nri -  1.96 * y$se.nri
  out[4,10] <- y$nri + 1.96 * y$se.nri
  out[4,11] <- pnorm(-abs(y$z.nri))
  out[4,12] <- y$idi
  out[4,13] <- y$idi -  1.96 * y$se.idi
  out[4,14] <- y$idi + 1.96 * y$se.idi
  out[4,15] <- pnorm(-abs(y$z.idi))
  out[5,1] <- outcome[i]
  out[5,2] <- "Age>60"
  out[5,3] <- "age+sex+CKD"
  out[5,4] <- ROC5$auc
  out[5,5] <- ci.auc(ROC5, conf.level=0.95)[1]
  out[5,6] <- ci.auc(ROC5, conf.level=0.95)[3]
  out[5,7] <- roc.test(ROC1, ROC5, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk5))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk5, cutoff)
  out[5,8] <- y$nri
  out[5,9] <- y$nri -  1.96 * y$se.nri
  out[5,10] <- y$nri + 1.96 * y$se.nri
  out[5,11] <- pnorm(-abs(y$z.nri))
  out[5,12] <- y$idi
  out[5,13] <- y$idi -  1.96 * y$se.idi
  out[5,14] <- y$idi + 1.96 * y$se.idi
  out[5,15] <- pnorm(-abs(y$z.idi))
  out[6,1] <- outcome[i]
  out[6,2] <- "Age>60"
  out[6,3] <- "age+sex+COPD"
  out[6,4] <- ROC6$auc
  out[6,5] <- ci.auc(ROC6, conf.level=0.95)[1]
  out[6,6] <- ci.auc(ROC6, conf.level=0.95)[3]
  out[6,7] <- roc.test(ROC1, ROC6, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk6))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk6, cutoff)
  out[6,8] <- y$nri
  out[6,9] <- y$nri -  1.96 * y$se.nri
  out[6,10] <- y$nri + 1.96 * y$se.nri
  out[6,11] <- pnorm(-abs(y$z.nri))
  out[6,12] <- y$idi
  out[6,13] <- y$idi -  1.96 * y$se.idi
  out[6,14] <- y$idi + 1.96 * y$se.idi
  out[6,15] <- pnorm(-abs(y$z.idi))
  out[7,1] <- outcome[i]
  out[7,2] <- "Age>60"
  out[7,3] <- "age+sex+transplantation"
  out[7,4] <- ROC7$auc
  out[7,5] <- ci.auc(ROC7, conf.level=0.95)[1]
  out[7,6] <- ci.auc(ROC7, conf.level=0.95)[3]
  out[7,7] <- roc.test(ROC1, ROC7, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk7))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk7, cutoff)
  out[7,8] <- y$nri
  out[7,9] <- y$nri -  1.96 * y$se.nri
  out[7,10] <- y$nri + 1.96 * y$se.nri
  out[7,11] <- pnorm(-abs(y$z.nri))
  out[7,12] <- y$idi
  out[7,13] <- y$idi -  1.96 * y$se.idi
  out[7,14] <- y$idi + 1.96 * y$se.idi
  out[7,15] <- pnorm(-abs(y$z.idi))
  out[8,1] <- outcome[i]
  out[8,2] <- "Age>60"
  out[8,3] <- "age+sex+CHF"
  out[8,4] <- ROC8$auc
  out[8,5] <- ci.auc(ROC8, conf.level=0.95)[1]
  out[8,6] <- ci.auc(ROC8, conf.level=0.95)[3]
  out[8,7] <- roc.test(ROC1, ROC8, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk8))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk8, cutoff)
  out[8,8] <- y$nri
  out[8,9] <- y$nri -  1.96 * y$se.nri
  out[8,10] <- y$nri + 1.96 * y$se.nri
  out[8,11] <- pnorm(-abs(y$z.nri))
  out[8,12] <- y$idi
  out[8,13] <- y$idi -  1.96 * y$se.idi
  out[8,14] <- y$idi + 1.96 * y$se.idi
  out[8,15] <- pnorm(-abs(y$z.idi))
  out[9,1] <- outcome[i]
  out[9,2] <- "Age>60"
  out[9,3] <- "age+sex+diabetes"
  out[9,4] <- ROC9$auc
  out[9,5] <- ci.auc(ROC9, conf.level=0.95)[1]
  out[9,6] <- ci.auc(ROC9, conf.level=0.95)[3]
  out[9,7] <- roc.test(ROC1, ROC9, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk9))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk9, cutoff)
  out[9,8] <- y$nri
  out[9,9] <- y$nri -  1.96 * y$se.nri
  out[9,10] <- y$nri + 1.96 * y$se.nri
  out[9,11] <- pnorm(-abs(y$z.nri))
  out[9,12] <- y$idi
  out[9,13] <- y$idi -  1.96 * y$se.idi
  out[9,14] <- y$idi + 1.96 * y$se.idi
  out[9,15] <- pnorm(-abs(y$z.idi))
  out[10,1] <- outcome[i]
  out[10,2] <- "Age>60"
  out[10,3] <- "age+sex+rs10490770"
  out[10,4] <- ROC10$auc
  out[10,5] <- ci.auc(ROC10, conf.level=0.95)[1]
  out[10,6] <- ci.auc(ROC10, conf.level=0.95)[3]
  out[10,7] <- roc.test(ROC1, ROC10, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk10))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk10, cutoff)
  out[10,8] <- y$nri
  out[10,9] <- y$nri -  1.96 * y$se.nri
  out[10,10] <- y$nri + 1.96 * y$se.nri
  out[10,11] <- pnorm(-abs(y$z.nri))
  out[10,12] <- y$idi
  out[10,13] <- y$idi -  1.96 * y$se.idi
  out[10,14] <- y$idi + 1.96 * y$se.idi
  out[10,15] <- pnorm(-abs(y$z.idi))
  out[11,1] <- outcome[i]
  out[11,2] <- "Age>60"
  out[11,3] <- "age+sex+allriskfactors"
  out[11,4] <- ROC11$auc
  out[11,5] <- ci.auc(ROC11, conf.level=0.95)[1]
  out[11,6] <- ci.auc(ROC11, conf.level=0.95)[3]
  out[11,7] <- roc.test(ROC1, ROC11, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk11))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk1, predrisk2=predRisk11, cutoff)
  out[11,8] <- y$nri
  out[11,9] <- y$nri -  1.96 * y$se.nri
  out[11,10] <- y$nri + 1.96 * y$se.nri
  out[11,11] <- pnorm(-abs(y$z.nri))
  out[11,12] <- y$idi
  out[11,13] <- y$idi -  1.96 * y$se.idi
  out[11,14] <- y$idi + 1.96 * y$se.idi
  out[11,15] <- pnorm(-abs(y$z.idi))
  out[12,1] <- outcome[i]
  out[12,2] <- "Age>60"
  out[12,3] <- "age+sex+allriskfactors+rs10490770"
  out[12,4] <- ROC12$auc
  out[12,5] <- ci.auc(ROC12, conf.level=0.95)[1]
  out[12,6] <- ci.auc(ROC12, conf.level=0.95)[3]
  out[12,7] <- roc.test(ROC11, ROC12, method="delong", alternative="two.sided")$p
  value<-as.vector(quantile(predRisk12))
  cutoff <- c(0,value[2],value[3],value[4],1)
  y <- reclassification_rev(data=age60, cOutcome=11,
                            predrisk1=predRisk11, predrisk2=predRisk12, cutoff)
  out[12,8] <- y$nri
  out[12,9] <- y$nri -  1.96 * y$se.nri
  out[12,10] <- y$nri + 1.96 * y$se.nri
  out[12,11] <- pnorm(-abs(y$z.nri))
  out[12,12] <- y$idi
  out[12,13] <- y$idi -  1.96 * y$se.idi
  out[12,14] <- y$idi + 1.96 * y$se.idi
  out[12,15] <- pnorm(-abs(y$z.idi))
  out %>% write.table("results/AUC.tsv",row.names = F, col.names = F, append=T,sep="\t")
}

out <- fread("results/AUC.tsv")
colnames(out) <- c("outcome","agerange","model","AUC","AUC.LL","AUC.UL","AUC.p-value","NRI","NRI.LL", "NRI.UL","NRI.p-value","IDI","IDI.LL","IDI.UL","IDI.p-value")
write.table(out, "results/AUC.tsv",row.names = F, col.names = T, append=F,sep="\t")
write.xlsx(out, file="results/AUC.xlsx")

