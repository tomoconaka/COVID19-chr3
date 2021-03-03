setwd("/home/tomoko/05.survival/")
demog_data <-  readRDS("../01.data.QC/curated_clinical/demog_data.rds")
hosp_data <- readRDS("../01.data.QC/curated_clinical/hosp_data.rds")
covid_data <- readRDS("../01.data.QC/curated_clinical/covid_data.rds")
comorb_data <- readRDS("../01.data.QC/curated_clinical/comorb_data.rds") %>% select(-"study") 
geno <- readRDS("../covid19-hgi-clinical-values/bgen/all_ans_pca_raw.rds")
data1 <- hosp_data %>% inner_join(demog_data[,-7], by ="anonymized_patient_id") %>% 
  inner_join(covid_data[,-6], by ="anonymized_patient_id") %>% inner_join(geno, by ="anonymized_patient_id") %>% inner_join(comorb_data, by ="anonymized_patient_id") %>% 
  mutate_at(.vars=c(colnames(comorb_data)[grepl("com_",colnames(comorb_data))]), .funs=funs(as.numeric(.)))
ukb <- readRDS("../01.data.QC/ukbb/clinical_value_20210217.rds")
ukb$date_of_death[is.na(ukb$date_of_death)] <- max(as.Date(ukb$date_of_death), na.rm=T)
genoukb <- readRDS("../01.data.QC/ukbb/all_chr.rds")
pcukb <- fread("~/11.pca/cases_controls_pca_sup_pops_0.5_probs.txt")
genoukb <- genoukb %>% merge(pcukb, by.x="FID",by.y="s")

ukb <- ukb %>% inner_join(genoukb, by=c("anonymized_patient_id"="FID")) %>% mutate(anonymized_patient_id = as.character(anonymized_patient_id))
ukb <- ukb  %>% mutate(anonymized_patient_id = paste0("UKB_",anonymized_patient_id))
ukb <- ukb %>% mutate(age2 = age_at_diagnosis^2)
ukb <- ukb %>% select(any_of(colnames(data1))) %>% mutate(anonymized_patient_id = as.character(anonymized_patient_id))
ukb <- ukb %>% mutate(study = "UKB")
final <- bind_rows(ukb, data1)

final <- final %>% mutate(date_of_death = case_when(!is.na(date_of_death) ~ date_of_death,
                                                    study == "belgium_rahmouni" & hospitalization == 0 ~ as.Date(covid19_test_date + 30),
                                                    TRUE ~ hospitalization_end),
                          covid19_test_date = case_when(!is.na(covid19_test_date) ~ covid19_test_date,
                                                        hospitalization_start > as.Date("2020-03-01") ~ hospitalization_start,
                                                        !is.na(covid19_first_symptoms_date) ~ covid19_first_symptoms_date))

final <- final %>% mutate(time = as.numeric(date_of_death - covid19_test_date))

final <- final %>% mutate(smoke = case_when(smoking == 0 ~ "current",
                                            smoking == 1 ~ "ex",
                                            smoking == 2 ~ "0none",
                                            TRUE ~ "NA"),
                          BMI = case_when(weight/((height/100)^2) >= 40 ~ "class4",
                                          weight/((height/100)^2) < 40 & weight/((height/100)^2) >=35 ~ "class3",
                                          weight/((height/100)^2) < 35 & weight/((height/100)^2) >=30 ~ "class2",
                                          weight/((height/100)^2) < 30 ~ "0none",
                                          TRUE ~ "NA"))


final <- final %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))
final$study[is.na(final$study)] <- "UKB"

final <- final %>% mutate_at(.vars = c("death"), .funs = funs(ifelse(. == -1, NA, .)))
final <- final %>% mutate_at(.vars = c("smoke","BMI"), .funs = funs(ifelse(. == "NA", NA, .)))

final <- final %>% mutate(time = ifelse(time > -30 & time < 0, 0, time)) %>% filter(time >= 0) %>%
  filter(covid19_test_date >= as.Date("2020-02-05"))# %>% mutate(death = ifelse(time > 30, 0, death)) 
final <- final %>% drop_na(death)
final <- final %>% drop_na(date_of_death) %>% drop_na(covid19_test_date) %>% drop_na(cause_of_death) %>% drop_na(PC1)

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
final <- final %>% mutate(hospitalization = ifelse(death == 1, 1, hospitalization))
library("survival")
library(ggplot2)
library(dplyr)
library(ggfortify)
library(frailtyHL)
library(coxme)
final <- final %>% mutate(age2 = age_at_diagnosis*age_at_diagnosis)
data_EUR <- final %>% filter(pop=="EUR")
quantile(data_EUR$time)
fit <- survfit(Surv(time, death) ~ snp, data = data_EUR)
summary(fit)
summary(survfit(Surv(time, death) ~ snp, data = data_EUR), times = 30)

cox <- coxme(Surv(time, death) ~ snp + sex + age_at_diagnosis + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5, data = data_EUR)

save(data_EUR,fit,cox,file="all_cause_EUR.obj")

coxmeTable <- function (mod){
  if(!any(class(mod)=="coxme")){stop("Model not from mixed effects Cox model")}
  beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$variance) - nvar
  se <- sqrt(diag(as.matrix(mod$variance))[nfrail + 1:nvar])
  z <- beta/se
  p <- 2*(1-pnorm(abs(z)))
  #p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

out <- data.frame(matrix(0,10,9))
colnames(out) <- c("beta", "se", "pvalue", "pop", "type", "tpos", "fpos", "tneg", "fneg")
out[1,1] <- coxmeTable(cox)[1,1]
out[1,2] <- coxmeTable(cox)[1,2]
out[1,3] <- coxmeTable(cox)[1,4]
out[1,4] <- "EUR"
out[1,5] <- "all_cause"
out[1,6] <- sum(data_EUR$death == 1 & data_EUR$snp == 1,na.rm=T)
out[1,7] <- sum(data_EUR$death == 0 & data_EUR$snp == 1,na.rm=T)
out[1,8] <- sum(data_EUR$death == 0 & data_EUR$snp == 0,na.rm=T)
out[1,9] <- sum(data_EUR$death == 1 & data_EUR$snp == 0,na.rm=T)

data_SAS <- final %>% filter(pop=="SAS")
sum(data_SAS$death)
sum(data_SAS$snp)
cox <- coxme(Surv(time, death) ~ snp + sex + age_at_diagnosis + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5, data = data_SAS)
out[2,1] <- coxmeTable(cox)[1,1]
out[2,2] <- coxmeTable(cox)[1,2]
out[2,3] <- coxmeTable(cox)[1,4]
out[2,4] <- "SAS"
out[2,5] <- "all_cause"
out[2,6] <- sum(data_SAS$death == 1 & data_SAS$snp == 1)
out[2,7] <- sum(data_SAS$death == 0 & data_SAS$snp == 1)
out[2,8] <- sum(data_SAS$death == 0 & data_SAS$snp == 0)
out[2,9] <- sum(data_SAS$death == 1 & data_SAS$snp == 0)

data_AMR <- final %>% filter(pop=="AMR")
sum(data_AMR$death)
sum(data_AMR$snp)
cox <- coxme(Surv(time, death) ~ snp + sex + age_at_diagnosis + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5, data = data_AMR)
out[3,1] <- coxmeTable(cox)[1,1]
out[3,2] <- coxmeTable(cox)[1,2]
out[3,3] <- coxmeTable(cox)[1,4]
out[3,4] <- "AMR"
out[3,5] <- "all_cause"
out[3,6] <- sum(data_AMR$death == 1 & data_AMR$snp == 1,na.rm=T)
out[3,7] <- sum(data_AMR$death == 0 & data_AMR$snp == 1,na.rm=T)
out[3,8] <- sum(data_AMR$death == 0 & data_AMR$snp == 0,na.rm=T)
out[3,9] <- sum(data_AMR$death == 1 & data_AMR$snp == 0,na.rm=T)

data_AFR <- final %>% filter(pop=="AFR")
sum(data_AFR$death)
sum(data_AFR$snp)#2
cox <- coxme(Surv(time, death) ~ snp + sex + age_at_diagnosis + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5, data = data_AFR)
out[4,1] <- NA
out[4,2] <- NA
out[4,3] <- NA
out[4,4] <- "AFR"
out[4,5] <- "all_cause"
out[4,6] <- sum(data_AFR$death == 1 & data_AFR$snp == 1)
out[4,7] <- sum(data_AFR$death == 0 & data_AFR$snp == 1)
out[4,8] <- sum(data_AFR$death == 0 & data_AFR$snp == 0)
out[4,9] <- sum(data_AFR$death == 1 & data_AFR$snp == 0)

data_EAS <- final %>% filter(pop=="EAS")
sum(data_EAS$death)#2
sum(data_EAS$snp)
out[5,1] <- NA
out[5,2] <- NA
out[5,3] <- NA
out[5,4] <- "EAS"
out[5,5] <- "all_cause"
out[5,6] <- sum(data_EAS$death == 1 & data_EAS$snp == 1)
out[5,7] <- sum(data_EAS$death == 0 & data_EAS$snp == 1)
out[5,8] <- sum(data_EAS$death == 0 & data_EAS$snp == 0)
out[5,9] <- sum(data_EAS$death == 1 & data_EAS$snp == 0)

#compatitive risk
library(cmprsk)
final <- final %>% mutate(death = ifelse(cause_of_death == 1 & death == 1, 2, death))
data_EUR <- final %>% filter(pop == "EUR") %>% drop_na(PC1)
ci_fit <- 
  cuminc(
    ftime = data_EUR$time, 
    fstatus = data_EUR$death,
    group = data_EUR$snp,
    cencode = 0
  )

fit5 <- survfit(Surv(time, death, type = "mstate") ~ snp, data = data_EUR)

covs1 <- model.matrix(~ snp + age_at_diagnosis + sex + nation + PC1 + PC2 + PC3 + PC4 + PC5, data = data_EUR)[, -1]

shr_fit <- 
  crr(
    ftime = data_EUR$time,
    fstatus = data_EUR$death,
    cov1 = covs1,
    cencode = 0, failcode = 2
  )

res <- summary.crr(shr_fit, conf.int = 0.95)

save(res, shr_fit, ci_fit, fit5,covs1,data_EUR,file="covid_cause_EUR.obj")
png("plots/covid_cause_mortality_EUR.png", width=500, height=500)
plot(ci_fit[3:4], xlim=c(0,30), ylim=c(0,0.2), xlab="time (days)", ylab="COVID-19 related mortality rate", 
     main="COVID-19 related mortality (European ancestry)", lty=1:1, color=c("pink", "lightblue"), lw=5)
legend(0.5,0.2, title = "rs10490770", 
       legend=c(paste0("TT (N =",dim(data_EUR[data_EUR$snp==0,])[1],")"), paste0("TC or CC (N =",dim(data_EUR[data_EUR$snp==1,])[1],")")),
       col=c("pink", "lightblue"), lty=1:1, cex=1.5, lw=5)
text(10,0.13, cex=1.5, paste0("HR: ",round(c(res$conf.int[1,c(1,3,4)], res$coef[1,5]) ,2)[1]," [",round(c(res$conf.int[1,c(1,3,4)], res$coef[1,5]) ,2)[2],", ",round(c(res$conf.int[1,c(1,3,4)], res$coef[1,5]) ,2)[3],"]\np = ",round(c(res$conf.int[1,c(1,3,4)], res$coef[1,5]) ,5)[4]))
dev.off()
out[6,1] <- res$coef[1,1]
out[6,2] <- res$coef[1,3]
out[6,3] <- res$coef[1,5]
out[6,4] <- "EUR"
out[6,5] <- "covid_related"
out[6,6] <- sum(data_EUR$death == 2 & data_EUR$snp == 1, na.rm=T)
out[6,7] <- sum(data_EUR$death == 0 & data_EUR$snp == 1, na.rm=T)
out[6,8] <- sum(data_EUR$death == 0 & data_EUR$snp == 0, na.rm=T)
out[6,9] <- sum(data_EUR$death == 2 & data_EUR$snp == 0, na.rm=T)

data_SAS <- final %>% filter(pop=="SAS")
sum(data_SAS$death)
sum(data_SAS$snp == 1)

ci_fit <- 
  cuminc(
    ftime = data_SAS$time, 
    fstatus = data_SAS$death,
    group = data_SAS$snp,
    cencode = 0
  )

covs1 <- model.matrix(~ snp + age_at_diagnosis + sex + nation + PC1 + PC2 + PC3 + PC4 + PC5, data = data_SAS)[, -1]

shr_fit <- 
  crr(
    ftime = data_SAS$time,
    fstatus = data_SAS$death,
    cov1 = covs1,
    cencode = 0, failcode = 2
  )

res <- summary.crr(shr_fit, conf.int = 0.95)
out[7,1] <- res$coef[1,1]
out[7,2] <- res$coef[1,3]
out[7,3] <- res$coef[1,5]
out[7,4] <- "SAS"
out[7,5] <- "covid_related"
out[7,6] <- sum(data_SAS$death == 2 & data_SAS$snp == 1, na.rm=T)
out[7,7] <- sum(data_SAS$death == 0 & data_SAS$snp == 1, na.rm=T)
out[7,8] <- sum(data_SAS$death == 0 & data_SAS$snp == 0, na.rm=T)
out[7,9] <- sum(data_SAS$death == 2 & data_SAS$snp == 0, na.rm=T)

data_AMR <- final %>% filter(pop=="AMR")
sum(data_AMR$death == 2)
sum(data_AMR$snp == 1)
ci_fit <- 
  cuminc(
    ftime = data_AMR$time, 
    fstatus = data_AMR$death,
    group = data_AMR$snp,
    cencode = 0
  )
covs1 <- model.matrix(~ snp + age_at_diagnosis + sex + nation + PC1 + PC2 + PC3 + PC4 + PC5, data = data_AMR)[, -1]

shr_fit <- 
  crr(
    ftime = data_AMR$time,
    fstatus = data_AMR$death,
    cov1 = covs1,
    cencode = 0, failcode = 2
  )

res <- summary.crr(shr_fit, conf.int = 0.95)
out[8,1] <- res$coef[1,1]
out[8,2] <- res$coef[1,3]
out[8,3] <- res$coef[1,5]
out[8,4] <- "AMR"
out[8,5] <- "covid_related"
out[8,6] <- sum(data_AMR$death == 2 & data_AMR$snp == 1, na.rm=T)
out[8,7] <- sum(data_AMR$death == 0 & data_AMR$snp == 1, na.rm=T)
out[8,8] <- sum(data_AMR$death == 0 & data_AMR$snp == 0, na.rm=T)
out[8,9] <- sum(data_AMR$death == 2 & data_AMR$snp == 0, na.rm=T)


data_AFR <- final %>% filter(pop=="AFR")
sum(data_AFR$death == 2)
sum(data_AFR$snp == 1)
ci_fit <- 
  cuminc(
    ftime = data_AFR$time, 
    fstatus = data_AFR$death,
    group = data_AFR$snp,
    cencode = 0
  )
covs1 <- model.matrix(~ snp + age_at_diagnosis + sex, data = data_AFR)[, -1]

shr_fit <- 
  crr(
    ftime = data_AFR$time,
    fstatus = data_AFR$death,
    cov1 = covs1,
    cencode = 0, failcode = 2
  )

res <- summary.crr(shr_fit, conf.int = 0.95)
out[9,1] <- NA
out[9,2] <- NA
out[9,3] <- NA
out[9,4] <- "AFR"
out[9,5] <- "covid_related"
out[9,6] <- sum(data_AFR$death == 2 & data_AFR$snp == 1, na.rm=T)
out[9,7] <- sum(data_AFR$death == 0 & data_AFR$snp == 1, na.rm=T)
out[9,8] <- sum(data_AFR$death == 0 & data_AFR$snp == 0, na.rm=T)
out[9,9] <- sum(data_AFR$death == 2 & data_AFR$snp == 0, na.rm=T)

data_EAS <- final %>% filter(pop=="EAS")
out[10,1] <- NA
out[10,2] <- NA
out[10,3] <- NA
out[10,4] <- "EAS"
out[10,5] <- "covid_related"
out[10,6] <- sum(data_EAS$death == 2 & data_EAS$snp == 1, na.rm=T)
out[10,7] <- sum(data_EAS$death == 0 & data_EAS$snp == 1, na.rm=T)
out[10,8] <- sum(data_EAS$death == 0 & data_EAS$snp == 0, na.rm=T)
out[10,9] <- sum(data_EAS$death == 2 & data_EAS$snp == 0, na.rm=T)

write.table(out, file="survival.meta.tsv", sep="\t", quote=F, col.names = T,row.names = F)
