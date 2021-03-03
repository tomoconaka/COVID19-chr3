setwd("/home/tomoko/02.associations")
final <- readRDS("../01.data.QC/curated_clinical/complication_dat_final.rds")
snps <- colnames(final)[grepl("chr", colnames(final))][3:12]
outcome <- c("resp_severe","resp_mild", "icu_admit_rev" ,"vte", "cardiac",
             "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")

final <- final %>% mutate_at(.vars = c(outcome), .funs = funs(ifelse(. == -1, NA, .)))
final <- final %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))
data_EUR <- final %>% filter(pop == "EUR")

data_EUR <- data_EUR %>% mutate(snp = round(`chr3:45823240:T:C_C`),
                                        snp1 = ifelse(snp >= 1, 1, 0),
                                        snp2 = ifelse(snp == 2, 1, 0))
LM <- glmer(formula = as.formula(paste0("a1 ~ snp + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")),
    family = "binomial", data = data_EUR)

summary(LM)
#AIC      BIC   logLik deviance df.resid 
#5534.6   5605.7  -2757.3   5514.6     8973 

LM <- glmer(formula = as.formula(paste0("a1 ~ snp1 + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")),
            family = "binomial", data = data_EUR)

summary(LM)
#     AIC      BIC   logLik deviance df.resid 
# 5534.3   5605.3  -2757.1   5514.3     8973 

LM <- glmer(formula = as.formula(paste0("a1 ~ snp2 + age_at_diagnosis + sex + (1 | nation) + PC1 + PC2 + PC3 + PC4 + PC5")),
            family = "binomial", data = data_EUR)

summary(LM)
#     AIC      BIC   logLik deviance df.resid 
#5567.2   5638.2  -2773.6   5547.2     8973

