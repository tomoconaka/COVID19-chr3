setwd("/home/tomoko/02.associations/")

hgi_EUR <- fread("results/clinical_EUR.results") %>% mutate(pop="COVID19-HGI-EUR")
hgi_SAS <- fread("results/clinical_SAS.results")  %>% mutate(pop="COVID19-HGI-SAS")
hgi_AMR <- fread("results/clinical_AMR.results")  %>% mutate(pop="COVID19-HGI-AMR")
hgi_AFR <- fread("results/clinical_AFR.results")  %>% mutate(pop="COVID19-HGI-AFR")
hgi_EAS <- fread("results/clinical_EAS.results")  %>% mutate(pop="COVID19-HGI-EAS")

hgi <- bind_rows(hgi_EUR,hgi_SAS,hgi_AMR,hgi_AFR,hgi_EAS) %>%
  mutate(beta = beta_age_snp, se=se_age_snp, pvalue = pvalue_age_snp)

fingen <- read.table("~/covid19-hgi-clinical-values/HGI/Finngen_20210219_age.csv", sep=";", header=T)
fingen <- fingen %>% mutate(beta = as.numeric(gsub(",",".",beta)),
                            se = as.numeric(gsub(",",".",se)),
                            pvalue = as.numeric(gsub(",",".",pvalue)),
                            MAF = as.numeric(gsub(",",".",MAF))
)
fingen <- fingen %>% mutate(Outcome = case_when(group == "resp_severe" ~ "resp_severe",
                                                group == "high_who_score" ~ "a1",
                                                group == "hospitalized" ~ "hospitalization",
                                                group == "icu_admitted" ~ "icu_admit_rev"))
fingen <- fingen %>% mutate(pop = "FinnGen-EUR", N_case=N)
fingen <- fingen %>% mutate(caseMAC = MAF*N)
cubb <- readxl::read_excel("~/covid19-hgi-clinical-values/HGI/Chr_3_analysis_CUBB.v3.xlsx", sheet = 2)
cubb <- cubb %>% mutate(Outcome = case_when(outcome == "resp_severe1" ~ "resp_severe",
                                            outcome == "resp_severe2" ~ "a1",
                                            outcome == "hospitalization" ~ "hospitalization",
                                            outcome == "icu_admission" ~ "icu_admit"))
cubb <- cubb %>% mutate(pop = case_when(pop == "EUR" ~ "CUB-EUR",
                                        pop == "AFR" ~ "CUB-AFR",
                                        pop == "Hispanic" ~ "CUB-AMR"))

cubb <- cubb %>% mutate(caseMAC = (AF_Cases*N_case))

all <- bind_rows(hgi,fingen,cubb)
all$caseMAC[is.na(all$caseMAC)] <- all$MAC[is.na(all$caseMAC)]
all <- all %>% mutate(beta = ifelse(caseMAC < 3 | N_case < 10, NA, beta),
                      se = ifelse(caseMAC < 3 | N_case < 10, NA, se),
                      pvalue = ifelse(caseMAC < 3 | N_case < 10, NA, pvalue))

tmp <- all %>% filter(Outcome == "hospitalization")
tmp <- tmp[c(1,6,2,3,4,5),]
colnames(tmp)[c(5,13)] <- c("N cases", "Study")
m1 <- metagen(beta,
              se,
              data=tmp,
              studlab=paste(Study),
              comb.fixed = TRUE,
              comb.random = TRUE,
              prediction=FALSE,
              sm="MD")
png("plots/SuppleFig7_1.png", width=700, height=300)
forest(m1, smlab="",xlim=c(-5,5),leftcols=c("studlab", "N cases"),
       rightcols=c("effect", "ci"),print.I2.ci = TRUE,print.tau2 = FALSE,
       leftlabs = c("Ancestry", "N cases"),rightlabs = c("Coefficients","95%CI"),zero.pval = TRUE)
dev.off()

tmp <- all %>% filter(Outcome == "icu_admit_rev")
tmp <- tmp[c(1,6,2,3,4,5),]
colnames(tmp)[c(5,13)] <- c("N cases", "Study")
m1 <- metagen(beta,
              se,
              data=tmp,
              studlab=paste(Study),
              comb.fixed = TRUE,
              comb.random = TRUE,
              prediction=FALSE,
              sm="MD")
png("plots/SuppleFig7_2.png", width=700, height=300)
forest(m1, smlab="",xlim=c(-5,5),leftcols=c("studlab", "N cases"),
       rightcols=c("effect", "ci"),print.I2.ci = TRUE,print.tau2 = FALSE,
       leftlabs = c("Ancestry", "N cases"),rightlabs = c("Coefficients","95%CI"),zero.pval = TRUE)
dev.off()

tmp <- all %>% filter(Outcome == "a1")
tmp <- tmp[c(1,6,7,3,9,4,8,2,5),]
colnames(tmp)[c(5,13)] <- c("N cases", "Study")
m1 <- metagen(beta,
              se,
              data=tmp,
              studlab=paste(Study),
              comb.fixed = TRUE,
              comb.random = TRUE,
              prediction=FALSE,
              sm="MD")
png("plots/SuppleFig7_3.png", width=700, height=300)
forest(m1, smlab="",xlim=c(-5,5),leftcols=c("studlab", "N cases"),
       rightcols=c("effect", "ci"),print.I2.ci = TRUE,print.tau2 = FALSE,
       leftlabs = c("Ancestry", "N cases"),rightlabs = c("Coefficients","95%CI"),zero.pval = TRUE)
dev.off()
