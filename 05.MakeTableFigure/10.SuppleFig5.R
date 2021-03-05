setwd("/home/tomoko/02.associations/")

hgi_EUR <- fread("results/clinical_EUR.results") %>% mutate(pop="COVID19-HGI-EUR")
hgi_SAS <- fread("results/clinical_SAS.results")  %>% mutate(pop="COVID19-HGI-SAS")
hgi_AMR <- fread("results/clinical_AMR.results")  %>% mutate(pop="COVID19-HGI-AMR")
hgi_AFR <- fread("results/clinical_AFR.results")  %>% mutate(pop="COVID19-HGI-AFR")
hgi_EAS <- fread("results/clinical_EAS.results")  %>% mutate(pop="COVID19-HGI-EAS")


fingen <- readxl::read_excel("~/covid19-hgi-clinical-values/HGI/Finngen_20210126_complications.xlsx", sheet = 1)
fingen <- fingen %>% mutate(Outcome = case_when(outcome == "resp_severe" ~ "resp_severe",
                                                outcome == "high_who_score" ~ "a1",
                                                outcome == "hospitalization" ~ "hospitalization",
                                                outcome == "icu_admission" ~ "icu_admit_rev"))
fingen <- fingen %>% filter(age_group == "ALL") %>% mutate(pop = "FinnGen-EUR")
fingen <- fingen %>% mutate(MAC = 30)
cubb <- readxl::read_excel("~/covid19-hgi-clinical-values/HGI/Chr_3_analysis_CUBB.v2.xlsx", sheet = 1)
cubb <- cubb %>% mutate(Outcome = case_when(outcome == "resp_severe1" ~ "resp_severe",
                                            outcome == "resp_severe2" ~ "a1",
                                            outcome == "hospitalization" ~ "hospitalization",
                                            outcome == "icu_admission" ~ "icu_admit"))
cubb <- cubb %>% mutate(pop = case_when(pop == "EUR" ~ "CUBB-EUR",
                                        pop == "AFR" ~ "CUBB-AFR",
                                        pop == "Hispanic" ~ "CUBB-AMR"))

cubb <- cubb %>% mutate(MAC = (AF_Cases*N_case + AF_Controls*N_control))

all <- bind_rows(hgi_EUR,hgi_SAS,hgi_AMR,hgi_AFR,hgi_EAS,fingen)
all <- all %>% mutate(beta = ifelse(MAC < 10 | N_case < 10 | N_control < 10, NA, beta),
                      se = ifelse(MAC < 10 | N_case < 10 | N_control < 10, NA, se),
                      pvalue = ifelse(MAC < 10 | N_case < 10 | N_control < 10, NA, pvalue))

tmp <- all %>% filter(Outcome == "hospitalization")
tmp <- tmp[c(1,6,2,3,4,5),]
colnames(tmp)[c(5:6,13)] <- c("N cases", "N controls","Study")
m1 <- metagen(beta,
              se,
              data=tmp,
              studlab=paste(Study),
              comb.fixed = TRUE,
              comb.random = TRUE,
              prediction=FALSE,
              sm="OR")
png("plots/SuppleFig5_1.png", width=700, height=300)
forest(m1, smlab="",xlim=c(0.5,8),leftcols=c("studlab", "N cases", "N controls"),
       rightcols=c("effect", "ci"),print.I2.ci = TRUE,print.tau2 = FALSE,
       leftlabs = c("Ancestry", "N cases", "N controls"),zero.pval = TRUE)
dev.off()

tmp <- all %>% filter(Outcome == "icu_admit_rev")
tmp <- tmp[c(1,6,2,3,4,5),]
colnames(tmp)[c(5:6,13)] <- c("N cases", "N controls","Study")
m1 <- metagen(beta,
              se,
              data=tmp,
              studlab=paste(Study),
              comb.fixed = TRUE,
              comb.random = TRUE,
              prediction=FALSE,
              sm="OR")
png("plots/SuppleFig5_2.png", width=700, height=300)
forest(m1, smlab="",xlim=c(0.5,8),leftcols=c("studlab", "N cases", "N controls"),
       rightcols=c("effect", "ci"),print.I2.ci = TRUE,print.tau2 = FALSE,
       leftlabs = c("Ancestry", "N cases", "N controls"),zero.pval = TRUE)
dev.off()

tmp <- all %>% filter(Outcome == "a1")
tmp <- tmp[c(1,6,2,3,4,5),]
colnames(tmp)[c(5:6,13)] <- c("N cases", "N controls","Study")
m1 <- metagen(beta,
              se,
              data=tmp,
              studlab=paste(Study),
              comb.fixed = TRUE,
              comb.random = TRUE,
              prediction=FALSE,
              sm="OR")
png("plots/SuppleFig5_3.png", width=700, height=300)
forest(m1, smlab="",xlim=c(0.5,8),leftcols=c("studlab", "N cases", "N controls"),
       rightcols=c("effect", "ci"),print.I2.ci = TRUE,print.tau2 = FALSE,
       leftlabs = c("Ancestry", "N cases", "N controls"),zero.pval = TRUE)
dev.off()
