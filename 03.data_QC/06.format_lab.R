setwd("/home/tomoko/01.data.QC")
demog_data <-  readRDS("../01.data.QC/curated_clinical/demog_data.rds")
hosp_data <- readRDS("../01.data.QC/curated_clinical/hosp_data.rds")
comorb_data <- readRDS("../01.data.QC/curated_clinical/comorb_data.rds")
covid_data <- readRDS("../01.data.QC/curated_clinical/covid_data.rds")

path <- "../covid19-hgi-clinical-values/"

identification <- c("anonymized_patient_id")

demographics <- c("age_at_diagnosis","sex")

hospital_info <- 
  c("highest_respiratory_support",
    "highest_who_score")

covid19 <-
  c("covid19_test_date",
    "covid19_first_symptoms_date",
    "covid19_test")

blood_values <- c("lab_result_date",
                  "lab_wbc",
                  "lab_lymphocytes",
                  "lab_cd4",
                  "lab_cd8",
                  "lab_neutrophils",
                  "lab_monocytes",
                  "lab_platelets",
                  "lab_eosinophils",
                  "lab_basophils",
                  "lab_crp",
                  "lab_trop_t",
                  "lab_trop_i",
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
                  "lab_aptt",
                  "lab_inr",
                  "lab_creatinine",
                  "lab_nk", "lab_n_l_ratio","lab_na","lab_k", "lab_mcv","lab_lipase", "lab_amylase")

parse_date_covid19 <- function(date,order_type)
{
  date_out <- as_date(rep(NA,length(date)))
  for (i in 1:length(date))
  {
    i <- as.character(date[i])
    if(i=="-1" | is.na(i) | i=="1899-12-29" | i=="1899-12-30") {
      o <- NA}
    else if (nchar(i)==5) {
      o <- as_date(as.numeric(i),origin="1899-12-30") 
    }
    else if (nchar(i)>5) {
      o <- as_date(parse_date_time(i,order_type))
    }
    else {o <- NA}
    date_out[which(i==date)] <- o
  }
  return(date_out)
}

a <- readRDS("../01.data.QC/hgi_belgium/belgium_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("BM_",anonymized_patient_id)) %>% mutate(study ="BelCovid_1")
b <- readRDS("../01.data.QC/hgi_sweden/swe_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("SW_",anonymized_patient_id)) %>% mutate(study ="Sweden")
c <- readRDS("../01.data.QC/hgi_spain_bujanda/spain_bujanda_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("SBuj_",anonymized_patient_id)) %>% mutate(study ="Hostage_1")
d <- readRDS("../01.data.QC/hgi_canada/canada_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("CA_",anonymized_patient_id)) %>% mutate(study ="BQC19")
e <- readRDS("../01.data.QC/hgi_spain_butti/spain_butti_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("SBut_",anonymized_patient_id)) %>% mutate(study ="Hostage_2")
f <- readRDS("../01.data.QC/hgi_germany_schulte/germany_schulte_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("GS_",anonymized_patient_id)) %>% mutate(study ="COMRI")
g <- readRDS("../01.data.QC/hgi_brasil/brazil_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("BZ_",anonymized_patient_id)) %>% mutate(study ="BRACOVID")
i <- readRDS("../01.data.QC/hgi_italy_valenti/italy_valenti_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("ITV_",anonymized_patient_id)) %>% mutate(study ="FoGS")
j <- readRDS("../01.data.QC/hgi_germany_ludwig/germany_ludwig_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("GL_",anonymized_patient_id)) %>% mutate(study ="BosCO")
l <- readRDS("../01.data.QC/hgi_spain_gomez/spain_gomez_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("SG_",anonymized_patient_id)) %>% mutate(study ="Hostage_3")
m <- readRDS("../01.data.QC/hgi_spain_planas/spain_planas_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("SP_",anonymized_patient_id)) %>% mutate(study ="INMUNGEN-CoV2")
n <- readRDS("../01.data.QC/hgi_belgium_rahmouni/belgium_rahmouni_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("BR_",anonymized_patient_id)) %>% mutate(lab_trop_t = -1) %>% mutate(study ="BelCovid_2")
o <- readRDS("../01.data.QC/hgi_norway/norway_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("NO_",anonymized_patient_id)) %>% mutate(study ="NorCoV2")
p <- readRDS("../01.data.QC/hgi_italy_duga/italy_duga_lab.rds") %>% select(any_of(c(identification,blood_values))) %>% mutate(anonymized_patient_id = paste0("ITD_",anonymized_patient_id)) %>% mutate(study ="Hostage_4")


## READ IN ITALIAN DATA
h <- fread("../01.data.QC/hgi_italy/Italy.lab_20201209TN.tsv") 

h <- h %>% dplyr::select(any_of(c(identification, blood_values))) %>% mutate_at(vars(c("anonymized_patient_id")),as.character) %>%
  mutate_at(.vars = vars(any_of(blood_values[-1])), .funs = funs(ifelse(. < 0,NA,.)))

h <- h %>% filter(rowSums(is.na(.)) != 15) 

tmp <- hosp_data %>% filter(study == "italy_renieri") %>% select("anonymized_patient_id","hospitalization_start")

h <- h[,-2] %>% merge(tmp, by="anonymized_patient_id", all.x=T) %>% rename(lab_result_date = hospitalization_start) %>% mutate(anonymized_patient_id = paste0("ITR_",anonymized_patient_id)) %>% mutate(study ="GEN-COVID")

lab_data <- bind_rows(list(a,b,c,d,e,f,g,h,i,j,l,m,n,o,p))

lab_data <- lab_data %>% mutate(lab_result_date = as.Date(lab_result_date)) 
lab_data$lab_result_date[lab_data$lab_result_date < as.Date("2010-01-01")] <- NA
geno <- readRDS("~/covid19-hgi-clinical-values/bgen/all_ans_pca_raw.rds")
lab_data <- lab_data %>% filter(anonymized_patient_id %in% covid_data$anonymized_patient_id) %>%
  filter(anonymized_patient_id %in% geno$anonymized_patient_id)

saveRDS(lab_data, file="curated_clinical/all_lab_data.rds")

lab_data <- readRDS("curated_clinical/all_lab_data.rds")
lab_data <-lab_data[order(lab_data$study),]

outlier <- function (x,method="median",addthres=TRUE){
  
  ID=obs=NULL
med <- median(x)
MAD <-median(abs(med-x))
dtf <- data.frame(ID=seq.int(length(x)), obs=x, outlier=abs(x-med)>5*(MAD/0.6745))
midp <- med
lower <- med-5*(MAD/0.6745)
upper <- med+5*(MAD/0.6745)
outliern <- length(which(dtf=="TRUE"))
if (addthres==TRUE) {
  p <- ggplot(dtf, aes(x=ID, y=obs, label=ID)) + geom_point(aes(colour=outlier)) + geom_text_repel(data = subset(dtf, outlier=="TRUE"), aes(label = ID), size = 2.7, colour="black", box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) + labs(x=paste0("observation ID number\n number of outliers detected=", outliern, "\n(outlier detection method=", method, ")"), y="observation value") + geom_hline(yintercept = midp, colour="black", linetype = "longdash") + geom_hline(yintercept = lower, colour="black", linetype = "longdash") + geom_hline(yintercept = upper, colour="black", linetype = "longdash")
} else {
  p <- ggplot(dtf, aes(x=ID, y=obs, label=ID)) + geom_point(aes(colour=outlier)) + geom_text_repel(data = subset(dtf, outlier=="TRUE"), aes(label = ID), size = 2.7, colour="black", box.padding = unit(0.35, "lines"), point.padding = unit(0.3, "lines")) + labs(x=paste0("observation ID number\n number of outliers detected=", outliern, "\n(outlier detection method=", method, ")"), y="observation value")
}
print(p)
results <- (list(method=method, midpoint=midp, lowerbound=lower, upperbound=upper, outlN=outliern, flaggedData=dtf))
return(results)
}

lab_data$index <- seq(1,dim(lab_data)[1])

png("plots/labs/lab_wbc.png", height=700)
lab_data <- lab_data %>% mutate(lab_wbc = case_when(lab_wbc == 0 & !is.na(lab_wbc) ~ min(lab_data$lab_wbc[lab_data$lab_wbc != 0], na.rm=T)/2,
                                          TRUE ~ lab_wbc))
tmp <- lab_data %>% drop_na(lab_wbc) %>% filter(lab_wbc > 0 & !is.infinite(lab_wbc))
tmp$lab_wbc <- log(tmp$lab_wbc)
LM <- glm(lab_wbc ~ study, data=tmp, family="gaussian")
tmp$lab_wbc <- residuals(LM)
tmp1 <- outlier(tmp$lab_wbc,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_wbc, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_wbc (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c((dim(tmp)[2]-1),dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_wbc = ifelse(outlier==TRUE, -1, lab_wbc))
tmp <- lab_data %>% drop_na(lab_wbc) %>% filter(lab_wbc > 0 & !is.infinite(lab_wbc))
p2 <- ggplot(tmp, aes(x=study, y=lab_wbc, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_wbc")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_lymphocytes.png", height=700)
lab_data <- lab_data %>% mutate(lab_lymphocytes = case_when(lab_lymphocytes == 0 & !is.na(lab_lymphocytes) ~ min(lab_data$lab_lymphocytes[lab_data$lab_lymphocytes != 0], na.rm=T)/2,
                                                    TRUE ~ lab_lymphocytes))
tmp <- lab_data %>% drop_na(lab_lymphocytes) %>% filter(lab_lymphocytes > 0 & !is.infinite(lab_lymphocytes))
tmp$lab_lymphocytes <- log(tmp$lab_lymphocytes)
LM <- glm(lab_lymphocytes ~ study, data=tmp, family="gaussian")
tmp$lab_lymphocytes <- residuals(LM)
tmp1 <- outlier(tmp$lab_lymphocytes,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_lymphocytes, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_lymphocytes (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_lymphocytes = ifelse(outlier==TRUE, -1, lab_lymphocytes))
tmp <- lab_data %>% drop_na(lab_lymphocytes) %>% filter(lab_lymphocytes > 0 & !is.infinite(lab_lymphocytes))
p2 <- ggplot(tmp, aes(x=study, y=lab_lymphocytes, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_lymphocytes")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_neutrophils.png", height=700)
lab_data <- lab_data %>% mutate(lab_neutrophils = case_when(lab_neutrophils == 0 & !is.na(lab_neutrophils) ~ min(lab_data$lab_neutrophils[lab_data$lab_neutrophils != 0], na.rm=T)/2,
                                                            TRUE ~ lab_neutrophils))
tmp <- lab_data %>% drop_na(lab_neutrophils) %>% filter(lab_neutrophils > 0 & !is.infinite(lab_neutrophils))
tmp$lab_neutrophils <- log(tmp$lab_neutrophils)
LM <- glm(lab_neutrophils ~ study, data=tmp, family="gaussian")
tmp$lab_neutrophils <- residuals(LM)
tmp1 <- outlier(tmp$lab_neutrophils,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_neutrophils, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_neutrophils (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_neutrophils = ifelse(outlier==TRUE, -1, lab_neutrophils))
tmp <- lab_data %>% drop_na(lab_neutrophils) %>% filter(lab_neutrophils > 0 & !is.infinite(lab_neutrophils))
p2 <- ggplot(tmp, aes(x=study, y=lab_neutrophils, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_neutrophils")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_monocytes.png", height=700)
lab_data <- lab_data %>% mutate(lab_monocytes = case_when(lab_monocytes == 0 & !is.na(lab_monocytes) ~ min(lab_data$lab_monocytes[lab_data$lab_monocytes != 0], na.rm=T)/2,
                                                            TRUE ~ lab_monocytes))
tmp <- lab_data %>% drop_na(lab_monocytes) %>% filter(lab_monocytes > 0 & !is.infinite(lab_monocytes))
tmp$lab_monocytes <- log(tmp$lab_monocytes)
LM <- glm(lab_monocytes ~ study, data=tmp, family="gaussian")
tmp$lab_monocytes <- residuals(LM)
tmp1 <- outlier(tmp$lab_monocytes,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_monocytes, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_monocytes (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_monocytes = ifelse(outlier==TRUE, -1, lab_monocytes))
tmp <- lab_data %>% drop_na(lab_monocytes) %>% filter(lab_monocytes > 0 & !is.infinite(lab_monocytes))
p2 <- ggplot(tmp, aes(x=study, y=lab_monocytes, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_monocytes")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_platelets.png", height=700)
lab_data <- lab_data %>% mutate(lab_platelets = case_when(lab_platelets == 0 & !is.na(lab_platelets) ~ min(lab_data$lab_platelets[lab_data$lab_platelets != 0], na.rm=T)/2,
                                                            TRUE ~ lab_platelets))
tmp <- lab_data %>% drop_na(lab_platelets) %>% filter(lab_platelets > 0 & !is.infinite(lab_platelets))
tmp$lab_platelets <- log(tmp$lab_platelets)
LM <- glm(lab_platelets ~ study, data=tmp, family="gaussian")
tmp$lab_platelets <- residuals(LM)
tmp1 <- outlier(tmp$lab_platelets,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_platelets, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_platelets (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_platelets = ifelse(outlier==TRUE, -1, lab_platelets))
tmp <- lab_data %>% drop_na(lab_platelets) %>% filter(lab_platelets > 0 & !is.infinite(lab_platelets))
p2 <- ggplot(tmp, aes(x=study, y=lab_platelets, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_platelets")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_eosinophils.png", height=700)
lab_data <- lab_data %>% mutate(lab_eosinophils = case_when(lab_eosinophils == 0 & !is.na(lab_eosinophils) ~ min(lab_data$lab_eosinophils[lab_data$lab_eosinophils != 0], na.rm=T)/2,
                                                            TRUE ~ lab_eosinophils))
tmp <- lab_data %>% drop_na(lab_eosinophils) %>% filter(lab_eosinophils > 0 & !is.infinite(lab_eosinophils))
tmp$lab_eosinophils <- log(tmp$lab_eosinophils)
LM <- glm(lab_eosinophils ~ study, data=tmp, family="gaussian")
tmp$lab_eosinophils <- residuals(LM)
tmp1 <- outlier(tmp$lab_eosinophils,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_eosinophils, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_eosinophils (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_eosinophils = ifelse(outlier==TRUE, -1, lab_eosinophils))
tmp <- lab_data %>% drop_na(lab_eosinophils) %>% filter(lab_eosinophils > 0 & !is.infinite(lab_eosinophils))
p2 <- ggplot(tmp, aes(x=study, y=lab_eosinophils, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_eosinophils")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_basophils.png", height=700)
lab_data <- lab_data %>% mutate(lab_basophils = case_when(lab_basophils == 0 & !is.na(lab_basophils) ~ min(lab_data$lab_basophils[lab_data$lab_basophils != 0], na.rm=T)/2,
                                                            TRUE ~ lab_basophils))
tmp <- lab_data %>% drop_na(lab_basophils) %>% filter(lab_basophils > 0 & !is.infinite(lab_basophils))
tmp$lab_basophils <- log(tmp$lab_basophils)
LM <- glm(lab_basophils ~ study, data=tmp, family="gaussian")
tmp$lab_basophils <- residuals(LM)
tmp1 <- outlier(tmp$lab_basophils,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_basophils, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_basophils (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_basophils = ifelse(outlier==TRUE, -1, lab_basophils))
tmp <- lab_data %>% drop_na(lab_basophils) %>% filter(lab_basophils > 0 & !is.infinite(lab_basophils))
p2 <- ggplot(tmp, aes(x=study, y=lab_basophils, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_basophils")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_crp.png", height=700)
lab_data <- lab_data %>% mutate(lab_crp = case_when(lab_crp == 0 & !is.na(lab_crp) ~ min(lab_data$lab_crp[lab_data$lab_crp != 0], na.rm=T)/2,
                                                            TRUE ~ lab_crp))
tmp <- lab_data %>% drop_na(lab_crp) %>% filter(lab_crp > 0 & !is.infinite(lab_crp))
tmp$lab_crp <- log(tmp$lab_crp)
LM <- glm(lab_crp ~ study, data=tmp, family="gaussian")
tmp$lab_crp <- residuals(LM)
tmp1 <- outlier(tmp$lab_crp,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_crp, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_crp (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_crp = ifelse(outlier==TRUE, -1, lab_crp))
tmp <- lab_data %>% drop_na(lab_crp) %>% filter(lab_crp > 0 & !is.infinite(lab_crp))
p2 <- ggplot(tmp, aes(x=study, y=lab_crp, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_crp")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_trop_t.png", height=700)
lab_data <- lab_data %>% mutate(lab_trop_t = case_when(lab_trop_t == 0 & !is.na(lab_trop_t) ~ min(lab_data$lab_trop_t[lab_data$lab_trop_t != 0], na.rm=T)/2,
                                                            TRUE ~ lab_trop_t))
tmp <- lab_data %>% drop_na(lab_trop_t) %>% filter(lab_trop_t > 0 & !is.infinite(lab_trop_t))
tmp$lab_trop_t <- log(tmp$lab_trop_t)
LM <- glm(lab_trop_t ~ study, data=tmp, family="gaussian")
tmp$lab_trop_t <- residuals(LM)
tmp1 <- outlier(tmp$lab_trop_t,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_trop_t, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_trop_t (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_trop_t = ifelse(outlier==TRUE, -1, lab_trop_t))
tmp <- lab_data %>% drop_na(lab_trop_t) %>% filter(lab_trop_t > 0 & !is.infinite(lab_trop_t))
p2 <- ggplot(tmp, aes(x=study, y=lab_trop_t, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_trop_t")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_ast.png", height=700)
lab_data <- lab_data %>% mutate(lab_ast = case_when(lab_ast == 0 & !is.na(lab_ast) ~ min(lab_data$lab_ast[lab_data$lab_ast != 0], na.rm=T)/2,
                                                            TRUE ~ lab_ast))
tmp <- lab_data %>% drop_na(lab_ast) %>% filter(lab_ast > 0 & !is.infinite(lab_ast))
tmp$lab_ast <- log(tmp$lab_ast)
LM <- glm(lab_ast ~ study, data=tmp, family="gaussian")
tmp$lab_ast <- residuals(LM)
tmp1 <- outlier(tmp$lab_ast,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_ast, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_ast (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_ast = ifelse(outlier==TRUE, -1, lab_ast))
tmp <- lab_data %>% drop_na(lab_ast) %>% filter(lab_ast > 0 & !is.infinite(lab_ast))
p2 <- ggplot(tmp, aes(x=study, y=lab_ast, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_ast")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_alt.png", height=700)
lab_data <- lab_data %>% mutate(lab_alt = case_when(lab_alt == 0 & !is.na(lab_alt) ~ min(lab_data$lab_alt[lab_data$lab_alt != 0], na.rm=T)/2,
                                                            TRUE ~ lab_alt))
tmp <- lab_data %>% drop_na(lab_alt) %>% filter(lab_alt > 0 & !is.infinite(lab_alt))
tmp$lab_alt <- log(tmp$lab_alt)
LM <- glm(lab_alt ~ study, data=tmp, family="gaussian")
tmp$lab_alt <- residuals(LM)
tmp1 <- outlier(tmp$lab_alt,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_alt, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_alt (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_alt = ifelse(outlier==TRUE, -1, lab_alt))
tmp <- lab_data %>% drop_na(lab_alt) %>% filter(lab_alt > 0 & !is.infinite(lab_alt))
p2 <- ggplot(tmp, aes(x=study, y=lab_alt, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_alt")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_bilirubin.png", height=700)
lab_data <- lab_data %>% mutate(lab_bilirubin = case_when(lab_bilirubin == 0 & !is.na(lab_bilirubin) ~ min(lab_data$lab_bilirubin[lab_data$lab_bilirubin != 0], na.rm=T)/2,
                                                            TRUE ~ lab_bilirubin))
tmp <- lab_data %>% drop_na(lab_bilirubin) %>% filter(lab_bilirubin > 0 & !is.infinite(lab_bilirubin))
tmp$lab_bilirubin <- log(tmp$lab_bilirubin)
LM <- glm(lab_bilirubin ~ study, data=tmp, family="gaussian")
tmp$lab_bilirubin <- residuals(LM)
tmp1 <- outlier(tmp$lab_bilirubin,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_bilirubin, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_bilirubin (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_bilirubin = ifelse(outlier==TRUE, -1, lab_bilirubin))
tmp <- lab_data %>% drop_na(lab_bilirubin) %>% filter(lab_bilirubin > 0 & !is.infinite(lab_bilirubin))
p2 <- ggplot(tmp, aes(x=study, y=lab_bilirubin, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_bilirubin")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_ldh.png", height=700)
lab_data <- lab_data %>% mutate(lab_ldh = case_when(lab_ldh == 0 & !is.na(lab_ldh) ~ min(lab_data$lab_ldh[lab_data$lab_ldh != 0], na.rm=T)/2,
                                                            TRUE ~ lab_ldh))
tmp <- lab_data %>% drop_na(lab_ldh) %>% filter(lab_ldh > 0 & !is.infinite(lab_ldh))
tmp$lab_ldh <- log(tmp$lab_ldh)
LM <- glm(lab_ldh ~ study, data=tmp, family="gaussian")
tmp$lab_ldh <- residuals(LM)
tmp1 <- outlier(tmp$lab_ldh,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_ldh, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_ldh (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_ldh = ifelse(outlier==TRUE, -1, lab_ldh))
tmp <- lab_data %>% drop_na(lab_ldh) %>% filter(lab_ldh > 0 & !is.infinite(lab_ldh))
p2 <- ggplot(tmp, aes(x=study, y=lab_ldh, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_ldh")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_ggt.png", height=700)
lab_data <- lab_data %>% mutate(lab_ggt = case_when(lab_ggt == 0 & !is.na(lab_ggt) ~ min(lab_data$lab_ggt[lab_data$lab_ggt != 0], na.rm=T)/2,
                                                            TRUE ~ lab_ggt))
tmp <- lab_data %>% drop_na(lab_ggt) %>% filter(lab_ggt > 0 & !is.infinite(lab_ggt))
tmp$lab_ggt <- log(tmp$lab_ggt)
LM <- glm(lab_ggt ~ study, data=tmp, family="gaussian")
tmp$lab_ggt <- residuals(LM)
tmp1 <- outlier(tmp$lab_ggt,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_ggt, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_ggt (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_ggt = ifelse(outlier==TRUE, -1, lab_ggt))
tmp <- lab_data %>% drop_na(lab_ggt) %>% filter(lab_ggt > 0 & !is.infinite(lab_ggt))
p2 <- ggplot(tmp, aes(x=study, y=lab_ggt, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_ggt")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_alp.png", height=700)
lab_data <- lab_data %>% mutate(lab_alp = case_when(lab_alp == 0 & !is.na(lab_alp) ~ min(lab_data$lab_alp[lab_data$lab_alp != 0], na.rm=T)/2,
                                                            TRUE ~ lab_alp))
tmp <- lab_data %>% drop_na(lab_alp) %>% filter(lab_alp > 0 & !is.infinite(lab_alp))
tmp$lab_alp <- log(tmp$lab_alp)
LM <- glm(lab_alp ~ study, data=tmp, family="gaussian")
tmp$lab_alp <- residuals(LM)
tmp1 <- outlier(tmp$lab_alp,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_alp, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_alp (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_alp = ifelse(outlier==TRUE, -1, lab_alp))
tmp <- lab_data %>% drop_na(lab_alp) %>% filter(lab_alp > 0 & !is.infinite(lab_alp))
p2 <- ggplot(tmp, aes(x=study, y=lab_alp, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_alp")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_d_dimer.png", height=700)
lab_data <- lab_data %>% mutate(lab_d_dimer = case_when(lab_d_dimer == 0 & !is.na(lab_d_dimer) ~ min(lab_data$lab_d_dimer[lab_data$lab_d_dimer != 0], na.rm=T)/2,
                                                            TRUE ~ lab_d_dimer))
tmp <- lab_data %>% drop_na(lab_d_dimer) %>% filter(lab_d_dimer > 0 & !is.infinite(lab_d_dimer))
tmp$lab_d_dimer <- log(tmp$lab_d_dimer)
LM <- glm(lab_d_dimer ~ study, data=tmp, family="gaussian")
tmp$lab_d_dimer <- residuals(LM)
tmp1 <- outlier(tmp$lab_d_dimer,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_d_dimer, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_d_dimer (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_d_dimer = ifelse(outlier==TRUE, -1, lab_d_dimer))
tmp <- lab_data %>% drop_na(lab_d_dimer) %>% filter(lab_d_dimer > 0 & !is.infinite(lab_d_dimer))
p2 <- ggplot(tmp, aes(x=study, y=lab_d_dimer, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_d_dimer")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_il_6.png", height=700)
lab_data <- lab_data %>% mutate(lab_il_6 = case_when(lab_il_6 == 0 & !is.na(lab_il_6) ~ min(lab_data$lab_il_6[lab_data$lab_il_6 != 0], na.rm=T)/2,
                                                            TRUE ~ lab_il_6))
tmp <- lab_data %>% drop_na(lab_il_6) %>% filter(lab_il_6 > 0 & !is.infinite(lab_il_6))
tmp$lab_il_6 <- log(tmp$lab_il_6)
LM <- glm(lab_il_6 ~ study, data=tmp, family="gaussian")
tmp$lab_il_6 <- residuals(LM)
tmp1 <- outlier(tmp$lab_il_6,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_il_6, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_il_6 (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_il_6 = ifelse(outlier==TRUE, -1, lab_il_6))
tmp <- lab_data %>% drop_na(lab_il_6) %>% filter(lab_il_6 > 0 & !is.infinite(lab_il_6))
p2 <- ggplot(tmp, aes(x=study, y=lab_il_6, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_il_6")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_serum_ferritin.png", height=700)
lab_data <- lab_data %>% mutate(lab_serum_ferritin = case_when(lab_serum_ferritin == 0 & !is.na(lab_serum_ferritin) ~ min(lab_data$lab_serum_ferritin[lab_data$lab_serum_ferritin != 0], na.rm=T)/2,
                                                            TRUE ~ lab_serum_ferritin))
tmp <- lab_data %>% drop_na(lab_serum_ferritin) %>% filter(lab_serum_ferritin > 0 & !is.infinite(lab_serum_ferritin))
tmp$lab_serum_ferritin <- log(tmp$lab_serum_ferritin)
LM <- glm(lab_serum_ferritin ~ study, data=tmp, family="gaussian")
tmp$lab_serum_ferritin <- residuals(LM)
tmp1 <- outlier(tmp$lab_serum_ferritin,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_serum_ferritin, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_serum_ferritin (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_serum_ferritin = ifelse(outlier==TRUE, -1, lab_serum_ferritin))
tmp <- lab_data %>% drop_na(lab_serum_ferritin) %>% filter(lab_serum_ferritin > 0 & !is.infinite(lab_serum_ferritin))
p2 <- ggplot(tmp, aes(x=study, y=lab_serum_ferritin, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_serum_ferritin")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)


png("plots/labs/lab_procalcitonin.png", height=700)
lab_data <- lab_data %>% mutate(lab_procalcitonin = case_when(lab_procalcitonin == 0 & !is.na(lab_procalcitonin) ~ min(lab_data$lab_procalcitonin[lab_data$lab_procalcitonin != 0], na.rm=T)/2,
                                                            TRUE ~ lab_procalcitonin))
tmp <- lab_data %>% drop_na(lab_procalcitonin) %>% filter(lab_procalcitonin > 0 & !is.infinite(lab_procalcitonin))
tmp$lab_procalcitonin <- log(tmp$lab_procalcitonin)
LM <- glm(lab_procalcitonin ~ study, data=tmp, family="gaussian")
tmp$lab_procalcitonin <- residuals(LM)
tmp1 <- outlier(tmp$lab_procalcitonin,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_procalcitonin, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_procalcitonin (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_procalcitonin = ifelse(outlier==TRUE, -1, lab_procalcitonin))
tmp <- lab_data %>% drop_na(lab_procalcitonin) %>% filter(lab_procalcitonin > 0 & !is.infinite(lab_procalcitonin))
p2 <- ggplot(tmp, aes(x=study, y=lab_procalcitonin, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_procalcitonin")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_ck.png", height=700)
lab_data <- lab_data %>% mutate(lab_ck = case_when(lab_ck == 0 & !is.na(lab_ck) ~ min(lab_data$lab_ck[lab_data$lab_ck != 0], na.rm=T)/2,
                                                            TRUE ~ lab_ck))
tmp <- lab_data %>% drop_na(lab_ck) %>% filter(lab_ck > 0 & !is.infinite(lab_ck))
tmp$lab_ck <- log(tmp$lab_ck)
LM <- glm(lab_ck ~ study, data=tmp, family="gaussian")
tmp$lab_ck <- residuals(LM)
tmp1 <- outlier(tmp$lab_ck,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_ck, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_ck (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_ck = ifelse(outlier==TRUE, -1, lab_ck))
tmp <- lab_data %>% drop_na(lab_ck) %>% filter(lab_ck > 0 & !is.infinite(lab_ck))
p2 <- ggplot(tmp, aes(x=study, y=lab_ck, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_ck")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_fibrinogen.png", height=700)
lab_data <- lab_data %>% mutate(lab_fibrinogen = case_when(lab_fibrinogen == 0 & !is.na(lab_fibrinogen) ~ min(lab_data$lab_fibrinogen[lab_data$lab_fibrinogen != 0], na.rm=T)/2,
                                                            TRUE ~ lab_fibrinogen))
tmp <- lab_data %>% drop_na(lab_fibrinogen) %>% filter(lab_fibrinogen > 0 & !is.infinite(lab_fibrinogen))
tmp$lab_fibrinogen <- log(tmp$lab_fibrinogen)
LM <- glm(lab_fibrinogen ~ study, data=tmp, family="gaussian")
tmp$lab_fibrinogen <- residuals(LM)
tmp1 <- outlier(tmp$lab_fibrinogen,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_fibrinogen, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_fibrinogen (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_fibrinogen = ifelse(outlier==TRUE, -1, lab_fibrinogen))
tmp <- lab_data %>% drop_na(lab_fibrinogen) %>% filter(lab_fibrinogen > 0 & !is.infinite(lab_fibrinogen))
p2 <- ggplot(tmp, aes(x=study, y=lab_fibrinogen, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_fibrinogen")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_creatinine.png", height=700)
lab_data <- lab_data %>% mutate(lab_creatinine = case_when(lab_creatinine == 0 & !is.na(lab_creatinine) ~ min(lab_data$lab_creatinine[lab_data$lab_creatinine != 0], na.rm=T)/2,
                                                            TRUE ~ lab_creatinine))
tmp <- lab_data %>% drop_na(lab_creatinine) %>% filter(lab_creatinine > 0 & !is.infinite(lab_creatinine))
tmp$lab_creatinine <- log(tmp$lab_creatinine)
LM <- glm(lab_creatinine ~ study, data=tmp, family="gaussian")
tmp$lab_creatinine <- residuals(LM)
tmp1 <- outlier(tmp$lab_creatinine,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_creatinine, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_creatinine (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_creatinine = ifelse(outlier==TRUE, -1, lab_creatinine))
tmp <- lab_data %>% drop_na(lab_creatinine) %>% filter(lab_creatinine > 0 & !is.infinite(lab_creatinine))
p2 <- ggplot(tmp, aes(x=study, y=lab_creatinine, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_creatinine")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

png("plots/labs/lab_trop_i.png", height=700)
lab_data <- lab_data %>% mutate(lab_trop_i = case_when(lab_trop_i == 0 & !is.na(lab_trop_i) ~ min(lab_data$lab_trop_i[lab_data$lab_trop_i != 0], na.rm=T)/2,
                                                           TRUE ~ lab_trop_i))
tmp <- lab_data %>% drop_na(lab_trop_i) %>% filter(lab_trop_i > 0 & !is.infinite(lab_trop_i))
tmp$lab_trop_i <- log(tmp$lab_trop_i)
LM <- glm(lab_trop_i ~ study, data=tmp, family="gaussian")
tmp$lab_trop_i <- residuals(LM)
tmp1 <- outlier(tmp$lab_trop_i,method="median",addthres=TRUE)
tmp <- cbind(tmp, tmp1$flaggedData[,3])
colnames(tmp)[dim(tmp)[2]] <- "outlier"
p1 <- ggplot(tmp, aes(x=anonymized_patient_id,y=lab_trop_i, color=study, shape=outlier)) + geom_point() +
  geom_hline(yintercept=tmp1$lowerbound, linetype = "dashed") + xlab("") +
  geom_hline(yintercept=tmp1$upperbound, linetype = "dashed") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_shape_manual(values = c(1,8)) + ggtitle(paste0("lab_trop_i (outliers=", tmp1$outlN,"[",round(tmp1$outlN/dim(tmp1$flaggedData)[1]*100,2),"%])"))
lab_data <- lab_data %>% merge(tmp[,c(1,dim(tmp)[2])], by="index", all.x=T) 
lab_data <- lab_data %>% mutate(lab_trop_i = ifelse(outlier==TRUE, -1, lab_trop_i))
tmp <- lab_data %>% drop_na(lab_trop_i) %>% filter(lab_trop_i > 0 & !is.infinite(lab_trop_i))
p2 <- ggplot(tmp, aes(x=study, y=lab_trop_i, color=study)) + geom_jitter() + geom_violin(trim=FALSE)+
  theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ggtitle("Post-QC lab_trop_i")
grid.arrange(p1,p2)
dev.off()

lab_data <- lab_data %>% select(-outlier)

lab_data %>% select(c(-index,-study)) %>% saveRDS(file="curated_clinical/all_lab_data.rds")

lab_data <- readRDS("curated_clinical/all_lab_data.rds")

blood_values <- c("lab_result_date",
                  "lab_wbc",
                  "lab_lymphocytes",
                  "lab_cd4",
                  "lab_cd8",
                  "lab_neutrophils",
                  "lab_monocytes",
                  "lab_platelets",
                  "lab_eosinophils",
                  "lab_basophils",
                  "lab_crp",
                  "lab_trop_t",
                  "lab_trop_i",
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
                  "lab_aptt",
                  "lab_inr",
                  "lab_creatinine",
                  "lab_nk", "lab_n_l_ratio","lab_na","lab_k", "lab_mcv","lab_lipase", "lab_amylase")

lab_data <- lab_data %>% mutate_at(.vars = c(blood_values), .funs = funs(ifelse(.==-1,NA,.)))

lab_data <- readRDS("curated_clinical/all_lab_data.rds")
