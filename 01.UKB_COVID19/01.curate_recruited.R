setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS")
library(data.table)
infile <- "~/scratch/09.COVID19/data/01.UKBB/ukb27449_20688_fetch_data.tab.gz"
t <- fread(infile, sep="\t")
d <- data.frame(t$f.eid, t$f.31.0.0, t$f.54.0.0, t$f.53.0.0, t$f.21003.0.0, t$f.22000.0.0)
d$GTARRAY <- NA
d$GTARRAY[d$t.f.22000.0.0 < 0] <- "UKBiLEVE"
d$GTARRAY[d$t.f.22000.0.0 > 0] <- "UKAxiom"
d$GTARRAY <- ifelse(d$GTARRAY == "UKAxiom", 1, 0)
d1 <- d[-6]
colnames(d1)[1:5] <- c("ID", "SEX", "CENTRE","1stDATE","AGE")
qc <- fread("/mnt/RICHARDS_JBOD1/SHARE/DATA/GENOTYPE_DATA/UKB/PHENOVAR-27449/full_release/v2/genotyped/ukb2744_cal_v2_s488374_sqc_v2.txt")
qc1 <- qc[,-1]
colnames(qc1)[1] <- "ID"
dpqc <- merge(d1, qc1, by="ID")#488363
dpqc1 <- dpqc[dpqc$het.missing.outliers == 0,]
dpqc1 <- dpqc1[dpqc1$putative.sex.chromosome.aneuploidy == 0,]
dpqc1 <- dpqc1[dpqc1$Submitted.Gender==dpqc1$Inferred.Gender,]
dpqc1 <- dpqc1[dpqc1$in.Phasing.Input.chr1_22==1,]
dpqc1 <- dpqc1[dpqc1$in.Phasing.Input.chrX==1,]
dpqc1 <- dpqc1[dpqc1$in.Phasing.Input.chrXY==1,]#291
#486238

#withdrawal
w <- fread("/scratch/richards/tomoko.nakanishi/09.COVID19/data/01.UKBB/w27449_20200820.csv")
colnames(w) <- "ID"
dpqc1$withdrawal <- 0
dpqc1$withdrawal[dpqc1$ID %in% w$ID] <- 1
dpqc2 <- dpqc1[dpqc1$withdrawal == 0,]#486128
dpqc3 <- dpqc2 %>% mutate(SEX = SEX.x) %>% select(c("ID", "SEX", "CENTRE", "1stDATE", "AGE", "GTARRAY"))

## dead
death <- fread("../../../data/01.UKBB/death_20201021.txt.gz")
colnames(death)[1] <- "ID"
death <- death %>% mutate(date_of_death = as.Date(date_of_death, format="%d/%m/%Y"))
res1 <- merge(dpqc3, death, by="ID", all.x=T)
res2 <- res1 %>% filter(as.Date(date_of_death) > as.Date("2020-03-16") | is.na(date_of_death))#457941

saveRDS(res2, file="../../../data/01.UKBB/Recruited20201023.rds")
