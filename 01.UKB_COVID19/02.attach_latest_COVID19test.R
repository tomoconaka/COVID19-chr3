setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS")
data <- readRDS("../../../data/01.UKBB/Recruited20201023.rds")
test1 <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/covid19_result_20210104.txt.gz")

t1 <- test1
t2 <- data.frame(unique(t1$eid))
colnames(t2) <- "eid"
t2$status <- 0
t2$COVIDdate <- 0
t2$hospital <- 0
for(i in seq(1,dim(t2)[1])){
  tmp <- t1[t1$eid == t2$eid[i],]
  a <- ifelse(any(tmp$result == 1), 1, 0)
  b <- ifelse(any(tmp$result == 1), min(as.Date(tmp$specdate[tmp$result == 1], format="%d/%m/%Y")), max(as.Date(tmp$specdate[tmp$result == 0], format="%d/%m/%Y")))
  c <- ifelse(any(tmp$hosaq[tmp$result == 1] == 1) | any(tmp$reqorg[tmp$result == 1] == 1), 1, 0) 
  t2$status[i] <- a
  t2$COVIDdate[i] <- as.Date(b, origin="1970-01-01")
  t2$hospital[i] <- c
}

t2$COVIDdate <- as.Date(t2$COVIDdate, origin="1970-01-01")
colnames(t2)[1] <- "ID"
res <- merge(data, t2, by="ID", all.x=T)
res$status[is.na(res$status)] <- 0
res$COVIDdate[is.na(res$COVIDdate)] <- "2020-12-21"
res$hospital[is.na(res$hospital)] <- 0
res$Duration <- round(as.numeric(as.Date(res$COVIDdate) - as.Date(res$`1stDATE`))/365)
res$newage <- res$AGE + res$Duration

res <- res %>% select(c("ID","SEX","CENTRE","date_of_death","COVIDdate","hospital", "newage", "status")) %>% rename(AGE = newage)
saveRDS(res,file="/home/richards/tomoko.nakanishi/09.COVID19/src/01.UKBB/01.GWAS/data/UKBall20210104.457941.rds")

