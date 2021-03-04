setwd("~/02.associations")
out <- fread("results/AUC.tsv")
tmp <- out %>% filter(outcome %in% c("a1")) %>%
  filter(agerange %in% c("ALL", "Age≤60")) %>%
  filter(model %in% c("age+sex+allriskfactors", "age+sex+allriskfactors+rs10490770"))
tmp$AUC <- as.numeric(tmp$AUC)

changeSciNot <- function(n) {
  output <- format(n, scientific = TRUE) #Transforms the number into scientific notation even if small
  output <- sub("e", "10x", output) #Replace e with 10^
  output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
  output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
  output
}

tmp <- tmp %>% mutate(AUC = paste0(round(AUC,2)," [",round(AUC.LL,2),"; ",round(AUC.UL,2),"]"),
                      `AUC p-value` = changeSciNot(`AUC.p-value`))
tmp <- tmp %>% mutate(NRI = paste0(round(NRI,2)," [",round(NRI.LL,2),"; ",round(NRI.UL,2),"]"),
                      `NRI p-value` = changeSciNot(`NRI.p-value`))
tmp$`AUC p-value`[c(1,3)] <- NA
tmp$`NRI p-value`[c(1,3)] <- NA
tmp$`NRI`[c(1,3)] <- NA

tmp <- tmp %>% mutate(Outcome = ifelse(outcome=="a1","Death or severe respiratory failure", "Severe respiratory failure"))

tmp <- tmp %>% select(c("agerange","model","AUC","AUC p-value","NRI","NRI p-value"))
tmp$model <- rep(c("Baseline","Baseline and rs10490770"),2)
tmp %>% write.xlsx("results/Table3.xlsx",row.names = F)
