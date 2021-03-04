setwd("/home/tomoko/02.associations/")
out3 <- fread("results/clinical_EUR_age_young60.results") %>% mutate(Group = "All") %>% mutate(AgeGroup = "Age≤60")
out4 <- fread("results/clinical_EUR_age_young60_hosp.results") %>% mutate(Group = "Hospitalized only") %>% mutate(AgeGroup = "Age≤60")
out5 <- fread("results/clinical_EUR_age_old60.results") %>%  mutate(Group = "All") %>% mutate(AgeGroup = "Age>60")
out6 <- fread("results/clinical_EUR_age_old60_hosp.results") %>%  mutate(Group = "Hospitalized only") %>% mutate(AgeGroup = "Age>60")

OUT <- bind_rows(out3, out4, out5, out6)

out <- OUT %>% filter(Outcome %in% c("a1", "hospitalization","icu_admit","icu_admit_rev"))
out$Outcome[out$Outcome == "icu_admit_rev"] <- "icu_admit"
tmp <- out %>% mutate(order = case_when(Outcome == "hospitalization" ~ 1,
                                        Outcome == "icu_admit" ~ 2,
                                        Outcome == "a1" ~ 3
))

tmp <- tmp %>% mutate(Outcome = case_when(Outcome == "a1" ~ "Death or severe respiratory failure",
                                          Outcome == "hospitalization" ~ "Hospitalization",
                                          Outcome == "icu_admit" ~ "ICU admission"
))
tmp <- tmp[order(tmp$order, decreasing=F),]
tmp <- tmp %>% mutate(OR = exp(beta),
                      LL = exp(beta + qnorm(0.025)*se),
                      UL = exp(beta + qnorm(0.975)*se))

#tmp$AgeGroup <- factor(tmp$AgeGroup, levels = rev(levels(tmp$AgeGroup)))
#levels(tmp$AgeGroup)[levels(tmp$AgeGroup)=="2"] <- "Age>60"


require(gridExtra)
require(ggrepel)
tmp1 <- tmp %>% filter(Outcome =="Death or severe respiratory failure") 
plot1 <- ggplot(tmp1, aes(x=rev(tmp1$AgeGroup), y=OR, ymin=LL, ymax=UL, color=rev(tmp1$AgeGroup))) +
  geom_pointrange(aes(col=AgeGroup), lwd=0.8) + geom_hline(aes(fill=AgeGroup), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") + geom_text(aes(label=round(OR,1), y=OR, col=AgeGroup, hjust = 1.2, vjust = 0), size=8,show.legend = FALSE) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=AgeGroup), width=0.1, cex=1) +
  facet_wrap(~Group, strip.position="bottom", nrow=1) + theme_classic() +
  scale_color_brewer(palette = "Dark2") + labs(color="Age Group", size=FALSE) + ggtitle("") +
  theme(text = element_text(size=20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

final <- readRDS("../01.data.QC/curated_clinical/complication_dat_final.rds")
snps <- colnames(final)[grepl("chr", colnames(final))][3:12]
outcome <- c("resp_severe","resp_mild", "icu_admit_rev" ,"vte", "cardiac",
             "hosp_bleeding", "hospitalization", "a1","aki","hepatic","hosp_thrombo")

final <- final %>% mutate_at(.vars = c(outcome), .funs = funs(ifelse(. == -1, NA, .)))
final <- final %>% mutate(age2 = age_at_diagnosis*age_at_diagnosis)
final <- final %>% mutate(Agegroup = case_when(age_at_diagnosis < 40 ~ "Age<40",
                                               age_at_diagnosis >= 40 & age_at_diagnosis < 50 ~ "Age40-50",
                                               age_at_diagnosis >= 50 & age_at_diagnosis < 60 ~ "0Age50-60",
                                               age_at_diagnosis >= 60 & age_at_diagnosis < 70 ~ "Age60-70",
                                               age_at_diagnosis >= 70 & age_at_diagnosis < 80 ~ "Age70-80",
                                               age_at_diagnosis >= 80 ~ "Age>80"),
                          Male = ifelse(sex == 0, 1, 0))

final <- final %>% mutate(smoke = case_when(smoking == 0 ~ "current",
                                            smoking == 1 ~ "ex",
                                            smoking == 2 ~ "0none",
                                            TRUE ~ "NA"),
                          BMI = case_when(weight/((height/100)^2) >= 40 ~ "class4",
                                          weight/((height/100)^2) < 40 & weight/((height/100)^2) >=35 ~ "class3",
                                          weight/((height/100)^2) < 35 & weight/((height/100)^2) >=30 ~ "class2",
                                          weight/((height/100)^2) < 30 ~ "0none",
                                          TRUE ~ "NA"))

final <- final %>% mutate_at(.vars = c("smoke","BMI",outcome), .funs = funs(ifelse(. == "NA", NA, .)))


#final <- final %>% filter(hospitalization == 1)
table(final[,c("study","icu_admit")])
table(final[,c("study","vte")])
final <- final %>% mutate(snp = ifelse(round(`chr3:45823240:T:C_C`) >= 1, 1, 0))
#final <- final %>% mutate_at(.vars=c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"), .funs = funs(scale(.)))
data_EUR <- final %>% filter(pop == "EUR")

data_EUR <- data_EUR %>% filter(a1 == 1) 
data_EUR <- data_EUR %>%
  mutate(snp = ifelse(snp == 1, paste0("carrier (N=",sum(data_EUR$snp == 1),")"), paste0("non-carrier (N=",sum(data_EUR$snp == 0),")")))
plot2 <- ggplot(data_EUR, aes(x=snp,y=age_at_diagnosis)) + geom_violin(aes(fill=snp), trim = FALSE) + 
  geom_boxplot(alpha=0.6,trim = FALSE) + 
  scale_fill_manual(values=c("#00AFBB", "#E7B800")) + ggtitle("") + labs(fill="rs10490770\nrisk allele") +
  theme_classic() + xlab("") + ylab("Age") + geom_text(aes(label=paste0("p=",round(OUT1$pvalue_age_snp[OUT1$Outcome == "a1"],5)),x=1.5, y=130),size=10) +
  geom_segment(aes(x = 1, y = 120, xend = 2, yend = 120)) +
  scale_y_continuous(breaks=c(30,40,50,60,70,80,90,100)) +
  theme(text = element_text(size=20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))

out1 <- fread("results/clinical_EUR_comorbidities.results") %>% mutate(Group="All")
out2 <- fread("results/clinical_EUR_comorbidities_hosp.results") %>% mutate(Group="Hospitalized only")
out <- bind_rows(out1, out2) %>% filter(Outcome %in% c("a1"))

tmp <- out %>% mutate(order = case_when(Outcome == "a1" ~ 1,
                                        Outcome == "resp_severe" ~ 2
))

tmp <- tmp %>% mutate(Outcome = case_when(Outcome == "a1" ~ "Death or severe respiratory failure",
                                          Outcome == "resp_severe" ~ "Severe respiratory failure"
))
tmp <- tmp[order(tmp$order, decreasing=F),]
tmp$N_comorb <- as.factor(tmp$N_comorb)

tmp <- tmp %>% mutate(OR = exp(beta),
                      LL = exp(beta + qnorm(0.025)*se),
                      UL = exp(beta + qnorm(0.975)*se))

levels(tmp$N_comorb)[levels(tmp$N_comorb)=="0"] <- "0"
levels(tmp$N_comorb)[levels(tmp$N_comorb)=="1"] <- "1"
levels(tmp$N_comorb)[levels(tmp$N_comorb)=="2"] <- "≥2"

plot3 <- ggplot(tmp, aes(x=N_comorb, y=OR, ymin=LL, ymax=UL, color=N_comorb)) +
  geom_pointrange(aes(col=N_comorb), lwd=0.8) + geom_hline(aes(fill=N_comorb), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") + geom_text(aes(label=round(OR,1), y=OR, col=N_comorb, hjust = 1.2,vjust = 0), size=8,show.legend = FALSE) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=N_comorb), width=0.1, cex=1) +
  facet_wrap(~Group, strip.position="bottom", nrow=1) + theme_classic() +
  scale_color_brewer(palette = "Set1") + labs(color="Number of risk factors") + ggtitle("") +
  theme(text = element_text(size=20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
grid.arrange(arrangeGrob(plot1,plot3,ncol=1),plot2,nrow=1)
dev.off()
png("plots/Figure3B.png", width=700, height=800)
layout1 <- rbind(c(1, 1, 1, 1,1,NA),
                 c(2, 2, 2, 2, 2,2))
grid.arrange(plot1,plot3,layout_matrix = layout1)
dev.off()
png("plots/Figure3A.png", width=800, height=800)
plot2
dev.off()
