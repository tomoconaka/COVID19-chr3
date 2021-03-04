setwd("/home/tomoko/02.associations/")

out1 <- fread("results/multivariate.table.tsv")
colnames(out1) <- c("outcome","agerange","riskfactor","OR","LL","UL","p","freq")

tmp <- out1 %>% filter(outcome %in% c("a1")) 
tmp <- tmp %>% filter(agerange %in% c("ALL","Age≤60")) %>% filter(riskfactor != "age_at_diagnosis")

tmp <- tmp %>% mutate(agerange = c(rep("All",10),rep("Age≤60",10)),
                      riskfactor = case_when(riskfactor=="snp" ~ paste0("rs10490770 risk allele carrier"),
                                             riskfactor=="Male" ~ paste0("Male"),
                                             riskfactor=="smoke" ~ paste0("Smoking"),
                                             riskfactor=="BMIhigh" ~ paste0("BMI≥30"),
                                             riskfactor=="com_cancer" ~ paste0("Cancer"),
                                             riskfactor=="com_chronic_kidney" ~ paste0("CKD"),
                                             riskfactor=="com_chronic_pulm" ~ paste0("COPD"),
                                             riskfactor=="com_transplant" ~ paste0("Transplantation"),
                                             riskfactor=="com_heart_failure" ~ paste0("CHF"),
                                             riskfactor=="com_diabetes" ~ paste0("DM")
                                             ))

tmp$riskfactor <- factor(tmp$riskfactor, levels=c("Male", "Smoking","BMI≥30","Cancer","CKD","COPD","Transplantation","CHF","DM","rs10490770 risk allele carrier"))
tmp <- tmp %>% mutate(group = case_when(riskfactor=="rs10490770 risk allele carrier" ~ "genotype",
                                        riskfactor %in% c("Male","Smoking","BMI≥30") ~ "demographics",
                                        TRUE ~ "Comorbidities"
                                        ))

tmp1 <- tmp %>% filter(outcome=="a1")
png("plots/Figure4A.png", width=700, height=900)
p1 <- ggplot(tmp1, aes(x=riskfactor, y=OR, ymin=LL, ymax=UL, color=group)) +
  geom_pointrange(aes(col=group, fatten=5*log10(100*freq)), lwd=0.8) + geom_hline(aes(fill=group), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") + geom_text(aes(label=round(OR,1),y=OR, hjust = 0.5, vjust = 1.5, color=group),size=7) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=group), width=0.1, cex=1) +
  facet_wrap(~agerange, strip.position="bottom", nrow=2) + theme_classic() +
  labs(color="") + ggtitle("") + scale_color_manual(values=c("#4E84C4","#52854C","#D16103")) + 
  theme(text = element_text(size=15,face="bold"),
        axis.text.x = element_text(face="bold", 
                                   size=15),
        axis.text.y = element_text(face="bold", 
                                   size=15),
        legend.position = "none") + scale_y_log10() + coord_flip()
grid.arrange(p1)
dev.off()

out2 <- fread("results/AUC.tsv")
tmp <- out2 %>% filter(outcome %in% c("a1")) 
tmp <- tmp %>% filter(agerange %in% c("ALL","Age≤60")) %>% filter(model != "age+sex+allriskfactors" & model != "age+sex+allriskfactors+rs10490770")
tmp <- tmp %>% mutate(modelgroup=case_when(model %in% c("age+sex", "age+sex+smoking","age+sex+BMI") ~ "demographics",
                                      model == "age+sex+rs10490770" ~ "genotype",
                                      TRUE ~ "comorbidities"))
tmp$agerange[tmp$agerange == "ALL"] <- "All"
tmp$model[tmp$model  == "age+sex+diabetes"] <- "age+sex+DM"
tmp$model <- factor(tmp$model, levels=unique(tmp$model))
png("plots/Figure4B.png", width=1000, height=800)
tmp1 <- tmp %>% filter(outcome %in% c("a1")) 
p1 <- ggplot(tmp1, aes(x=model, y=AUC, ymin=AUC.LL, ymax=AUC.UL, fill=modelgroup)) +
  geom_bar(stat="identity", position=position_dodge()) + facet_wrap(~agerange, strip.position="bottom", nrow=2) +
  geom_errorbar(width=.2,
                position=position_dodge(.9)) + coord_cartesian(ylim=c(0.5,0.9)) + 
  xlab("") + ylab("AUC (95% Confidence Interval)") + geom_text(aes(y=AUC.UL+0.03, label=round(AUC,2), hjust = 0.5,vjust=0),size=6, position=position_dodge(.9))+
  theme_classic() + ggtitle("") +
  scale_fill_manual(values=c("#4E84C4","#52854C","#D16103")) + labs(fill="Model group")  +
  theme(text = element_text(size=15, face="bold"),legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1, size=15, face="bold")) + guides(color = guide_legend(legend.position = "None")) 
p1
dev.off()
