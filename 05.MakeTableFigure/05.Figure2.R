setwd("/home/tomoko/02.associations/")

out1 <- fread("results/clinical_EUR.results") %>% mutate(Group = "All")
out1$Outcome[3] <- "icu_admit"
out2 <- fread("results/clinical_EUR_hosp.results") %>% mutate(Group = "Hospitalized only")

OUT <- bind_rows(out1, out2)

out <- OUT %>% filter(Outcome %in% c("hospitalization","a1", "resp_severe", "icu_admit", "vte", "cardiac", "aki", "hepatic"))

library(meta)

tmp <- out %>% mutate(order = case_when(Outcome == "hospitalization" ~ 1,
                                        Outcome == "icu_admit" ~ 2,
                                        Outcome == "a1" ~ 3,
                                        Outcome == "resp_severe" ~ 4,
                                        Outcome == "vte" ~ 8,
                                        Outcome == "hepatic" ~ 5,
                                        Outcome == "cardiac" ~ 6,
                                        Outcome == "aki" ~ 7))

tmp <- tmp %>% mutate(Outcome = case_when(Outcome == "hospitalization" ~ paste0("Hospitalization (Ncase=",N_case,")"),
                                          Outcome == "icu_admit" ~  paste0("ICU admission (Ncase=",N_case,")"),
                                          Outcome == "a1" ~  paste0("Death or severe respiratory failure (Ncase=",N_case,")"),
                                          Outcome == "resp_severe" ~  paste0("Severe respiratory failure (Ncase=",N_case,")"),
                                          Outcome == "vte" ~  paste0("Venous thromboembolism (Ncase=",N_case,")"),
                                          Outcome == "hepatic" ~  paste0("Hepatic injury (Ncase=",N_case,")"),
                                          Outcome == "cardiac" ~  paste0("Cardiovascular complications (Ncase=",N_case,")"),
                                          Outcome == "aki" ~  paste0("Kidney injury (Ncase=",N_case,")")))
tmp$Outcome[tmp$Outcome == "Cardiovascular complications (Ncase=694)"] <- "Cardiovascular complications (Ncase=697)"
tmp <- tmp[order(tmp$order, decreasing=F),]
tmp <- tmp %>% mutate(OR = exp(beta),
                      LL = exp(beta + qnorm(0.025)*se),
                      UL = exp(beta + qnorm(0.975)*se))

tmp$Outcome <- factor(tmp$Outcome, levels=c(unique(tmp$Outcome)))

png("plots/Figure2.png", width=1000, height=600)
p1 <- ggplot(tmp, aes(x=Group, y=OR, ymin=LL, ymax=UL, color=Group)) +
  geom_pointrange(aes(col=Group), lwd=0.8) + geom_hline(aes(fill=Group), yintercept =1, linetype=2) +
  xlab("") + ylab("Odds Ratio (95% Confidence Interval)") + geom_text(aes(label=round(OR,1), y=OR, col=Group, hjust = 0.5, vjust = 1.3), size=5) +
  geom_errorbar(aes(ymin=LL, ymax=UL, col=Group), width=0.1, cex=1) + ylim(0.8,3.5)+
  facet_wrap(~Outcome,  strip.position = 'left', nrow = 10) + theme_minimal() +
  scale_color_brewer(palette = "Set1") + labs(color="Group") +
  theme(plot.title=element_text(size=15),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15,face="bold"), 
        legend.title = element_text(size=15,face="bold"), 
        strip.text.y.left = element_text(hjust=0.5,vjust = 0.5,angle=0,size=15,face="bold"))+
  coord_flip() +guides(col = guide_legend(reverse = TRUE))
gt = ggplotGrob(p1)
gt$heights[16] = unit(0.5, "cm")
grid.newpage()
grid.draw(gt)
dev.off()
