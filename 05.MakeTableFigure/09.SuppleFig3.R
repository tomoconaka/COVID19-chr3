setwd("/home/tomoko/02.associations")

out1 <- fread("results/worst_EUR_biomarker.results") %>% filter(SNP == "chr3:45823240:T:C_C")
out1 <- out1 

out2 <- fread("results/worst_EUR_biomarker_clinical.results")

OUT <- merge(out1, out2, by="biomarker")

OUT <- OUT %>% mutate(sig = ifelse(pvalue.x < 0.05, TRUE, FALSE))

OUT <- OUT %>% mutate(LL.x=beta.x + qnorm(0.025)*se.x,
                      UL.x=beta.x + qnorm(0.975)*se.x,
                      LL.y=beta.y + qnorm(0.025)*se.y,
                      UL.y=beta.y + qnorm(0.975)*se.y
)

OUT$sig <- factor(OUT$sig, levels=c(rev(unique(OUT$sig))))

OUT <- OUT %>% mutate(biomarker=case_when(biomarker=="lab_alp" ~ "ALP",
                                          biomarker=="lab_alt" ~ "ALT",
                                          biomarker=="lab_ast" ~ "AST",
                                          biomarker=="lab_bilirubin" ~ "T-bil",
                                          biomarker=="lab_ck" ~ "CK",
                                          biomarker=="lab_creatinine" ~ "Cre",
                                          biomarker=="lab_crp" ~ "CRP",
                                          biomarker=="lab_d_dimer" ~ "D-dimer",
                                          biomarker=="lab_fibrinogen" ~ "Fibrinogen",
                                          biomarker=="lab_ggt" ~ "GGT",
                                          biomarker=="lab_il_6" ~ "IL-6",
                                          biomarker=="lab_ldh" ~ "LDH",
                                          biomarker=="lab_lymphocytes" ~ "Lymph",
                                          biomarker=="lab_monocytes" ~ "Mono",
                                          biomarker=="lab_neutrophils" ~ "Neut",
                                          biomarker=="lab_platelets" ~ "Plt",
                                          biomarker=="lab_procalcitonin" ~ "PCT",
                                          biomarker=="lab_serum_ferritin" ~ "Ferritin",
                                          biomarker=="lab_trop_t" ~ "Trop-T",
                                          biomarker=="lab_wbc" ~ "WBC"
                                          ))

OUT1 <- OUT %>% filter(pvalue.y < 0.05)

png("plots/Figure1C.png", width=900, height = 500)
p <- ggplot(OUT, aes(x=beta.y, xmin=LL.y, xmax=UL.y, y=beta.x,ymin=LL.x, ymax=UL.x)) + 
  theme_classic() +
  scale_alpha_discrete(range = c(1, 0.3),guide = 'none') + 
  geom_pointrange(aes(alpha=sig),lwd=0.5,size=5) + scale_color_brewer(palette = "Set1") +
  geom_errorbarh(aes(height = 0,alpha=sig),lwd=0.5) + 
  geom_smooth(method='lm',se=FALSE,color="#e4007f",alpha=0.3) + labs(alpha="p<0.05") +
  geom_label_repel(aes(label=biomarker), size=5,color="#0068b7") + 
  xlab("Coefficients for the associations\nbetween biomarkers and death and severe respiratory failure") + 
  ylab("Coefficients for the associations\nbetween rs10490770 and biomarkers") +
  theme(plot.title=element_text(size=15),
        axis.text.x=element_text(size=15,face="bold"),
        axis.text.y=element_text(size=15,face="bold"),
        axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text(size=15), 
        legend.title = element_text(size=15,face="bold"))
print(p)
dev.off()

