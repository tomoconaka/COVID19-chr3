setwd("/home/tomoko/02.associations/")
tmp <- readRDS("~/01.data.QC/curated_clinical/complication_dat_final.rds")

unique(tmp$study)
tmp <- tmp %>% mutate(study = case_when(study == unique(tmp$study)[2] ~ "BelCovid_1",
                                        study == unique(tmp$study)[3] ~ "BelCovid_2",
                                        study == unique(tmp$study)[4] ~ "BRACOVID",
                                        study == unique(tmp$study)[5] ~ "BQC19",
                                        study == unique(tmp$study)[6] ~ "BosCo",
                                        study == unique(tmp$study)[7] ~ "COMRI",
                                        study == unique(tmp$study)[8] ~ "Hostage_4",
                                        study == unique(tmp$study)[9] ~ "GEN-COVID",
                                        study == unique(tmp$study)[10] ~ "FoGS",
                                        study == unique(tmp$study)[11] ~ "NorCov2",
                                        study == unique(tmp$study)[12] ~ "SPGRX",
                                        study == unique(tmp$study)[13] ~ "Hostage_1",
                                        study == unique(tmp$study)[14] ~ "Hostage_2",
                                        study == unique(tmp$study)[15] ~ "Hostage_3",
                                        study == unique(tmp$study)[16] ~ "INMUNGEN-CoV2",
                                        study == unique(tmp$study)[17] ~ "SweCovid",
                                        TRUE ~ "UKB"
                                        ))

p1 <- ggplot(tmp, aes(x=PC1, y=PC2, col=pop)) + geom_point() + theme_classic() +
  scale_color_brewer(palette="Set1") + 
  theme(axis.title =element_text(size=20)) + labs(color="Population") + guides(color=FALSE)
p2 <- ggplot(tmp, aes(x=PC2, y=PC3, col=pop)) + geom_point() + theme_classic() +
  scale_color_brewer(palette="Set1") + 
  theme(axis.title =element_text(size=20)) + labs(color="Population") + guides(color=FALSE)
p3 <- ggplot(tmp, aes(x=PC3, y=PC4, col=pop)) + geom_point() + theme_classic() +
  scale_color_brewer(palette="Set1") + 
  theme(axis.title =element_text(size=20),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10)) + labs(color="Population") 
p4 <- ggplot(tmp, aes(x=PC1, y=PC2, col=study)) + geom_point() + theme_classic() +
  theme(axis.title =element_text(size=20)) + labs(color="Study") + guides(color=FALSE)
p5 <- ggplot(tmp, aes(x=PC2, y=PC3, col=study)) + geom_point() + theme_classic() +
  theme(axis.title =element_text(size=20)) + labs(color="Study") + guides(color=FALSE)
p6 <- ggplot(tmp, aes(x=PC3, y=PC4, col=study)) + geom_point() + theme_classic() + 
  theme(axis.title =element_text(size=20),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10)) + labs(color="Study") 

png("plots/SuppleFig1.png", width=1200, height=800)
grid.arrange(arrangeGrob(p1,p2,p4,p5, ncol=2, nrow=2), arrangeGrob(p3,p6, ncol=1, nrow=2), widths=c(3,2))
dev.off()
