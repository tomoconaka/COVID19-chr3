attach("all_cause_EUR.obj")
library(survminer)
ggsurv <- ggsurvplot(fit,conf.int=TRUE,pval=FALSE, risk.table=TRUE,fun = "event", xlim = c(0,30),
           ylim = c(0,0.15), legend.title = "rs10490770 risk allele", 
           break.time.by = 5,xlab = "Time in days",palette = c("#E7B800", "#2E9FDF"),
           legend.labs = 
             c(paste0("Non-carrier (N=",dim(data_EUR[data_EUR$snp==0,])[1],")"), paste0("Carrier (N =",dim(data_EUR[data_EUR$snp==1,])[1],")")),
           font.x = 20,risk.table.fontsize = 10,font.tickslab = c(20),
           font.y = 20,font.legend = list(size = 20, color = "black"),
           tables.theme = clean_theme())

ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate("text", x=10,  y=0.12, size=15, label=paste0("HR: 1.4 [1.2; 1.6]"))+
  ggplot2::theme(legend.title = element_text(size = 20, color = "black"),
                )

      
ggsurv$table <- ggpubr::ggpar(ggsurv$table, 
                              font.title = list(size = 20, color = "black", face = "bold"),
                              font.main = list(size=20),
                              font.ytickslab = list(size=20))
  
ggsurv
png("all_cause_mortality_EUR.png", width=900, height=700)
ggsurv
dev.off()

attach("covid_cause_EUR.obj")
ggsurv <- ggcompetingrisks(ci_fit,xlim = c(0,30),palette = "Dark2",
                           conf.int=FALSE,pval=FALSE, risk.table=TRUE,fun = "event", xlim = c(0,30),
                     ylim = c(0,0.15), legend.title = "rs10490770", 
                     break.time.by = 5,xlab = "Time in days",palette = c("#E7B800", "#2E9FDF"),
                     legend.labs = 
                       c(paste0("TT (N =",dim(data_EUR[data_EUR$snp==0,])[1],")"), paste0("TC or CC (N =",dim(data_EUR[data_EUR$snp==1,])[1],")")),
                     font.x = 20,font.tickslab = c(20),
                     font.y = 20) 
ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate("text", x=5,  y=0.12, size=10, label=paste0("HR: ",round(c(res$conf.int[1,c(1,3,4)], res$coef[1,5]) ,2)[1]," [",round(c(res$conf.int[1,c(1,3,4)], res$coef[1,5]) ,2)[2],", ",round(c(res$conf.int[1,c(1,3,4)], res$coef[1,5]) ,2)[3],"]"))

fit5$states <- c("Alive", "Non-COVID-19 death", "COVID-19 death")
png("covid_cause_mortality_EUR.png", width=800, height=700)
ggcompetingrisks(fit5, palette = "Dark2",xlim=c(0,30),ylim=c(0,0.15),
                 legend = "top",xlab = "Time in days",legend.title = "",
                 font.x = 20,main="",font.legend = list(size = 20, color = "black"),font.tickslab = c(20),
                 font.y = 20)
dev.off()

png("all_cause_mortality_EUR.png", width=800, height=700)
ggsurv
dev.off()
