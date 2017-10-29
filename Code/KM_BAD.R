TCGA_data_EGFR_sub_GNPNAT1<-TCGA_data_w_clinic_wt
cox_ph<-coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_data_EGFR_sub_GNPNAT1[,'GNPNAT1']),data = TCGA_data_EGFR_sub_GNPNAT1)
cox.predict <- predict(cox_ph, TCGA_data_EGFR_sub_GNPNAT1)
TCGA_data_EGFR_sub_GNPNAT1$risk<-factor(cox.predict>0)
levels(TCGA_data_EGFR_sub_GNPNAT1$risk)<-c('High','Low')
GNPNAT1.TCGA.surv.fit<-survfit(Surv(time = Date/30, event = Relapse) ~ risk,data = TCGA_data_EGFR_sub_GNPNAT1)
logr<-survdiff(Surv(time = Date/30, event = Relapse) ~ risk,data = TCGA_data_EGFR_sub_GNPNAT1)
pdf('GNPNAT1\'(TCGA).pdf',width = 6,height = 5)
plot(GNPNAT1.TCGA.surv.fit, col=c("blue", "red"), lty=1:2, xlab='Time (month)',main = 'GNPNAT1 in EGFR- (TCGA)')
legend("bottomleft",col=c('blue','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('logrank p-value= %.3e',1 - pchisq(logr$chisq, 1)))
dev.off()


GSE_GNPNAT1<-exp_temp2
colnames(GSE_GNPNAT1)[6:length(colnames(GSE_GNPNAT1))]<-exp_header[6:length(colnames(GSE_GNPNAT1))]
d<-which(colnames(GSE_GNPNAT1)=='GNPNAT1')[which.min(cox.p.relapse[which(colnames(GSE_GNPNAT1)=='GNPNAT1')])]
cox_ph<-coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE_GNPNAT1[,d]), data = GSE_GNPNAT1)
cox.predict <- predict(cox_ph,GSE_GNPNAT1)
GSE_GNPNAT1$risk<-factor(cox.predict>0)
levels(GSE_GNPNAT1$risk)<-c('High','Low')
GNPNAT1.GSE.surv.fit<-survfit(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE_GNPNAT1)
logr<-survdiff(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE_GNPNAT1)
pdf('GNPNAT1\'(GSE).pdf',width = 6,height = 5)
plot(GNPNAT1.GSE.surv.fit, col=c("blue", "red"), lty=1:2, xlab='Time (month)',main = 'GNPNAT1 in EGFR- (GSE31210)')
legend("bottomleft",col=c('blue','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('logrank p-value= %.3e',1 - pchisq(logr$chisq, 1)))
dev.off()

rm(list=c('TCGA_data_EGFR_sub_GNPNAT1','GSE_GNPNAT1'))
