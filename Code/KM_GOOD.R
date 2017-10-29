setwd('/Users/NeilWu/Github/Research_EGFRmut_survival/KMplot/')
library(Hmisc)
TCGA_data_EGFR_sub_S1PR4<-TCGA_wt
GSE_S1PR4<-GSE31210_wt
GSE2_S1PR4<-GSE13213_wt

cox.predict <- predict(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_data_EGFR_sub_S1PR4[,'S1PR4']), data = TCGA_data_EGFR_sub_S1PR4),
                       TCGA_data_EGFR_sub_S1PR4)
TCGA_data_EGFR_sub_S1PR4$risk<-cut2(cox.predict,g=2)
levels(TCGA_data_EGFR_sub_S1PR4$risk)<-c('High','Low')
S1PR4.TCGA.surv.fit<-survfit(Surv(time = Date/30, event = Relapse) ~ risk,data = TCGA_data_EGFR_sub_S1PR4)
logr<-survdiff(Surv(time = Date/30, event = Relapse) ~ risk,data = TCGA_data_EGFR_sub_S1PR4)
pdf('S1PR4(TCGA-).pdf',width = 6,height = 5)
plot(S1PR4.TCGA.surv.fit, col=c("blue", "red"), lty=1:2, xlab='Time (month)',
     main = 'S1PR4 in EGFR- (TCGA)',xlim = c(0,60))
legend("bottomleft",col=c('blue','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('logrank p-value= %.3e',1 - pchisq(logr$chisq, 1)))
dev.off()



#colnames(GSE_S1PR4)[6:length(colnames(GSE_S1PR4))]<-exp_header[6:length(colnames(GSE_S1PR4))]
d<-which(colnames(GSE_S1PR4)=='S1PR4')#[which.min(cox.p.relapse[which(colnames(GSE_S1PR4)=='S1PR4')])]
cox.predict <- predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE_S1PR4[,d[1]]), data = GSE_S1PR4),
                       GSE_S1PR4)
GSE_S1PR4$risk<-cut2(cox.predict,g=2)
levels(GSE_S1PR4$risk)<-c('Low','High')
S1PR4.GSE.surv.fit<-survfit(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE_S1PR4)
logr<-survdiff(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE_S1PR4)
pdf('S1PR4(GSE31210-).pdf',width = 6,height = 5)
plot(S1PR4.GSE.surv.fit, col=c("blue", "red"), lty=1:2, xlab='Time (month)',
     main = 'S1PR4 in EGFR- (GSE31210)',xlim = c(0,60))
legend("bottomleft",col=c('blue','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('logrank p-value= %.3e',1 - pchisq(logr$chisq, 1)))
dev.off()




#colnames(GSE2_S1PR4)[6:length(colnames(GSE2_S1PR4))]<-exp_header[6:length(colnames(GSE2_S1PR4))]
d<-which(colnames(GSE2_S1PR4)=='S1PR4')#[which.min(cox.p.relapse[which(colnames(GSE2_S1PR4)=='S1PR4')])]
cox.predict <- predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE2_S1PR4[,d]), data = GSE2_S1PR4),
                       GSE2_S1PR4)
GSE2_S1PR4$risk<-cut2(cox.predict,g=2)
levels(GSE2_S1PR4$risk)<-c('Low','High')
S1PR4.GSE2.surv.fit<-survfit(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE2_S1PR4)
logr<-survdiff(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE2_S1PR4)
pdf('S1PR4(GSE13213-).pdf',width = 6,height = 5)
plot(S1PR4.GSE2.surv.fit, col=c("blue", "red"), lty=1:2, xlab='Time (month)',
     main = 'S1PR4 in EGFR- (GSE13213)',xlim = c(0,60))
legend("bottomleft",col=c('blue','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('logrank p-value= %.3e',1 - pchisq(logr$chisq, 1)))
dev.off()


#############################

TCGA_data_EGFR_sub_S1PR4<-TCGA_mut
GSE_S1PR4<-GSE31210_mut
GSE2_S1PR4<-GSE13213_mut

cox.predict <- predict(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_data_EGFR_sub_S1PR4[,'S1PR4']), data = TCGA_data_EGFR_sub_S1PR4),
                       TCGA_data_EGFR_sub_S1PR4)
TCGA_data_EGFR_sub_S1PR4$risk<-cut2(cox.predict,g=2)
levels(TCGA_data_EGFR_sub_S1PR4$risk)<-c('High','Low')
S1PR4.TCGA.surv.fit<-survfit(Surv(time = Date/30, event = Relapse) ~ risk,data = TCGA_data_EGFR_sub_S1PR4)
logr<-survdiff(Surv(time = Date/30, event = Relapse) ~ risk,data = TCGA_data_EGFR_sub_S1PR4)
pdf('S1PR4(TCGA+).pdf',width = 6,height = 5)
plot(S1PR4.TCGA.surv.fit, col=c("blue", "red"), lty=1:2, xlab='Time (month)',
     main = 'S1PR4 in EGFR+ (TCGA)',xlim = c(0,60))
legend("bottomleft",col=c('blue','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('logrank p-value= %.3e',1 - pchisq(logr$chisq, 1)))
dev.off()



#colnames(GSE_S1PR4)[6:length(colnames(GSE_S1PR4))]<-exp_header[6:length(colnames(GSE_S1PR4))]
d<-which(colnames(GSE_S1PR4)=='S1PR4')#[which.min(cox.p.relapse[which(colnames(GSE_S1PR4)=='S1PR4')])]
cox.predict <- predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE_S1PR4[,d[1]]), data = GSE_S1PR4),
                       GSE_S1PR4)
GSE_S1PR4$risk<-cut2(cox.predict,g=2)
levels(GSE_S1PR4$risk)<-c('Low','High')
S1PR4.GSE.surv.fit<-survfit(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE_S1PR4)
logr<-survdiff(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE_S1PR4)
pdf('S1PR4(GSE31210+).pdf',width = 6,height = 5)
plot(S1PR4.GSE.surv.fit, col=c("blue", "red"), lty=1:2, xlab='Time (month)',
     main = 'S1PR4 in EGFR+ (GSE31210)',xlim = c(0,60))
legend("bottomleft",col=c('blue','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('logrank p-value= %.3e',1 - pchisq(logr$chisq, 1)))
dev.off()




#colnames(GSE2_S1PR4)[6:length(colnames(GSE2_S1PR4))]<-exp_header[6:length(colnames(GSE2_S1PR4))]
d<-which(colnames(GSE2_S1PR4)=='S1PR4')#[which.min(cox.p.relapse[which(colnames(GSE2_S1PR4)=='S1PR4')])]
cox.predict <- predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE2_S1PR4[,d]), data = GSE2_S1PR4),
                       GSE2_S1PR4)
GSE2_S1PR4$risk<-cut2(cox.predict,g=2)
levels(GSE2_S1PR4$risk)<-c('Low','High')
S1PR4.GSE2.surv.fit<-survfit(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE2_S1PR4)
logr<-survdiff(Surv(time = relapse_day/30, event = relapse_status) ~ risk,data = GSE2_S1PR4)
pdf('S1PR4(GSE13213+).pdf',width = 6,height = 5)
plot(S1PR4.GSE2.surv.fit, col=c("blue", "red"), lty=1:2, xlab='Time (month)',
     main = 'S1PR4 in EGFR+ (GSE13213)',xlim = c(0,60))
legend("bottomleft",col=c('blue','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('logrank p-value= %.3e',1 - pchisq(logr$chisq, 1)))
dev.off()

rm(list = grep('GSE13213_mut|GSE13213_wt|TCGA_wt|TCGA_mut|GSE31210_wt|GSE31210_mut|Data_table_val',
               ls(),invert = T,value = T))