temp<-read.table('/Users/NeilWu/Github/Research_EGFRmut_survival/GSE31210/GSE31210_series_matrix.txt',
                 skip = 25,fill=T,stringsAsFactors = F)
U133plus2<-read.csv('/Users/NeilWu/Github/Research_lncRNA-AML/HG-U133_Plus_2.na34.annot.csv',stringsAsFactors=F)
idx_EGFR_mut<-grep('EGFR mutation +',temp[17,])
idx_EGFR_wt<-which(grepl('gene alteration status:',temp[17,]) & !grepl('EGFR mutation +',temp[17,]))
temp_sub_wt<-temp[,c(1,idx_EGFR_wt)]
temp_sub_mut<-temp[,c(1,idx_EGFR_mut)]

z_trans<-function(l){
  l<-as.numeric(l)
  x<-l[!is.na(l)]
  y<-(x-mean(x))/sd(x)
  l[!is.na(l)]<-y
  return(l)
}

############### ER - clinical&merge
relapse_status<-c()
for(i in 1:dim(temp_sub_wt)[2]){
  relapse_status_idx<-grep('relapse: ',temp_sub_wt[1:30,i])
  relapse_status[i]<-!grepl('not',temp_sub_wt[relapse_status_idx,i])
}
relapse_day<-c()
for(i in 1:dim(temp_sub_wt)[2]){
  relapse_day_idx<-grep('days before relapse/censor:',temp_sub_wt[1:35,i])
  relapse_day[i]<-unname(as.numeric(sapply(temp_sub_wt[relapse_day_idx,i], function(x) strsplit(x,'days before relapse/censor: ')[[1]][2])))
}
sample_id<-as.character(temp_sub_wt[49,])
clinical_data_wt<-data.frame(sample_id,relapse_status,relapse_day)
expression_data<-temp_sub_wt[49:(dim(temp_sub_wt)[1]-1),] ## row49: GSM ID

### probe-gene replacement
matched_probe_raw<-as.character(U133plus2$Gene.Symbol[match(expression_data[,1],as.character(U133plus2$Probe.Set.ID))])
matched_probe_done<-sapply(matched_probe_raw,function(x) strsplit(x,' /// ')[[1]][1],USE.NAMES = F)
expression_data[,1]<-matched_probe_done
###

exp_temp_EGFR_wt<-data.frame(t(expression_data),stringsAsFactors = F)
exp_temp_EGFR_wt<-cbind(clinical_data_wt,exp_temp_EGFR_wt[,-1])
exp_header<-as.vector(t(exp_temp_EGFR_wt[1,]))
exp_temp_EGFR_wt<-exp_temp_EGFR_wt[-c(1,2),]
###




############### ER + clinical&merge
relapse_status<-c()
for(i in 1:dim(temp_sub_mut)[2]){
  relapse_status_idx<-grep('relapse: ',temp_sub_mut[1:30,i])
  relapse_status[i]<-!grepl('not',temp_sub_mut[relapse_status_idx,i])
}
relapse_day<-c()
for(i in 1:dim(temp_sub_mut)[2]){
  relapse_day_idx<-grep('days before relapse/censor:',temp_sub_mut[1:35,i])
  relapse_day[i]<-as.numeric(sapply(temp_sub_mut[relapse_day_idx,i], function(x) strsplit(x,'days before relapse/censor: ')[[1]][2],USE.NAMES = F))
}
sample_id<-as.character(temp_sub_mut[49,])
clinical_data_mut<-data.frame(sample_id,relapse_status,relapse_day)
expression_data<-temp_sub_mut[49:(dim(temp_sub_mut)[1]-1),]

### probe-gene replacement
matched_probe_raw<-as.character(U133plus2$Gene.Symbol[match(expression_data[,1],as.character(U133plus2$Probe.Set.ID))])
matched_probe_done<-matched_probe_raw ## Test for multiple gene probeset 15/10/07
#matched_probe_done<-sapply(matched_probe_raw,function(x) strsplit(x,' /// ')[[1]][1],USE.NAMES = F)
expression_data[,1]<-matched_probe_done
###

exp_temp_EGFR_mut<-data.frame(t(expression_data),stringsAsFactors = F)
exp_temp_EGFR_mut<-cbind(clinical_data_mut,exp_temp_EGFR_mut[,-1])
exp_header<-as.vector(t(exp_temp_EGFR_mut[1,]))
exp_temp_EGFR_mut<-exp_temp_EGFR_mut[-1,]
##################


################## Header add
colnames(exp_temp_EGFR_wt)[-c(1:3)]<-exp_header[-c(1:3)]
colnames(exp_temp_EGFR_mut)[-c(1:3)]<-exp_header[-c(1:3)]
##################


################# Log transform & Z-score
for (i in 4:dim(exp_temp_EGFR_wt)[2]){
  logscale_temp<-log(as.numeric(exp_temp_EGFR_wt[,i]),2)
  ztrans_temp<-(logscale_temp-mean(logscale_temp))/sd(logscale_temp)
  exp_temp_EGFR_wt[,i]<-ztrans_temp
}

for (i in 4:dim(exp_temp_EGFR_mut)[2]){
  logscale_temp<-log(as.numeric(exp_temp_EGFR_mut[,i]),2)
  ztrans_temp<-(logscale_temp-mean(logscale_temp))/sd(logscale_temp)
  exp_temp_EGFR_mut[,i]<-ztrans_temp
}
#################





################# EGFR wt cox loop
library(survival)
cox.p.relapse_wt<-c()
pb<-txtProgressBar(min=3,max=dim(exp_temp_EGFR_wt)[2],style=3)
for( i in 4:dim(exp_temp_EGFR_wt)[2]){
  cox.summary<-summary(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(exp_temp_EGFR_wt[,i]), data = exp_temp_EGFR_wt))
  cox.p.relapse_wt[i]<-cox.summary$coefficient[5]
  setTxtProgressBar(pb,i)
}
write.table(cox.p.relapse_wt,'../cox_p_relapse_EGFR_wt.txt')
#################


################# EGFR Mut cox loop
cox.p.relapse_Mut<-c()
pb<-txtProgressBar(min=3,max=dim(exp_temp_EGFR_mut)[2],style=3)
for( i in 4:dim(exp_temp_EGFR_mut)[2]){
  cox.summary<-summary(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(exp_temp_EGFR_mut[,i]), data = exp_temp_EGFR_mut))
  cox.p.relapse_Mut[i]<-cox.summary$coefficient[5]
  setTxtProgressBar(pb,i)
}
write.table(cox.p.relapse_Mut,'../cox_p_relapse_EGFR_mut.txt')
#################

cox.sig.relapse_wt<-which(cox.p.relapse_wt<0.05)
cox.sig.relapse_wt.s<-cox.sig.relapse_wt[-grep('---',exp_header[cox.sig.relapse_wt])]
gene.sig.relapse_wt.s<-exp_header[cox.sig.relapse_wt.s]
gene.sig.relapse_wt.dr<-gene.sig.relapse_wt.s[!duplicated(gene.sig.relapse_wt.s)]

cox.sig.relapse_Mut<-which(cox.p.relapse_Mut<0.05)
cox.sig.relapse_Mut.s<-cox.sig.relapse_Mut[-grep('---',exp_header[cox.sig.relapse_Mut])]
gene.sig.relapse_Mut.s<-exp_header[cox.sig.relapse_Mut.s]
gene.sig.relapse_Mut.dr<-gene.sig.relapse_Mut.s[!duplicated(gene.sig.relapse_Mut.s)]

idx_espg<-cox.sig.relapse_Mut.s[!cox.sig.relapse_Mut.s %in% cox.sig.relapse_wt.s]
gene.sig.espr2<-exp_header[idx_espg][!duplicated(exp_header[idx_espg])]
idx_espr2<-idx_espg[!duplicated(exp_header[idx_espg])]
gene.sig.espr3<-gene.sig.espr2[!gene.sig.espr2 %in% gene.sig.relapse_wt.dr]
idx_espr3<-idx_espr2[!gene.sig.espr2 %in% gene.sig.relapse_wt.dr]

#gene.sig.espr<-gene.sig.relapse_Mut.dr[!(gene.sig.relapse_Mut.dr %in% gene.sig.relapse_wt.dr)]


################## Data creation
Data_table<-data.frame(Gene_symbol=gene.sig.espr3,EGFRp_cox_beta=NA,EGFRp_cox_se=NA,EGFRp_cox_p=NA, 
                       EGFRn_cox_beta=NA,EGFRn_cox_se=NA,EGFRn_cox_p=NA,stringsAsFactors =F)
pb<-txtProgressBar(min=1,max=dim(Data_table)[1],style=3)
for( i in 1:dim(Data_table)[1]){
  now_gene<-Data_table$Gene_symbol[i]
  ###EGFRp
  cox.summary.p<-summary(coxph(Surv(time = relapse_day, event = relapse_status) ~ z_trans(exp_temp_EGFR_mut[,idx_espr3[i]]),
                               data = exp_temp_EGFR_mut))
  Data_table$EGFRp_cox_p[i]<-cox.summary.p$coefficient[5]
  Data_table$EGFRp_cox_beta[i]<-cox.summary.p$coefficient[1]
  Data_table$EGFRp_cox_se[i]<-cox.summary.p$coefficient[3]
  ###EGFRn
  cox.summary.n<-summary(coxph(Surv(time = relapse_day, event = relapse_status) ~ z_trans(exp_temp_EGFR_wt[,idx_espr3[i]]),
                               data = exp_temp_EGFR_wt))
  Data_table$EGFRn_cox_p[i]<-cox.summary.n$coefficient[5]
  Data_table$EGFRn_cox_beta[i]<-cox.summary.n$coefficient[1]
  Data_table$EGFRn_cox_se[i]<-cox.summary.n$coefficient[3]
  ###
  setTxtProgressBar(pb,i)
}