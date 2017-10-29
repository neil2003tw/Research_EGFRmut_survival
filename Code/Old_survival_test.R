temp<-read.table('../GSE/GSE31210_series_matrix.txt',skip = 25,fill=T,stringsAsFactors = F)
idx_EGFR_mut<-grep('EGFR mutation',temp[17,])
idx_EGFR_free<-which(!grepl('EGFR mutation',temp[17,]))
temp_subER<-temp[,c(1,idx_EGFR_free)]

relapse_status<-!grepl('not',temp_subER[20,])
relapse_day<-unname(as.numeric(sapply(temp_subER[21,], function(x) strsplit(x,'days before relapse/censor: ')[[1]][2])))
event_dead<-grepl('dead',temp_subER[23,])
event_day<-unname(as.numeric(sapply(temp_subER[24,], function(x) strsplit(x,'days before death/censor: ')[[1]][2])))
sample_id<-as.character(temp_subER[49,])
relapse_day[is.na(relapse_day)]<-event_day[is.na(relapse_day)]

clinical_data<-data.frame(sample_id,relapse_status,relapse_day,event_dead,event_day)

expression_data<-temp_subER[49:(dim(temp_subER)[1]-1),]

U133plus2<-read.csv('/Users/NeilWu/Github/Research_lncRNA-AML/HG-U133_Plus_2.na34.annot.csv',stringsAsFactors=F)
expression_data[,1]<-as.character(U133plus2$Gene.Symbol[match(expression_data[,1],as.character(U133plus2$Probe.Set.ID))])

exp_temp4<-data.frame(t(expression_data),stringsAsFactors = F)
exp_temp4<-cbind(clinical_data,exp_temp4[,-1])
exp_header<-as.vector(t(exp_temp[1,]))
exp_temp4<-exp_temp4[-1,]
for(i in 6:dim(exp_temp)[2]){exp_temp[,i]<-as.numeric(exp_temp[,i])}

library(survival)

pb<-txtProgressBar(min=6,max=dim(exp_temp)[2],style=3)
cox.p <- vector()
for( i in 6:dim(exp_temp)[2]){
  cox.p[i] <- summary(coxph(Surv(time = event_day, event = event_dead) ~ exp_temp[,i], data = exp_temp))$coefficient[5]
  setTxtProgressBar(pb,i)
}

cox_sig<-which(cox.p<0.05)
cox_sig<-cox_sig[-grep('---',exp_header[which(cox.p<0.05)])]
temp<-table(exp_header[which(cox.p<0.00001)])

cox_sig_abovehalf<-cox_sig[c(2,6,12,13,15,16,17,18,19,22,23,26,27,28,29,31)]



pb<-txtProgressBar(min=6,max=dim(exp_temp)[2],style=3)
cox.p.relapse <- vector()
for( i in 6:dim(exp_temp)[2]){
  cox.p.relapse[i] <- summary(coxph(Surv(time = relapse_day, event = relapse_status) ~ exp_temp[,i], data = exp_temp))$coefficient[5]
  setTxtProgressBar(pb,i)
}


idx_mut<-which(cox.p.relapse<0.05)
idx_free<-which(cox.p.relapse_free<0.05)
cox_sig<-idx_mut[!idx_mut %in% idx_free]
cox_sig<-cox_sig[-grep('---',exp_header[cox_sig])]



TCGA_p<-read.table('nationwidechildrens.org_clinical_patient_luad.txt',header = T,sep='\t',stringsAsFactors = F)
TCGA_p<-TCGA_p[-c(1,2),]
TCGA_f<-read.table('nationwidechildrens.org_clinical_follow_up_v1.0_luad.txt',header=T,sep='\t',stringsAsFactors = F)
TCGA_f<-TCGA_f[-c(1,2),]
TCGA_m<-read.table('mutation.txt',header=T,sep='\t',stringsAsFactors = F)
TCGA_srdf<-read.table('unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.1.15.0.sdrf.txt',header = T,sep='\t',stringsAsFactors = F)

#even_death<-grepl('Dead',TCGA_f$vital_status)
event_date<-c()
even_relapse<-!grepl('Not',TCGA_f$new_tumor_event_dx_days_to)
for(i in 1:length(even_relapse)){
  if(even_relapse[i]){
    d<-as.numeric(TCGA_f$new_tumor_event_dx_days_to[i])
  }else{
    d<-as.numeric(TCGA_f$last_contact_days_to[i])
    if(is.na(d)){
      d<-as.numeric(TCGA_f$death_days_to[i])
    }
  }
  event_date[i]<-d
}


TCGA_clinical<-data.frame(BCR_uuid=TCGA_f$bcr_patient_uuid,
                          BCR_barcode=TCGA_p$bcr_patient_barcode[match(TCGA_f$bcr_patient_uuid,TCGA_p$bcr_patient_uuid)],
                          #Event_dead=even_death,
                          Event_relapse=even_relapse,
                          Event_day=event_date,
                          Stage=TCGA_p$ajcc_pathologic_tumor_stage[match(TCGA_f$bcr_patient_uuid,TCGA_p$bcr_patient_uuid)])
TCGA_clinical$mutate<-TCGA_m$Hugo_Symbol[match(TCGA_clinical$BCR_barcode, substr(TCGA_m$Tumor_Sample_Barcode,1,12))]
TCGA_clinical$mutate[is.na(TCGA_clinical$mutate)]<-'No'

TCGA_srdf<-TCGA_srdf[grepl('genes.normalized_results',TCGA_srdf$Derived.Data.File),]
TCGA_clinical$Derived_file<-TCGA_srdf$Derived.Data.File[match(TCGA_clinical$BCR_barcode,substr(TCGA_srdf$Comment.TCGA.Barcode,1,12))]

TCGA_data<-read.csv('TCGA_AD_seq_summary_symbol_0626.csv')
TCGA_clinical$Derived_file_idx<-match(substr(TCGA_clinical$Derived_file,1,16),substr(colnames(TCGA_data),1,16))

duplicated(TCGA_clinical$Derived_file_idx)

TCGA_data_EGFR<-TCGA_data[,c(1,TCGA_clinical$Derived_file_idx[!is.na(TCGA_clinical$Derived_file_idx)])]
TCGA_data_EGFR<-data.frame(t(TCGA_data_EGFR),stringsAsFactors = F)

TCGA_sig_idx<-match(exp_header[cox_sig],TCGA_data_EGFR[1,])
TCGA_sig_idx<-c(TCGA_sig_idx[!is.na(TCGA_sig_idx)])
TCGA_data_EGFR_sub<-TCGA_data_EGFR[,TCGA_sig_idx[!duplicated(TCGA_sig_idx)]]
colnames(TCGA_data_EGFR_sub)<-TCGA_data_EGFR_sub[1,]
TCGA_data_EGFR_sub<-TCGA_data_EGFR_sub[-1,]
TCGA_data_EGFR_sub<-TCGA_data_EGFR_sub[grepl('results$',row.names(TCGA_data_EGFR_sub)),]
TCGA_data_EGFR_sub$Derived_file<-substr(row.names(TCGA_data_EGFR_sub),1,16)

TCGA_clinical_sub<-data.frame(Relapse=TCGA_clinical$Event_relapse,
                              #Dead=TCGA_clinical$Event_dead,
                              Date=TCGA_clinical$Event_day,
                              Mutate=TCGA_clinical$mutate,
                              Stage=TCGA_clinical$Stage,
                              Derived_file=substr(TCGA_clinical$Derived_file,1,16))

TCGA_data_EGFR_sub_w_clinic<-merge(TCGA_data_EGFR_sub,TCGA_clinical_sub,by='Derived_file')
TCGA_data_EGFR_sub_w_clinic$Date<-as.numeric(as.character(TCGA_data_EGFR_sub_w_clinic$Date))
TCGA_data_EGFR_sub_w_clinic<-TCGA_data_EGFR_sub_w_clinic[!is.na(TCGA_data_EGFR_sub_w_clinic$Date),]
TCGA_data_EGFR_sub_w_clinic<-TCGA_data_EGFR_sub_w_clinic[as.character(TCGA_data_EGFR_sub_w_clinic$Mutate)=='EGFR',]
TCGA_data_EGFR_sub_w_clinic<-TCGA_data_EGFR_sub_w_clinic[!grepl('Stage III',as.character(TCGA_data_EGFR_sub_w_clinic$Stage)),]
TCGA_data_EGFR_sub_w_clinic$Date<-as.numeric(TCGA_data_EGFR_sub_w_clinic$Date)

cox.validate<-c()
pb<-txtProgressBar(min=2,max=dim(TCGA_data_EGFR_sub_w_clinic)[2]-3,style=3)
for(i in 2:(dim(TCGA_data_EGFR_sub_w_clinic)[2]-3)){
  setTxtProgressBar(pb,i)
  cox.validate[i] <- summary(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_data_EGFR_sub_w_clinic[,i]), data = TCGA_data_EGFR_sub_w_clinic))$coefficient[5]
}


gene_list005<-exp_header[cox.p.relapse<0.05][exp_header[cox.p.relapse<0.05] %in% colnames(TCGA_data_EGFR_sub_w_clinic)[cox.validate<0.05]]
gene_list005<-gene_list005[!duplicated(gene_list005)][-1]
gene_list001<-exp_header[cox.p.relapse<0.01][exp_header[cox.p.relapse<0.01] %in% colnames(TCGA_data_EGFR_sub_w_clinic)[cox.validate<0.01]]
gene_list001<-gene_list001[!duplicated(gene_list001)][-1]
gene_list0005<-exp_header[cox.p.relapse<0.005][exp_header[cox.p.relapse<0.005] %in% colnames(TCGA_data_EGFR_sub_w_clinic)[cox.validate<0.005]]
gene_list0005<-gene_list0005[!duplicated(gene_list0005)][-1]
gene_list0001<-exp_header[cox.p.relapse<0.001][exp_header[cox.p.relapse<0.001] %in% colnames(TCGA_data_EGFR_sub_w_clinic)[cox.validate<0.001]]
gene_list0001<-gene_list0001[!duplicated(gene_list0001)][-1]

colnames(TCGA_data_EGFR_sub_w_clinic)[cox.validate<0.001][colnames(TCGA_data_EGFR_sub_w_clinic)[cox.validate<0.001] %in% exp_header[cox.p.relapse<0.01]]
exp_header[cox.p.relapse<0.05]



################KMplotter

TCGA_data_EGFR_sub_CCT6A<-TCGA_data_EGFR_sub_w_clinic
cox.predict <- predict(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_data_EGFR_sub_CCT6A[,'CCT6A']), data = TCGA_data_EGFR_sub_CCT6A),
                       TCGA_data_EGFR_sub_CCT6A)
TCGA_data_EGFR_sub_CCT6A$risk<-factor(cox.predict>0)
levels(TCGA_data_EGFR_sub_CCT6A$risk)<-c('High','Low')
CCT6A.TCGA.surv.fit<-survfit(Surv(time = Date/30, event = Relapse) ~ risk,data = TCGA_data_EGFR_sub_CCT6A)
pdf('CCT6A(TCGA).pdf',width = 6,height = 5)
plot(CCT6A.TCGA.surv.fit, col=c("green", "red"), lty=1:2, xlab='Time (month)',main = 'CCT6A (TCGA)')
legend("bottomleft",col=c('green','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('cox p-value= %.3e',cox.validate[which(colnames(TCGA_data_EGFR_sub_w_clinic)=='CCT6A')]))
dev.off()


GSE_CCT6A<-exp_temp
colnames(GSE_CCT6A)[6:length(colnames(GSE_CCT6A))]<-exp_header[6:length(colnames(GSE_CCT6A))]
d<-which(colnames(GSE_CCT6A)=='CCT6A')[which.min(cox.p.relapse[which(colnames(GSE_CCT6A)=='CCT6A')])]
cox.predict <- predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE_CCT6A[,d]), data = GSE_CCT6A),
                       GSE_CCT6A)
GSE_CCT6A$risk<-factor(cox.predict>0)
levels(GSE_CCT6A$risk)<-c('High','Low')
CCT6A.GSE.surv.fit<-survfit(Surv(time = event_day/30, event = event_dead) ~ risk,data = GSE_CCT6A)
pdf('CCT6A(GSE).pdf',width = 6,height = 5)
plot(CCT6A.GSE.surv.fit, col=c("green", "red"), lty=1:2, xlab='Time (month)',main = 'CCT6A (GSE31210)')
legend("bottomleft",col=c('green','red'),lty=1:2,legend=c("Low risk","High risk"))
legend('bottomright',sprintf('cox p-value= %.3e',cox.p.relapse[d]))
dev.off()

rm(list=c('TCGA_data_EGFR_sub_CCT6A','GSE_CCT6A'))

########

write.table(gene_list0001,'Sig_gene_list_p0001.txt',quote = F,col.names = F,row.names = F,sep='\n')
write.table(gene_list0005,'Sig_gene_list_p0005.txt',quote = F,col.names = F,row.names = F,sep='\n')
write.table(gene_list001,'Sig_gene_list_p001.txt',quote = F,col.names = F,row.names = F,sep='\n')
write.table(gene_list005,'Sig_gene_list_p005.txt',quote = F,col.names = F,row.names = F,sep='\n')

#######

Data_table<-data.frame(Gene=gene_list005,GSE_beta=NA,GSE_p=NA,TCGA_beta=NA,TCGA_p=NA,stringsAsFactors = F)

for(i in 1:dim(Data_table)[1]){
  ### GSE_cox
  idx<-which(exp_header==Data_table$Gene[i])[which.min(cox.p.relapse[which(exp_header==Data_table$Gene[i])])]
  GSE_sub<-exp_temp[-1,c(1:5,idx)]
  cox_summary_G<-summary(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE_sub[,6]), data = GSE_sub))
  Data_table$GSE_beta[i]<-cox_summary_G$coefficients[1]
  Data_table$GSE_p[i]<-cox_summary_G$coefficients[5]
  
  ### TCGA_cox
  idx<-which(colnames(TCGA_data_EGFR_sub_w_clinic)==Data_table$Gene[i])
  TCGA_sub<-TCGA_data_EGFR_sub_w_clinic[,c(6505:6508,idx)]
  cox_summary_T<-summary(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_sub[,5]), data = TCGA_sub))
  Data_table$TCGA_beta[i]<-cox_summary_T$coefficients[1]
  Data_table$TCGA_p[i]<-cox_summary_T$coefficients[5]
}
  
  
Data_table<-Data_table[order(Data_table$GSE_p),]  
write.table()