Data_table<-read.table('/Users/NeilWu/Github/Research_EGFRmut_survival/Data_0317.txt',sep='\t',header = T)

### GSE13213
validate<-read.table('/Users/NeilWu/Github/Research_EGFRmut_survival/GSE13213/GSE13213_series_matrix.txt',
                     skip = 36,fill=T,stringsAsFactors = F)
GPL6480_annot<-read.table('/Users/NeilWu/Github/Research_EGFRmut_survival/GSE13213/GPL6480-9577.txt',
                          skip=17,fill=T,header=T,sep='\t',stringsAsFactors = F)
idx_validate_EGFR_mut<-grep('EGFR status: Mut',validate[22,])
idx_validate_EGFR_wt<-grep('EGFR status: Wt',validate[22,])
validate_EGFR_mut<-validate[,c(1,idx_validate_EGFR_mut)]
validate_EGFR_wt<-validate[,c(1,idx_validate_EGFR_wt)]

#### ER +
relapse_status<-c()
for(i in 1:dim(validate_EGFR_mut)[2]){
  relapse_status_idx<-grep('Evidence of relapse: ',validate_EGFR_mut[1:30,i])
  relapse_status[i]<-grepl('Evidence of relapse: Y',validate_EGFR_mut[relapse_status_idx,i])
}

relapse_day<-c()
for(i in 2:dim(validate_EGFR_mut)[2]){
  relapse_day_idx<-grep('Survival',validate_EGFR_mut[1:35,i])
  relapse_day[i]<-unname(as.numeric(strsplit(validate_EGFR_mut[relapse_day_idx,i],'Survival \\(days\\): ')[[1]][2]))
}
sample_id<-as.character(validate_EGFR_mut[55,])
clinical_data_mut<-data.frame(sample_id,relapse_status,relapse_day)
expression_data<-validate_EGFR_mut[55:(dim(validate_EGFR_mut)[1]-1),]
expression_data[,1]<-as.character(GPL6480_annot$GENE_SYMBOL[match(expression_data[,1],as.character(GPL6480_annot$ID))])
exp_validate_EGFR_mut<-data.frame(t(expression_data),stringsAsFactors = F)
exp_validate_EGFR_mut<-cbind(clinical_data_mut,exp_validate_EGFR_mut[,-1])
exp_validate_header<-as.vector(t(exp_validate_EGFR_mut[1,]))
exp_validate_EGFR_mut<-exp_validate_EGFR_mut[-1,]
match_idx<-match(Data_table$Gene_symbol,exp_validate_header)
exp_validate_EGFR_mut_gsub<-exp_validate_EGFR_mut[,c(1:3,match_idx[!is.na(match_idx)])]
colnames(exp_validate_EGFR_mut_gsub)[-c(1:3)]<-exp_validate_header[match_idx[!is.na(match_idx)]]

### ER - 
relapse_status<-c()
for(i in 1:dim(validate_EGFR_wt)[2]){
  relapse_status_idx<-grep('Evidence of relapse: ',validate_EGFR_wt[1:30,i])
  relapse_status[i]<-grepl('Evidence of relapse: Y',validate_EGFR_wt[relapse_status_idx,i])
}
relapse_day<-c()
for(i in 2:dim(validate_EGFR_wt)[2]){
  relapse_day_idx<-grep('Survival',validate_EGFR_wt[1:35,i])
  relapse_day[i]<-unname(as.numeric(strsplit(validate_EGFR_wt[relapse_day_idx,i],'Survival \\(days\\): ')[[1]][2]))
}
sample_id<-as.character(validate_EGFR_wt[55,])
clinical_data_wt<-data.frame(sample_id,relapse_status,relapse_day)
expression_data<-validate_EGFR_wt[55:(dim(validate_EGFR_wt)[1]-1),]
expression_data[,1]<-as.character(GPL6480_annot$GENE_SYMBOL[match(expression_data[,1],as.character(GPL6480_annot$ID))])
exp_validate_EGFR_wt<-data.frame(t(expression_data),stringsAsFactors = F)
exp_validate_EGFR_wt<-cbind(clinical_data_wt,exp_validate_EGFR_wt[,-1])
exp_validate_header<-as.vector(t(exp_validate_EGFR_wt[1,]))
exp_validate_EGFR_wt<-exp_validate_EGFR_wt[-1,]

match_idx<-match(Data_table$Gene_symbol,exp_validate_header)
exp_validate_EGFR_wt_gsub<-exp_validate_EGFR_wt[,c(1:3,match_idx[!is.na(match_idx)])]
colnames(exp_validate_EGFR_wt_gsub)[-c(1:3)]<-exp_validate_header[match_idx[!is.na(match_idx)]]

####### WHALA
GSE13213_mut<-exp_validate_EGFR_mut_gsub
GSE13213_wt<-exp_validate_EGFR_wt_gsub
#######







######## TCGA
library(survival)
setwd('/Users/NeilWu/Github/Research_EGFRmut_survival/TCGA/')
TCGA_p<-read.table('nationwidechildrens.org_clinical_patient_luad.txt',header = T,sep='\t',stringsAsFactors = F)
TCGA_p<-TCGA_p[-c(1,2),]
TCGA_f<-read.table('nationwidechildrens.org_clinical_follow_up_v1.0_luad.txt',header=T,sep='\t',stringsAsFactors = F)
TCGA_f<-TCGA_f[-c(1,2),]
TCGA_m<-read.table('mutation.txt',header=T,sep='\t',stringsAsFactors = F)
TCGA_srdf<-read.table('unc.edu_LUAD.IlluminaHiSeq_RNASeqV2.1.15.0.sdrf.txt',header = T,sep='\t',stringsAsFactors = F)

#even_death<-grepl('Dead',TCGA_f$vital_status)
event_date<-c()
even_relapse<-!grepl('[Not Applicable]',TCGA_f$new_tumor_event_dx_days_to)
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

TCGA_data<-TCGA_data[,c(1,TCGA_clinical$Derived_file_idx[!is.na(TCGA_clinical$Derived_file_idx)])]
TCGA_data<-data.frame(t(TCGA_data),stringsAsFactors = F)

TCGA_sig_idx<-match(Data_table$Gene_symbol,TCGA_data[1,])
TCGA_sig_idx<-c(TCGA_sig_idx[!is.na(TCGA_sig_idx)])

TCGA_data_dr<-TCGA_data[,TCGA_sig_idx[!duplicated(TCGA_sig_idx)]]
colnames(TCGA_data_dr)<-TCGA_data_dr[1,]
TCGA_data_dr<-TCGA_data_dr[-1,]
TCGA_data_dr<-TCGA_data_dr[grepl('results$',row.names(TCGA_data_dr)),]
TCGA_data_dr$Derived_file<-substr(row.names(TCGA_data_dr),1,16)

TCGA_clinical_sub<-data.frame(Relapse=TCGA_clinical$Event_relapse,
                              Date=TCGA_clinical$Event_day,
                              Mutate=TCGA_clinical$mutate,
                              Stage=TCGA_clinical$Stage,
                              Derived_file=substr(TCGA_clinical$Derived_file,1,16))

TCGA_data_w_clinic<-merge(TCGA_data_dr,TCGA_clinical_sub,by='Derived_file')
TCGA_data_w_clinic$Date<-as.numeric(as.character(TCGA_data_w_clinic$Date))
TCGA_data_w_clinic<-TCGA_data_w_clinic[!is.na(TCGA_data_w_clinic$Date),]

TCGA_data_w_clinic_mut<-TCGA_data_w_clinic[as.character(TCGA_data_w_clinic$Mutate)=='EGFR',]
TCGA_data_w_clinic_mut<-TCGA_data_w_clinic_mut[!grepl('Stage III',as.character(TCGA_data_w_clinic_mut$Stage)),]

TCGA_data_w_clinic_wt<-TCGA_data_w_clinic[as.character(TCGA_data_w_clinic$Mutate)!='EGFR',]
TCGA_data_w_clinic_wt<-TCGA_data_w_clinic_wt[!grepl('Stage III',as.character(TCGA_data_w_clinic_wt$Stage)),]

####  WHALA
TCGA_wt<-TCGA_data_w_clinic_wt
TCGA_mut<-TCGA_data_w_clinic_mut
####








####### GSE31210
temp<-read.table('/Users/NeilWu/Github/Research_EGFRmut_survival/GSE31210/GSE31210_series_matrix.txt',
                 skip = 25,fill=T,stringsAsFactors = F)
U133plus2<-read.csv('/Users/NeilWu/Github/Research_lncRNA-AML/HG-U133_Plus_2.na34.annot.csv',stringsAsFactors=F)
idx_EGFR_mut<-grep('EGFR mutation +',temp[17,])
idx_EGFR_wt<-which(grepl('gene alteration status:',temp[17,]) & !grepl('EGFR mutation +',temp[17,]))
temp_sub_wt<-temp[,c(1,idx_EGFR_wt)]
temp_sub_mut<-temp[,c(1,idx_EGFR_mut)]


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
expression_data[,1]<-matched_probe_done
###

exp_temp_EGFR_mut<-data.frame(t(expression_data),stringsAsFactors = F)
exp_temp_EGFR_mut<-cbind(clinical_data_mut,exp_temp_EGFR_mut[,-1])
exp_header<-as.vector(t(exp_temp_EGFR_mut[1,]))
exp_temp_EGFR_mut<-exp_temp_EGFR_mut[-1,]


################## Header add
colnames(exp_temp_EGFR_wt)[-c(1:3)]<-exp_header[-c(1:3)]
colnames(exp_temp_EGFR_mut)[-c(1:3)]<-exp_header[-c(1:3)]
##################


###### WHALA
GSE31210_mut<-exp_temp_EGFR_mut
GSE31210_wt<-exp_temp_EGFR_wt
######


rm(list = grep('GSE13213_mut|GSE13213_wt|TCGA_wt|TCGA_mut|GSE31210_wt|GSE31210_mut',
               ls(),invert = T,value = T))


