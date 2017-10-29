library(survival)
setwd('/Users/NeilWu/Github/Research_EGFRmut_survival/TCGA/')
### TCGA validate
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

###

Data_table$TCGAp_p<-NA
Data_table$TCGAp_beta<-NA
Data_table$TCGAp_se<-NA
for(i in 1:dim(Data_table)[1]){
  #now_gene<-Data_table$Gene_symbol[i]
  now_gene_list<-strsplit(Data_table$Gene_symbol[i],split = ' /// ')[[1]]
  if(any(now_gene_list %in% colnames(TCGA_data_w_clinic_mut))){
    now_gene<-now_gene_list[now_gene_list %in% colnames(TCGA_data_w_clinic_mut)][1]
    cox.summary.t<-summary(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_data_w_clinic_mut[,now_gene]), data = TCGA_data_w_clinic_mut))
    Data_table$TCGAp_p[i]<-cox.summary.t$coefficient[5]
    Data_table$TCGAp_beta[i]<-cox.summary.t$coefficient[1]
    Data_table$TCGAp_se[i]<-cox.summary.t$coefficient[3]
  }
}

###
Data_table$TCGAn_p<-NA
Data_table$TCGAn_beta<-NA
Data_table$TCGAn_se<-NA
for(i in 1:dim(Data_table)[1]){
  now_gene_list<-strsplit(Data_table$Gene_symbol[i],split = ' /// ')[[1]]
  if(any(now_gene_list %in% colnames(TCGA_data_w_clinic_mut))){
    now_gene<-now_gene_list[now_gene_list %in% colnames(TCGA_data_w_clinic_mut)][1]
    cox.summary.t<-summary(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_data_w_clinic_wt[,now_gene]), data = TCGA_data_w_clinic_wt))
    Data_table$TCGAn_p[i]<-cox.summary.t$coefficient[5]
    Data_table$TCGAn_beta[i]<-cox.summary.t$coefficient[1]
    Data_table$TCGAn_se[i]<-cox.summary.t$coefficient[3]
  }
}

