
#################### GSE13213 validate
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

Data_table$GSE13213p_p<-NA
Data_table$GSE13213p_beta<-NA
Data_table$GSE13213p_se<-NA
for(i in 3:dim(Data_table)[1]){
  now_gene_list<-strsplit(Data_table$Gene_symbol[i],split = ' /// ')[[1]]
  if(any(now_gene_list %in% colnames(exp_validate_EGFR_mut_gsub))){
    now_gene<-now_gene_list[now_gene_list %in% colnames(exp_validate_EGFR_mut_gsub)][1]
    if(!all(is.na(as.numeric(exp_validate_EGFR_mut_gsub[,now_gene])))){
    cox.summary.v<-summary(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(exp_validate_EGFR_mut_gsub[,now_gene]), data = exp_validate_EGFR_mut_gsub))      
    Data_table$GSE13213p_p[i]<-cox.summary.v$coefficient[5]
    Data_table$GSE13213p_beta[i]<-cox.summary.v$coefficient[1]
    Data_table$GSE13213p_se[i]<-cox.summary.v$coefficient[3]
  }}
}




#### ER -
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

Data_table$GSE13213n_p<-NA
Data_table$GSE13213n_beta<-NA
Data_table$GSE13213n_se<-NA
for(i in 3:dim(Data_table)[1]){
  now_gene_list<-strsplit(Data_table$Gene_symbol[i],split = ' /// ')[[1]]
  if(any(now_gene_list %in% colnames(exp_validate_EGFR_mut_gsub))){
    now_gene<-now_gene_list[now_gene_list %in% colnames(exp_validate_EGFR_mut_gsub)][1]
    cox.summary.v<-summary(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(exp_validate_EGFR_wt_gsub[,now_gene]), data = exp_validate_EGFR_wt_gsub))      
    Data_table$GSE13213n_p[i]<-cox.summary.v$coefficient[5]
    Data_table$GSE13213n_beta[i]<-cox.summary.v$coefficient[1]
    Data_table$GSE13213n_se[i]<-cox.summary.v$coefficient[3]
  }
}
####