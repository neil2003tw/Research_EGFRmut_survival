setwd('/Users/NeilWu/Github/Research_EGFRmut_survival/')

z_trans<-function(l){
  l<-as.numeric(l)
  x<-l[!is.na(l)]
  y<-(x-mean(x))/sd(x)
  l[!is.na(l)]<-y
  return(l)
}

Data_table<-read.table('Data_0317.txt',sep='\t',header = T)

list_ESPG_TCGA<-(Data_table$TCGAp_p<0.05) & (Data_table$TCGAn_p>0.05)
list_ESPG_TCGA<-which(list_ESPG_TCGA)

list_ESPG_GSE13213<-(Data_table$GSE13213p_p<0.05) & (Data_table$GSE13213n_p>0.05)
list_ESPG_GSE13213<-which(list_ESPG_GSE13213)

Data_table_val<-Data_table[list_ESPG_GSE13213[list_ESPG_GSE13213 %in% list_ESPG_TCGA],]
Data_table_val[order(Data_table_val$EGFRp_cox_beta,decreasing = T),]
Data_table_val<-Data_table_val[,c(1,grep('p_',colnames(Data_table_val)))]
Data_table_beta<-Data_table_val[,c(1,grep('beta',colnames(Data_table_val)))]

Data_table_beta_z<-data.frame(Gene_symbol=Data_table_beta[,1],stringsAsFactors = F)
Data_table_beta_z$GSE31210_p<-NA
Data_table_beta_z$GSE31210_beta<-NA
Data_table_beta_z$TCGA_p<-NA
Data_table_beta_z$TCGA_beta<-NA
Data_table_beta_z$GSE13213_p<-NA
Data_table_beta_z$GSE13213_beta<-NA

for(i in 1:nrow(Data_table_beta_z)){
  gene_n<-as.character(Data_table_beta_z$Gene_symbol[i])
  cox_31210<-coxph(Surv(time = relapse_day/30, event = relapse_status) ~ z_trans(GSE31210_mut[,gene_n]),data=GSE31210_mut)
  cox_13213<-coxph(Surv(time = relapse_day/30, event = relapse_status) ~ z_trans(GSE13213_mut[,gene_n]),data=GSE13213_mut)
  cox_TCGA<-coxph(Surv(time = Date/30, event = Relapse) ~ z_trans(TCGA_mut[,gene_n]),data=TCGA_mut)
  Data_table_beta_z$GSE31210_p[i]<-summary(cox_31210)$coefficients[5]
  Data_table_beta_z$GSE31210_beta[i]<-summary(cox_31210)$coefficients[1]
  Data_table_beta_z$GSE13213_p[i]<-summary(cox_13213)$coefficients[5]
  Data_table_beta_z$GSE13213_beta[i]<-summary(cox_13213)$coefficients[1]
  Data_table_beta_z$TCGA_p[i]<-summary(cox_TCGA)$coefficients[5]
  Data_table_beta_z$TCGA_beta[i]<-summary(cox_TCGA)$coefficients[1]
}




