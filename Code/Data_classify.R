library(Hmisc)
Data_table<-read.table('/Users/NeilWu/Github/Research_EGFRmut_survival/Data/Data_21genes_z.txt',sep=' ')
Data_classtable<-data.frame(GeneSymbol=Data_table$Gene_symbol,GSE13213_mut=NA,
                            GSE13213_wt=NA,TCGA_mut=NA,TCGA_wt=NA,
                            GSE31210_wt=NA,GSE31210_mut=NA)

for(i in 1:nrow(Data_table)){
  n_gene<-as.character(Data_table$Gene_symbol[i])
  cox.predict<-predict(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_wt[,n_gene]), data = TCGA_wt),
                      TCGA_wt)
  Data_classtable$TCGA_wt[i]<-paste0(as.numeric(cut2(cox.predict,g=2))-1,collapse = ' ')
  
  cox.predict<-predict(coxph(Surv(time = Date, event = Relapse) ~ as.numeric(TCGA_mut[,n_gene]), data = TCGA_mut),
                       TCGA_mut)
  Data_classtable$TCGA_mut[i]<-paste0(as.numeric(cut2(cox.predict,g=2))-1,collapse = ' ')
  
  cox.predict<-predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE13213_wt[,n_gene]), data = GSE13213_wt),
                       GSE13213_wt)
  Data_classtable$GSE13213_wt[i]<-paste0(as.numeric(cut2(cox.predict,g=2))-1,collapse = ' ')
  
  cox.predict<-predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE13213_mut[,n_gene]), data = GSE13213_mut),
                       GSE13213_mut)
  Data_classtable$GSE13213_mut[i]<-paste0(as.numeric(cut2(cox.predict,g=2))-1,collapse = ' ')
  
  cox.predict<-predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE31210_wt[,n_gene]), data = GSE31210_wt),
                       GSE31210_wt)
  Data_classtable$GSE31210_wt[i]<-paste0(as.numeric(cut2(cox.predict,g=2))-1,collapse = ' ')

  cox.predict<-predict(coxph(Surv(time = relapse_day, event = relapse_status) ~ as.numeric(GSE31210_mut[,n_gene]), data = GSE31210_mut),
                       GSE31210_mut)
  Data_classtable$GSE31210_mut[i]<-paste0(as.numeric(cut2(cox.predict,g=2))-1,collapse = ' ')  
}

