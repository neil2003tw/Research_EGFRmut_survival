setwd('/Users/NeilWu/Github/Research_EGFRmut_survival/')
library(survival)
library(Hmisc)
source('/Users/NeilWu/Github/Research_EGFRmut_survival/Code/GSE31210.R')
source('/Users/NeilWu/Github/Research_EGFRmut_survival/Code/GSE13213_validate.R')
source('/Users/NeilWu/Github/Research_EGFRmut_survival/Code/TCGA_validate.R')

Data_table<-read.table('Data_0317.txt',sep='\t',header = T)
#write.table(Data_table,'/Users/NeilWu/Github/Research_EGFRmut_survival/Data_0317.txt',
#            sep = '\t',row.names=F)

list_ESPG_TCGA<-(Data_table$TCGAp_p<0.05) & (Data_table$TCGAn_p>0.05)
list_ESPG_TCGA<-which(list_ESPG_TCGA)

list_ESPG_GSE13213<-(Data_table$GSE13213p_p<0.05) & (Data_table$GSE13213n_p>0.05)
list_ESPG_GSE13213<-which(list_ESPG_GSE13213)

Data_table_val<-Data_table[list_ESPG_GSE13213[list_ESPG_GSE13213 %in% list_ESPG_TCGA],]
Data_table_val[order(Data_table_val$EGFRp_cox_beta,decreasing = T),]
#Data_table_val<-Data_table_val[Data_table_val$EGFRp_cox_p<0.05 & Data_table_val$EGFRn_cox_p>0.05,]
Data_table_val<-Data_table_val[,c(1,grep('p_',colnames(Data_table_val)))]
Data_table_beta<-Data_table_val[,c(1,grep('beta',colnames(Data_table_val)))]
write.csv(Data_table_beta,'Data_table_beta.csv')
##### 31210 to GO!!!
data_31210<-read.table('/Users/NeilWu/Github/Research_LungSCC_prognosis/GSE31210_main/Data_merged_31210_AD.txt',
                       sep='\t',stringsAsFactors = F)
data_31210sub<-data_31210[,c(as.character(Data_table_val$Gene_symbol),
                             'event_death','day_to_event_death')]
Module<-coxph(Surv(day_to_event_death,event_death) ~ .,data = data_31210sub)
data_31210sub$class<-as.numeric(cut2(predict(Module,data_31210sub),g=2))-1

data_31210_ready<-data_31210[,-c(ncol(data_31210)-1,ncol(data_31210))]
data_31210_ready$class<-data_31210sub$class
data_31210_ready<-t(data_31210_ready)
data_31210_ready<-cbind(row.names(data_31210_ready),data_31210_ready)
colnames(data_31210_ready)<-c('NAME',colnames(data_31210_ready[,-1]))
write.table(data_31210_ready,'GSE31210/Data_0318_GSEAuse.txt',sep='\t',row.names = F,quote = F)
write.table(t(data_31210sub$class),'GSE31210/Data_0318_GSEAuse.cls',sep='\t',
            row.names = F,quote=F,col.names = F)
