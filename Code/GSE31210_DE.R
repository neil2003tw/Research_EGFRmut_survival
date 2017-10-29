library(reshape)
setwd('/Users/NeilWu/Github/Research_EGFRmut_survival/Data/')
GSE31210_wt<-read.table('GSE31210_EGFRn.txt',sep='\t',header=T)
GSE31210_mut<-read.table('GSE31210_EGFRp.txt',sep='\t',header=T)

z_trans<-function(l){
  l<-as.numeric(l)
  x<-l[!is.na(l)]
  y<-(x-mean(x))/sd(x)
  l[!is.na(l)]<-y
  return(l)
}

t_pval<-c()
t_fc<-c()
pb<-txtProgressBar(1,54675,style = 3)
for(i in 1:54675){
  setTxtProgressBar(pb,i)
  t.sum<-t.test(GSE31210_mut[,i],GSE31210_wt[,i])
  t_pval[i]<-t.sum$p.value
  t_fc[i]<-t.sum$estimate[1]/t.sum$estimate[2]
}

deg_GSE31210<-colnames(GSE31210_mut)[t_pval<0.05]
GSE31210_mutsub<-GSE31210_mut[t_pval<0.05]
GSE31210_wtsub<-GSE31210_wt[t_pval<0.05]
GSE31210_desub<-rbind(GSE31210_wtsub,GSE31210_mutsub)
GSE31210_deflag<-c(rep('WT',nrow(GSE31210_wtsub)),rep('MT',nrow(GSE31210_mutsub)))

setwd('GSEA_DE/')
temp<-sapply(colnames(GSE31210_desub), function(x) strsplit(x,'\\.')[[1]][1],USE.NAMES = F)
GSE31210_desub<-GSE31210_desub[,!duplicated(temp)]
colnames(GSE31210_desub)<-temp[!duplicated(temp)]
write.table(t(GSE31210_desub)[-9825,],'GSE31210_DE_GSEA.txt',sep='\t',col.names = F,quote=F)
write.table(t(as.numeric(factor(GSE31210_deflag))-1),'GSE31210_DE_GSEA.cls',row.names = F,col.names = F,sep='\t')


t_pval[t_pval<0.05]<-p.adjust(t_pval[t_pval<0.05],method = 'bonferroni')

deg_GSE31210<-colnames(GSE31210_mut)[t_pval<0.05]
GSE31210_mutsub<-GSE31210_mut[t_pval<0.05]
GSE31210_wtsub<-GSE31210_wt[t_pval<0.05]
GSE31210_desub<-rbind(GSE31210_wtsub,GSE31210_mutsub)
GSE31210_deflag<-c(rep('WT',nrow(GSE31210_wtsub)),rep('MT',nrow(GSE31210_mutsub)))

setwd('GSEA_DE/')
temp<-sapply(colnames(GSE31210_desub), function(x) strsplit(x,'\\.')[[1]][1],USE.NAMES = F)
GSE31210_desub<-GSE31210_desub[,!duplicated(temp)]
colnames(GSE31210_desub)<-temp[!duplicated(temp)]
#write.table(t(GSE31210_desub)[-9825,],'GSE31210_DE_GSEA.txt',sep='\t',col.names = F,quote=F)
#write.table(t(as.numeric(factor(GSE31210_deflag))-1),'GSE31210_DE_GSEA.cls',row.names = F,col.names = F,sep='\t')


GSE31210_mutt10<-GSE31210_mut[,head(order(t_pval),13)[-c(5,10,11)]]
GSE31210_wtt10<-GSE31210_wt[,head(order(t_pval),13)[-c(5,10,11)]]
colnames(GSE31210_mutt10)<-c('ZDHHC11B','PABPCL1','ZDHHC11','ATP13A4','NPF10','PDIK1L','ATP13A5','MST1','TPAT1','URGCP')
colnames(GSE31210_wtt10)<-c('ZDHHC11B','PABPCL1','ZDHHC11','ATP13A4','NPF10','PDIK1L','ATP13A5','MST1','TPAT1','URGCP')

Big_frame<-rbind(
  data.frame( ExpressionValue=unlist(as.list(GSE31210_mutt10),use.names = F),
              Genes=unlist(lapply(colnames(GSE31210_mutt10),function(x) rep(x,127))),
              EGFR=rep('Mut',127*10)
      ),
  data.frame( ExpressionValue=unlist(as.list(GSE31210_wtt10),use.names = F),
              Genes=unlist(lapply(colnames(GSE31210_wtt10),function(x) rep(x,98))),
              EGFR=rep('Wt',98*10)
      )
)

ggplot(Big_frame, aes(x = Genes, y = ExpressionValue, fill = EGFR)) + geom_boxplot() + theme_bw()


