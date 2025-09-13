rm(list=ls());gc()
setwd('/project/CRI/DeBerardinis_lab/lcai/ForOthers/ZhengWu/06112024/data/')
library(openxlsx);library(data.table);library(ggplot2)
######## total ########
p.total<-read.xlsx('PCF-ZW-7426--05--2024_total_PD1.xlsx',sheet=1,skipEmptyRows = F);setDT(p.total)
table(p.total$Protein.FDR.Confidence.Combined,useNA = 'ifany')
which(is.na(p.total$Protein.FDR.Confidence.Combined))
p.total<-rbind(data.table(color='white',p.total[1:6604]),data.table(color='orange',p.total[6608:7339]))
dat.cols<-grep("Abundance",colnames(p.total))
p.total.features<-p.total[,1:(dat.cols[1]-1),with=F]
p.total.data<-p.total[,dat.cols,with=F]
p.total.data[p.total.data==0]<-NA
use.ind<-which(rowSums(is.na(p.total.data))==0)
norm.factor<-colSums(p.total.data[use.ind,])/mean(colSums(p.total.data[use.ind,]))
p.total.data<-as.matrix(p.total.data)
p.total.data[]<-t(apply(p.total.data,1,function(x) x/norm.factor))
table(table(p.total.features$Accession))
table(table(p.total.features$Gene.Symbol))
rownames(p.total.data)<-p.total.features$Accession
norm.factor.A<-norm.factor

######## phospho ########
p.ppt<-read.xlsx('PCF-ZW-7426--05--2024_phosphoprotein_PD1.xlsx',skipEmptyRows = F);setDT(p.ppt)
table(p.ppt$Protein.FDR.Confidence.Combined,useNA = 'ifany')
which(is.na(p.ppt$Protein.FDR.Confidence.Combined))
p.ppt<-rbind(data.table(color='white',p.ppt[1:700]),data.table(color='orange',p.ppt[704:1201]))
dat.cols<-grep("Abundance",colnames(p.ppt))
p.ppt.features<-p.ppt[,1:(dat.cols[1]-1),with=F]
p.ppt.data<-p.ppt[,dat.cols,with=F]
norm.factor<-colSums(p.ppt.data)/mean(colSums(p.ppt.data))
p.ppt.data<-as.matrix(p.ppt.data)
p.ppt.data[]<-t(apply(p.ppt.data,1,function(x) x/norm.factor))
table(table(p.ppt.features$Accession))
table(table(p.ppt.features$Gene.Symbol))
rownames(p.ppt.data)<-p.ppt.features$Accession
norm.factor.B<-norm.factor

p.ppd<-read.xlsx('PCF-ZW-7426--05--2024_phosphopeptide_PD1.xlsx');setDT(p.ppd)
dat.cols<-grep("Abundance",colnames(p.ppd))
p.ppd.data<-p.ppd[,dat.cols,with=F]
table(rowSums(is.na(p.ppd.data)))
table(rowSums(p.ppd.data==0))
p.ppd.data[which(rowSums(p.ppd.data==0)!=0),]
length(grep('Phospho',p.ppd$Modifications,invert = T)) # 64 before filtering by non-missing
non.missing.ind<-which(rowSums(is.na(p.ppd.data))==0)
p.ppd<-p.ppd[non.missing.ind]
p.ppd.features<-p.ppd[,1:(dat.cols[1]-1),with=F]
p.ppd.data<-p.ppd.data[non.missing.ind]
use.ind<-which(rowSums(is.na(p.ppd.data))==0)
norm.factor<-colSums(p.ppd.data[use.ind,])/mean(colSums(p.ppd.data[use.ind,]))
norm.factor.C<-norm.factor
p.ppd.data<-as.matrix(p.ppd.data)
p.ppd.data[]<-t(apply(p.ppd.data,1,function(x) x/norm.factor))
table(table(p.ppd.features$Master.Protein.Accessions))

pairs(do.call(cbind,list(norm.factor.A,norm.factor.B,norm.factor.C)))
# unify
identical(colnames(p.total.data),colnames(p.ppt.data))
identical(colnames(p.total.data),colnames(p.ppd.data))
unify.names<-gsub("Abundance.",'',colnames(p.total.data),fixed = T)
colnames(p.total.data)<-colnames(p.ppt.data)<-colnames(p.ppd.data)<-unify.names

length(grep('Phospho',p.ppd.features$Modifications,invert = T))
# do we need to log transform?
par(mfrow=c(5,5));for(i in sample(1:nrow(p.total.data),25))  plot(density(p.total.data[i,]),main=i)
par(mfrow=c(5,5));for(i in sample(1:nrow(p.ppt.data),25))  plot(density(p.ppd.data[i,]),main=i)
par(mfrow=c(5,5));for(i in sample(1:nrow(p.ppd.data),25))  plot(density(p.ppd.data[i,]),main=i)
# no

########### STATS ############
library(dplyr);library(tidyr)
library(effsize)
source('../script/beautify_df.R')
make.stat<-function(m.dat,simple=F){
  m.dat%>%separate(sample,into=c('genotype','treatment','replicate'),remove = F,sep = '-')%>%setDT->m.dat
  m.dat[,treatment:=factor(treatment,levels=c('MOCK','6TG'))][,genotype:=factor(genotype,levels=c('WT','NUDT5'))]
  if(!simple){
    # treatment within WT
    WT.stat<-m.dat[genotype=='WT',as.list(structure(tryCatch({c(with(t.test(value~treatment),c(p.value,estimate[2]/estimate[1])),with(cohen.d(value~treatment),-estimate))},error=function(e) rep(as.numeric(NA),3)),names=c('WT.trt.pv','WT.trt.FC','WT.trt.effsize'))),by=acc]
    KO.stat<-m.dat[genotype=='NUDT5',as.list(structure(tryCatch({c(with(t.test(value~treatment),c(p.value,estimate[2]/estimate[1])),with(cohen.d(value~treatment),-estimate))},error=function(e) rep(as.numeric(NA),3)),names=c('KO.trt.pv','KO.trt.FC','KO.trt.effsize'))),by=acc]
    untreat.stat<-m.dat[treatment=='MOCK',as.list(structure(tryCatch({c(with(t.test(value~genotype),c(p.value,estimate[2]/estimate[1])),with(cohen.d(value~genotype),-estimate))},error=function(e) rep(as.numeric(NA),3)),names=c('mock.KO.pv','mock.KO.FC','mock.KO.effsize'))),by=acc]
    treat.stat<-m.dat[treatment=='6TG',as.list(structure(tryCatch({c(with(t.test(value~genotype),c(p.value,estimate[2]/estimate[1])),with(cohen.d(value~genotype),-estimate))},error=function(e) rep(as.numeric(NA),3)),names=c('6TG.KO.pv','6TG.KO.FC','6TG.KO.effsize'))),by=acc]
    WT.stat[,WT.trt.padj:=p.adjust(WT.trt.pv,'BH')]
    KO.stat[,KO.trt.padj:=p.adjust(KO.trt.pv,'BH')]
    untreat.stat[,mock.KO.padj:=p.adjust(mock.KO.pv,'BH')]
    treat.stat[,`6TG.KO.padj`:=p.adjust(`6TG.KO.pv`,'BH')]
    
    # two-way ANOVA
    stat.df<-m.dat[,tryCatch({melt(coef(summary(lm(value~genotype*treatment)))[2:4,3:4])},error=function(e){print(acc)}),by=acc]
    stat.list<-list()
    for(l1 in levels(stat.df$Var1)){
      for(l2 in levels(stat.df$Var2)){
        stat.list[[paste(l1,l2,sep='.')]]<-stat.df[Var1==l1 & Var2==l2][['value']]
        if(grepl('Pr',l2)) stat.list[[paste(l1,'padj',sep='.')]]<-p.adjust(stat.df[Var1==l1 & Var2==l2][['value']],'BH')
      }
    }
    
    stat.df2<-data.table(acc=as.character(stat.df[Var1==l1 & Var2==l2][['acc']]),do.call(cbind,stat.list))
    colnames(stat.df2)[grep('Pr',colnames(stat.df2))]<-gsub('Pr(>|t|)','pv',colnames(stat.df2)[grep('Pr',colnames(stat.df2))])
    
    stat.df3<-merge(stat.df2,WT.stat,by='acc')
    stat.df3<-merge(stat.df3,KO.stat,by='acc')
    stat.df3<-merge(stat.df3,untreat.stat,by='acc')
    stat.df3<-merge(stat.df3,treat.stat,by='acc')
    stat.df3<-beautify.dt(stat.df3)
    list(m.dat,stat.df3)
  }else{
    m.dat
  }
}

# 1. protein differing between the two
m.total<-as.data.table(melt(p.total.data));colnames(m.total)[1:2]<-c('acc','sample')
total.results<-make.stat(m.total)
m.total<-total.results[[1]];total.stats<-total.results[[2]]

# 55 phosphoproteins not in total protein
table(rownames(p.ppt.data)%in%rownames(p.total.data))
p.ppt.data.normalized<-p.ppt.data[which(rownames(p.ppt.data)%in%rownames(p.total.data)),]
p.ppt.data.unnormalized<-p.ppt.data[which(!rownames(p.ppt.data)%in%rownames(p.total.data)),]
un.m.ppt<-rbind(data.table(normalization='no',melt(p.ppt.data.normalized)),data.table(normalization='no',melt(p.ppt.data.unnormalized)))
p.ppt.data.normalized<-p.ppt.data.normalized/p.total.data[match(rownames(p.ppt.data.normalized),rownames(p.total.data))]
m.ppt<-rbind(data.table(normalization='yes',melt(p.ppt.data.normalized)),data.table(normalization='no',melt(p.ppt.data.unnormalized)))
colnames(un.m.ppt)[2:3]<-colnames(m.ppt)[2:3]<-c('acc','sample')
un.m.ppt<-make.stat(un.m.ppt,T)
ppt.results<-make.stat(m.ppt)
m.ppt<-ppt.results[[1]];ppt.stats<-ppt.results[[2]]


# 68 phosphoproteins not in total protein
table(unique(p.ppd.features$Master.Protein.Accessions)%in%rownames(p.total.data))
matched.ind<-which(p.ppd.features$Master.Protein.Accessions%in%rownames(p.total.data))
unmatched.ind<-which(!p.ppd.features$Master.Protein.Accessions%in%rownames(p.total.data))

p.ppd.data.normalized<-p.ppd.data[matched.ind,]
p.ppd.data.unnormalized<-p.ppd.data[unmatched.ind,]
un.m.ppd<-rbind(data.table(normalization='no',melt(p.ppd.data.normalized)),data.table(normalization='no',melt(p.ppd.data.unnormalized)))
p.ppd.data.normalized<-p.ppd.data.normalized/p.total.data[match(p.ppd.features$Master.Protein.Accessions[matched.ind],rownames(p.total.data)),]
m.ppd<-rbind(data.table(normalization='yes',melt(p.ppd.data.normalized)),data.table(normalization='no',melt(p.ppd.data.unnormalized)))
p.ppd.features<-p.ppd.features[c(matched.ind,unmatched.ind)]
colnames(un.m.ppd)[2:3]<-colnames(m.ppd)[2:3]<-c('acc','sample')
m.ppd[,acc:=paste0('V',acc)];m.ppd[,acc:=factor(acc,levels=unique(acc))]
un.m.ppd[,acc:=paste0('V',acc)];un.m.ppd[,acc:=factor(acc,levels=unique(acc))]
ppd.results<-make.stat(m.ppd)
m.ppd<-ppd.results[[1]];ppd.stats<-ppd.results[[2]]
un.m.ppd<-make.stat(un.m.ppd,T)
total.stats<-data.table(p.total.features[match(total.stats$acc,Accession)],total.stats[,-1,with=F])
ppt.stats<-data.table(p.ppt.features[match(ppt.stats$acc,Accession)],normalization=m.ppt[match(ppt.stats$acc,acc)][['normalization']],ppt.stats[,-1,with=F])
ppd.stats<-data.table(p.ppd.features[as.integer(gsub('V','',ppd.stats$acc))],normalization=m.ppd[match(ppd.stats$acc,acc)][['normalization']],ppd.stats)
dir.create('../output')
setwd('../output/')


wb<-createWorkbook()
addWorksheet(wb,'total');writeDataTable(wb,'total',data.table(p.total.features[match(rownames(p.total.data),Accession)],p.total.data))
p.ppt.cbn<-rbind(p.ppt.data.normalized,p.ppt.data.unnormalized)
p.ppt.cbn<-data.table(p.ppt.features[match(rownames(p.ppt.cbn),Accession)],
                      normalization_to_total_protein=rep(c('yes','no'),c(nrow(p.ppt.data.normalized),nrow(p.ppt.data.unnormalized))),
                      p.ppt.cbn)
addWorksheet(wb,'phosphoprotein');writeDataTable(wb,'phosphoprotein',data.table(p.ppt.features,p.ppt.data))
addWorksheet(wb,'phosphoprotein_totalNormalized');writeDataTable(wb,'phosphoprotein_totalNormalized',p.ppt.cbn)
p.ppd.cbn<-rbind(p.ppd.data.normalized,p.ppd.data.unnormalized)
p.ppd.cbn<-data.table(p.ppd.features,
                      normalization_to_total_protein=rep(c('yes','no'),c(nrow(p.ppd.data.normalized),nrow(p.ppd.data.unnormalized))),
                      p.ppd.cbn)
addWorksheet(wb,'phosphopeptide');writeDataTable(wb,'phosphopeptide',data.table(p.ppd.features,p.ppd.data[c(matched.ind,unmatched.ind),]))
addWorksheet(wb,'phosphopeptide_totalNormalized');writeDataTable(wb,'phosphopeptide_totalNormalized',p.ppd.cbn)
saveWorkbook(wb,'normalized_data.xlsx',overwrite = T)

View(ppd.stats)
p.ppd.features[,Accession:=paste0('V',1:nrow(p.ppd.features))]
g.ref<-rbind(p.total.features[,c('Gene.Symbol','Accession')],p.ppt.features[,c('Gene.Symbol','Accession')])[!duplicated(Accession)]
mk.plot<-function(gene='',Acc=''){
  if(gene!='')  Acc<-g.ref[Gene.Symbol==gene][['Accession']]
  if(Acc!='')  gene<-g.ref[Accession%in%Acc][['Gene.Symbol']]
  if(length(Acc)!=0){
    if(Acc%in%m.total$acc){
      g1<-ggplot(m.total[acc%in%Acc],aes(x=treatment,y=value))+geom_point(alpha=.3)+facet_grid(acc~genotype)+theme_bw()+ggtitle('Total protein')
    }else{
      g1<-plot_spacer()
    }
    
    un.m.ppt.sub<-un.m.ppt[acc%in%Acc];un.m.ppt.sub[,acc:=paste0(acc,' (',ifelse(normalization=='yes','normalized','unnormalized'),')')]
    if(nrow(un.m.ppt.sub)>0){
      g2.1<-ggplot(un.m.ppt.sub,aes(x=treatment,y=value))+geom_point(alpha=.3)+facet_grid(acc~genotype)+theme_bw()+ggtitle('Phosphoprotein')
    }else{
      g2.1<-plot_spacer()
    }
    
    m.ppt.sub<-m.ppt[acc%in%Acc];m.ppt.sub[,acc:=paste0(acc,' (',ifelse(normalization=='yes','normalized','unnormalized'),')')]
    if(nrow(m.ppt.sub)>0){
      g2.2<-ggplot(m.ppt.sub,aes(x=treatment,y=value))+geom_point(alpha=.3)+facet_grid(acc~genotype)+theme_bw()+ggtitle('Phosphoprotein')
    }else{
      g2.2<-plot_spacer()
    }
    
    ppd.ind<-with(p.ppd.features,grep(paste(Acc,collapse =  '|'),Master.Protein.Accessions))
    un.m.ppd.sub<-un.m.ppd[acc%in%paste0('V',ppd.ind)]
    if(nrow(un.m.ppd.sub)>0){
      un.m.ppd.sub[,acc:=with(p.ppd.features[match(acc,Accession)],paste(Master.Protein.Accessions,Modifications,sep="\n"))]
      un.m.ppd.sub[,acc:=paste0(acc,"\n(",ifelse(normalization=='yes','normalized','unnormalized'),')')]
      g3.1<-ggplot(un.m.ppd.sub,aes(x=treatment,y=value))+geom_point(alpha=.3)+facet_grid(acc~genotype)+theme_bw()+ggtitle('Phosphopeptide')
    }else{
      g3.1<-plot_spacer()
    }
    
    m.ppd.sub<-m.ppd[acc%in%paste0('V',ppd.ind)]
    if(nrow(m.ppd.sub)>0){
      m.ppd.sub[,acc:=with(p.ppd.features[match(acc,Accession)],paste(Master.Protein.Accessions,Modifications,sep="\n"))]
      m.ppd.sub[,acc:=paste0(acc,"\n(",ifelse(normalization=='yes','normalized','unnormalized'),')')]
      g3.2<-ggplot(m.ppd.sub,aes(x=treatment,y=value))+geom_point(alpha=.3)+facet_grid(acc~genotype)+theme_bw()+ggtitle('Phosphopeptide')
    }else{
      g3.2<-plot_spacer()
    }
    g1+g2.1+g2.2+g3.1+g3.2+plot_annotation(title = paste(gene,Acc))
    
  }else{print('not found')}
}
dir.create('../work')
mk.plot('CAD')
mk.plot('TPR')
mk.plot(gene = 'TUFM') # interesting
mk.plot(Acc = 'O43175')

source('../script/HyperGeometricTest.R')
pw.list<-readRDS('pathway_list.rds')



ppd.stats[,Gene.Symbol:=g.ref[match(ppd.stats$Master.Protein.Accessions,g.ref$Accession)][['Gene.Symbol']]]
ppd.stats[,acc:=NULL]
annot.stat<-function(in.stat){
  in.stat[,trt.group:='']
  in.stat[WT.trt.padj<.05 & WT.trt.effsize<0  & WT.trt.effsize < KO.trt.effsize & ((WT.trt.padj<KO.trt.padj)|KO.trt.effsize>0),
          trt.group:=paste0(trt.group,', more down with treatment in WT')]
  in.stat[WT.trt.padj<.05 & WT.trt.effsize>0  & WT.trt.effsize > KO.trt.effsize & ((WT.trt.padj<KO.trt.padj)|KO.trt.effsize<0),
          trt.group:=paste0(trt.group,', more up with treatment in WT')]
  in.stat[KO.trt.padj<.05 & KO.trt.effsize<0 & KO.trt.effsize < WT.trt.effsize & ((KO.trt.padj<WT.trt.padj)|(WT.trt.effsize>0)),
          trt.group:=paste0(trt.group,', more down with treatment in KO')]
  in.stat[KO.trt.padj<.05 & KO.trt.effsize>0 & KO.trt.effsize > WT.trt.effsize & ((KO.trt.padj<WT.trt.padj)|(WT.trt.effsize<0)),
          trt.group:=paste0(trt.group,', more up with treatment in KO')]
  
  in.stat[,genetic.group:='']
  in.stat[mock.KO.padj<.05 & mock.KO.effsize<0  & mock.KO.effsize < `6TG.KO.effsize` & ((mock.KO.padj < `6TG.KO.padj`)|`6TG.KO.effsize`>0),
          genetic.group:=paste0(genetic.group,', more down with KO under mock treatment')]
  in.stat[mock.KO.padj<.05 & mock.KO.effsize>0  & mock.KO.effsize > `6TG.KO.effsize` & ((mock.KO.padj < `6TG.KO.padj`)|`6TG.KO.effsize`<0),
          genetic.group:=paste0(genetic.group,', more up with KO under mock treatment')]
  in.stat[`6TG.KO.padj`<.05 & `6TG.KO.effsize`<0 & `6TG.KO.effsize` < mock.KO.effsize & ((`6TG.KO.padj`<mock.KO.padj)|(mock.KO.effsize>0)),
          genetic.group:=paste0(genetic.group,', more down with KO under 6TG treatment')]
  in.stat[`6TG.KO.padj`<.05 & `6TG.KO.effsize`>0 & `6TG.KO.effsize` > mock.KO.effsize & ((`6TG.KO.padj`<mock.KO.padj)|(mock.KO.effsize<0)),
          genetic.group:=paste0(genetic.group,', more up with KO under 6TG treatment')]
  
  in.stat[,group:='']
  in.stat[genotypeNUDT5.padj<.05 & `genotypeNUDT5.t value`>0,group:=paste0(group,", up in KO")]
  in.stat[genotypeNUDT5.padj<.05 & `genotypeNUDT5.t value`<0,group:=paste0(group,", down in KO")]
  in.stat[treatment6TG.padj<.05 & `treatment6TG.t value`>0,group:=paste0(group,", up with trt")]
  in.stat[treatment6TG.padj<.05 & `treatment6TG.t value`<0,group:=paste0(group,", down with trt")]
  in.stat[`genotypeNUDT5:treatment6TG.padj`<.05 & `genotypeNUDT5:treatment6TG.t value`>0,group:=paste0(group,", up intr")]
  in.stat[`genotypeNUDT5:treatment6TG.padj`<.05 & `genotypeNUDT5:treatment6TG.t value`<0,group:=paste0(group,", down intr")]
  
  in.stat[,genetic.group:=gsub("^, ",'',genetic.group)][,trt.group:=gsub("^, ",'',trt.group)][,group:=gsub("^, ",'',group)]
  in.stat
}
total.stats<-annot.stat(total.stats)
ppt.stats<-annot.stat(ppt.stats)
ppd.stats<-annot.stat(ppd.stats)

wb<-createWorkbook()
addWorksheet(wb,'total');writeDataTable(wb,'total',total.stats)
addWorksheet(wb,'phosphoprotein');writeDataTable(wb,'phosphoprotein',ppt.stats)
addWorksheet(wb,'phosphopeptide');writeDataTable(wb,'phosphopeptide',ppd.stats)
saveWorkbook(wb,'stats.xlsx',overwrite = T)

total.stats[,.N,by=group]
total.stats[,.N,by=trt.group]
total.stats[,.N,by=genetic.group]
ggplot(total.stats[,.N,by=.(trt.group,genetic.group)],aes(x=trt.group,y=genetic.group))+geom_tile(mapping=aes(fill=N))+
  geom_text(color='white',mapping=aes(label=N))+ggtitle('total protein')+
  theme_classic()+scale_fill_viridis_c()+theme(legend.position = 'none',axis.text.x = element_text(angle = 30,vjust = 1,hjust=1))
total.stats[genetic.group=='more down with KO under mock treatment, more up with KO under 6TG treatment']
for(g in setdiff(unique(total.stats$trt.group),'')) {cat(paste0(g,'\n'));print(total.stats[trt.group==g][1:5][['Gene.Symbol']])}
for(g in setdiff(unique(total.stats$genetic.group),'')) {cat(paste0(g,'\n'));print(total.stats[genetic.group==g][1:5][['Gene.Symbol']])}

wrap_plots(lapply(c('NUDT4B','SMIM26','ACSL1','SNX9'),function(g){
  ggplot(m.dat[acc%in%g.ref[match(g,Gene.Symbol)][['Accession']]],aes(x=treatment,y=value))+geom_point(alpha=.3)+facet_grid(acc~genotype)+theme_bw()+
    ggtitle(g)
}),ncol = 2)

wrap_plots(lapply(c('NUDT4B','SMIM26','GATD3B','NUDT19','CNOT1','NQO1'),function(g){
  ggplot(m.dat[acc%in%g.ref[match(g,Gene.Symbol)][['Accession']]],aes(x=genotype,y=value))+geom_point(alpha=.3)+facet_grid(acc~treatment)+theme_bw()+
    ggtitle(g)
}),ncol = 3)

for(i in 1:3){
  print(i)
  in.stat<-switch(i,total.stats,ppt.stats,ppd.stats)
  for(big.gp in grep('group',colnames(in.stat),value = T)){
    print(big.gp)
    small.groups<-unique(unlist(sapply(unique(in.stat[[big.gp]]),function(x) unlist(strsplit(x,split = ', ')))))
    for(small.gp in small.groups){
      print(small.gp)
      gs<-in.stat[grepl(small.gp,in.stat[[big.gp]])][['Gene.Symbol']]
      if(length(gs)>=10){
        wb<-createWorkbook()
        for(pw in names(pw.list)[c(2,3,4,6,7,14:16,21)]){
          hyper.stat<-LC.hyper2(genesets = pw.list[[pw]],genes = unique(in.stat$Gene.Symbol),testset = gs,include.intersection = T)[order(pv)]
          hyper.stat[,padj:=p.adjust(pv,'BH')]
          hyper.stat<-hyper.stat[pv<.05]
          if(nrow(hyper.stat)>0){addWorksheet(wb,pw);writeDataTable(wb,pw,hyper.stat)}
        }
        saveWorkbook(wb,paste0(switch(i,'total_protein_','phosphoprotein_','phosphopeptide_'),big.gp,'_',small.gp,'.xlsx'),overwrite = T)
      }
    }
  }
}
save.image('../work/working.RData')
save(list=c('g.ref','m.total','un.m.ppt','m.ppt','un.m.ppd','m.ppd','p.ppd.features','total.stats','ppt.stats','ppd.stats'),file = '../work/plot_input.RData')
