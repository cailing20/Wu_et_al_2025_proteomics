rm(list=ls());gc()
library(data.table);library(ggplot2);library(patchwork);library(ComplexHeatmap)
setwd('/project/CRI/DeBerardinis_lab/shared/To_Zheng/06112024')
load('work/plot_input.RData')
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
mk.plot('TPR')
mk.plot(gene = 'HMGA1')
mk.plot(Acc = 'O43175')
########## HEATMAPS ###############
trt.grp.options<-c('more down with treatment in WT',
                   'more up with treatment in WT',
                   'more down with treatment in KO',
                   'more up with treatment in KO')
genetic.grp.options<-c('more down with KO under mock treatment',
                       'more up with KO under mock treatment',
                       'more down with KO under 6TG treatment',
                       'more up with KO under 6TG treatment')
grp.options<-c('down in KO','up in KO','down with trt','up with trt','down intr','up intr')

grouping<-c('trt.group','genetic.group','group')
dat.names<-c('total protein', 'phosphoprotein', 'phosphopeptide')
hm.plot<-function(data.choice=1,group=0,trend=0,genes=NULL,col.title=NULL){
  m.dat<-switch(data.choice,m.total,m.ppt,m.ppd)
  stats<-switch(data.choice,total.stats,ppt.stats,ppd.stats)
  genes<-unlist(strsplit(genes,split = ','))
  genes<-trimws(gsub('\n','',genes,fixed = T))
  accessions<-g.ref[match(genes,Gene.Symbol)][['Accession']]
  if(data.choice==3) accessions<-p.ppd.features[Master.Protein.Accessions%in%accessions][['Accession']]
  if(group%in%(1:3)&trend%in%(1:6)){
    trend<-switch(group,trt.grp.options,genetic.grp.options,grp.options)[[trend]]
    ind<-which(grepl(trend,stats[[grouping[[group]]]]))
    accessions<-intersect(accessions,stats[ind][['Accession']])
  }
  stats[match(accessions,Accession)]
  tmp<-acast(m.dat[acc%in%accessions],sample~acc,value.var = 'value')
  if(data.choice!=3) colnames(tmp)<-g.ref[match(colnames(tmp),Accession)][['Gene.Symbol']] else colnames(tmp)<-g.ref[match(p.ppd.features[match(colnames(tmp),Accession)][['Master.Protein.Accessions']],Accession)][['Gene.Symbol']]
  tmp[]<-apply(tmp,2,scale)
  Heatmap(t(tmp),cluster_columns = F,name='z score',column_title = paste(c(dat.names[data.choice],col.title),collapse = "\n"))
}
# 1 is total protein, 2 is phosphoprotein, 3 is phosphopeptide

hm.plot(data.choice = 1,genes = "PSMD12,IFRD1,PSMD14,YKT6,PHGDH,IDH1,PSMG1,LDHA,DHFR,GSR,PNP,HPRT1,PGK1,ALDOA,GAPDH,ENO1,GPI,EPRS1,ASNS,HMBS,LTA4H,HSPA5,SLC2A1,HSP90B1,GOT1,
        PFKL,WARS1,PSMA3,FKBP2,CALR,PSMB5,NMT1,RRM2,STIP1,CTH,HSPA4,PSMC2,PGM1,CCT6A,SLC1A4,NAMPT,PSMC4,MAP2K3,ME1,GCLC,PITPNB,GSK3B,ACLY,TPI1,ACTR3,ACTR2,UFM1,PPIA,TUBA4A,PSPH,
        SORD,GBE1,PRDX1,PDAP1,IDI1,PPA1,SLC1A5,TXNRD1,COPS5,ARPC5L,CACYBP,SDF2L1,TES,TBK1,SLC7A11,UCHL5,ATP6V1D,PSAT1",col.title='HALLMARK_MTORC1_SIGNALING')

hm.plot(data.choice = 2,genes="DNM1L,BIN1,SPAG9,FTH1,HMGB1,HSPA8,VCL,PGAM1,CTNNA1,PGM1,VCP,PDAP1,MTMR2,UBR4,CBARP,OSTF1,TMEM230,DBNL,MAGED2,PA2G4,DYNC1LI1",col.title = "GO_SECRETORY_VESICLE")
####### Please note filtering is sometimes necessary for phosphopeptides since there are multiple peptides under the same protein accession ID.
hm.plot(data.choice = 3,genes = "RIF1,KDM1A,SERBP1,PRRC2C,HSPA8,HMGA1,SSR1,MCM6",group=1,trend=2,col.title='MIDORIKAWA_AMPLIFIED_IN_LIVER_CANCER')
hm.plot(data.choice = 3,genes = "RIF1,KDM1A,SERBP1,PRRC2C,HSPA8,HMGA1,SSR1,MCM6",col.title='MIDORIKAWA_AMPLIFIED_IN_LIVER_CANCER')

