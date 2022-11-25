#!/bin/Rscript

require(plotrix)
require(scales)

system("mkdir /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES//AutoPerception/")
setwd("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES//AutoPerception/")
tab<-read.table("../AutoperceptionResults.txt",stringsAsFactors = F,header=T,sep="\t")
proportion<-read.table("../SummaryAncestryComponents_Admixture_SF.ByInd.tsv",stringsAsFactors = F,header=T,sep="\t")


tab<-merge(tab,proportion[,grepl("SFComb",names(proportion)) | grepl("AdmComb",names(proportion)) | names(proportion) =="Ind" ],by.x="Target",by.y="Ind")

names(tab)[names(tab)=="AdmComb.MiddleEastNorthAfricaCaucasus"]<-"AdmComb.Oriental"
names(tab)[names(tab)=="SFComb.EasternMediterranean"]<-"SFComb.Oriental"
for(meth in c("AdmComb","SFComb")){
  pdf(paste("Oriental_GeneticVsAutoperception_",meth,".pdf",sep=""),height=5,width=5)
  ###oriental
  ylim=c(0,max(tab[,paste(meth,".Oriental",sep="")]))
  #ylim=c(0,1)
  violin_plot(tab[tab$Oriental.self.perception==0,paste(meth,".Oriental",sep="")],at=1,xlim = c(0,3),ylim=ylim,axes=F,violin_width = 1,col="brown",ann=F)
  violin_plot(tab[tab$Oriental.self.perception==1,paste(meth,".Oriental",sep="")],at=2,add=T,violin_width = 1,col="brown")
  axis(2,at=seq(0,0.5,0.1))
  axis(1,at=c(0,1,2,3),labels = c(NA,"Without","With",NA),tick = F)
  wt<-wilcox.test(tab[,paste(meth,".Oriental",sep="")]~tab$Oriental.self.perception)$p.value
  #title(main=paste("Eastern Mediterranean Ancestry:\nGenetic (",meth,") VS Self-perceived\nWilcoxon-test P = ",round(wt,digits = 4),sep=""),
  #    xlab="Self-perceived Ancestry",
  #    ylab="Genetic ancestry")
  title(main=paste("Eastern Mediterranean Ancestry:\nWilcoxon-test P = ",round(wt,digits = 4),sep=""),
            xlab="Self-perceived Ancestry",
            ylab=paste("Genetic ancestry (",meth,")",sep=""))

  dev.off()

  ###European
  #ylim=range(tab$European.genomic.ancestry)
  pdf(paste("European_GeneticVsAutoperception_",meth,".pdf",sep=""),height=5,width=8)
  ylim=c(0,max(tab[,paste(meth,".Europe",sep="")]))
  #ylim=c(0,1)
  violin_plot(na.omit(tab[tab$European.self.perception==0.1,paste(meth,".Europe",sep="")]),at=1,xlim = c(0,6),ylim=ylim,axes=F,violin_width = 1,col="goldenrod",ann=F)
  x=1
  for(i in seq(0.3,0.9,0.2)){
    x=x+1
    violin_plot(na.omit(tab[tab$European.self.perception==i,paste(meth,".Europe",sep="")]),at=x,xlim = c(0,6),ylim=ylim,axes=F,violin_width = 1,col="goldenrod",ann=F,add=T)
  }
  axis(2)
  axis(1,c(0:6),labels = c(NA,paste(seq(0,80,20),"-",seq(20,100,20),"%",sep=""),NA),tick=F)
  sp<-cor.test(tab[,paste(meth,".Europe",sep="")],tab$European.self.perception,method = "spearman")
  #title(main=paste("European Ancestry:\nGenetic (",meth,") VS Self-perceived\nSpearman correlation rho = ",round(sp$estimate,digits=4)," (P = ",scientific(sp$p.value,digits = 3),")",sep=""),
  #    xlab="Self-perceived Ancestry",
  #    ylab="Genetic ancestry")
  title(main=paste("European Ancestry:\nSpearman correlation rho = ",round(sp$estimate,digits=4)," (P = ",scientific(sp$p.value,digits = 3),")",sep=""),
            xlab="Self-perceived Ancestry",
            ylab=paste("Genetic ancestry (",ifelse(meth=="AdmComb","Admixture","Sourcefind"),")",sep=""))
  points(c(0.5,1:5,5.5),c(0,seq(0.1,0.9,0.2),1),type = "l",lwd=2,lty=2)
  dev.off()
  ###Native American
  #ylim=range(tab$Native.american.genomic.ancestry)
  #ylim=c(0,1)
  pdf(paste("NativeAmerican_GeneticVsAutoperception_",meth,".pdf",sep=""),height=5,width=8)
  ylim=c(0,max(tab[,paste(meth,".NativeAmerican",sep="")]))
  violin_plot(na.omit(tab[tab$Native.american.self.percepcion==0.1,paste(meth,".NativeAmerican",sep="")]),at=1,xlim = c(0,6),ylim=ylim,axes=F,violin_width = 1,col="cadetblue",ann=F)
  x=1
  for(i in seq(0.3,0.9,0.2)){
    x=x+1
    violin_plot(na.omit(tab[tab$Native.american.self.percepcion==i,paste(meth,".NativeAmerican",sep="")]),at=x,xlim = c(0,6),ylim=ylim,axes=F,violin_width = 1,col="cadetblue",ann=F,add=T)
  }
  axis(2)
  axis(1,c(0:6),labels = c(NA,paste(seq(0,80,20),"-",seq(20,100,20),"%",sep=""),NA),tick=F)
  sp<-cor.test(tab[,paste(meth,".NativeAmerican",sep="")],tab$Native.american.self.percepcion,method = "spearman")
  #title(main=paste("Native American Ancestry:\nGenetic (",meth,") VS Self-perceived\nSpearman correlation rho = ",round(sp$estimate,digits=4)," (P = ",scientific(sp$p.value,digits = 3),")",sep=""),
  #    xlab="Self-perceived Ancestry",
  #    ylab="Genetic ancestry")
  title(main=paste("Native American Ancestry:\nSpearman correlation rho = ",round(sp$estimate,digits=4)," (P = ",scientific(sp$p.value,digits = 3),")",sep=""),
        xlab="Self-perceived Ancestry",
        ylab=paste("Genetic ancestry (",ifelse(meth=="AdmComb","Admixture","Sourcefind"),")",sep=""))
  points(c(0.5,1:5,5.5),c(0,seq(0.1,0.9,0.2),1),type = "l",lwd=2,lty=2)
  dev.off()
}

#####

for(meth in c("AdmComb","SFComb")){
  pdf(paste("Ttest_",meth,".pdf",sep=""),height=5,width=5)
  ylim=1.2*range(c(tab[,paste(meth,".Europe",sep="")]-tab$European.self.perception,
             tab[,paste(meth,".NativeAmerican",sep="")]-tab$Native.american.self.percepcion),na.rm = T)
  #ylim=c(-1,1)
  tstudent<-t.test(tab[,paste(meth,".Europe",sep="")],tab$European.self.perception,paired = T)$p.value
  violin_plot(na.omit(tab[,paste(meth,".Europe",sep="")]-tab$European.self.perception),col = "goldenrod",
            axes=F,
            ylab="Genetic - Self-perceived",
            at=1,violin_width = 1,
            xlim=c(0,3),
            ylim=ylim,
            main=paste("Genetic (",ifelse(meth=="AdmComb","Admixture","Sourcefind"),") VS Self-perceived Ancestry\nP-values : paired Student's t-test",sep=""))
  #text(x=1,y=ylim[1],labels=scientific(tstudent,digits=3))
  axis(2,at=seq(-1,1,0.2),las=2)
  axis(1,at=1,labels = paste("European\nP = ",scientific(tstudent,digits=3)),tick = F)
  #axis(1,at=1,labels = "European",tick=F)
  #axis(3,at=1,labels = paste("P = ",scientific(tstudent,digits=3)),tick = F)

  tstudent<-t.test(tab[,paste(meth,".NativeAmerican",sep="")],tab$Native.american.self.percepcion,paired = T)$p.value
  violin_plot(na.omit(tab[,paste(meth,".NativeAmerican",sep="")]-tab$Native.american.self.percepcion),col = "cadetblue",
            at=2,violin_width = 1,
            add=T)
#text(x=2,y=ylim[1],labels=scientific(tstudent,digits=3))

  axis(1,at=2,labels = paste("Native American\nP = ",scientific(tstudent,digits=3)),tick = F)
  #axis(1,at=2,labels = "Native American",tick=F)
  #axis(3,at=2,labels = paste("P = ",scientific(tstudent,digits=3)),tick = F)

  abline(h=0,lwd=2)
  dev.off()
}

