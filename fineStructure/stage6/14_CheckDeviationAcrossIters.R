#!/bin/Rscript


require(stringr)
require(scales)
c25 <- c(
  "Caucasus"="#FF7F00",
  "NorthAfricaLevant"="#E31A1C",
  "Levant"="#FDBF6F",
  "CentralSouthAsia"="goldenrod",
  "NorthAfricaCentralSouthAsia"="brown",

  "Spain"="dodgerblue2",
  "Italy"="skyblue2",
  "NorthEurope"="blue1",
  "NorthWestEurope"="steelblue4",
  "Basque"="blue3",
  "Sardinia"="darkblue",
    
    

  "SouthernAfrica"="green1",
  "WesternAfrica2"="darkturquoise",
  "EasternAfrica"="darkgreen",
  "GuineanGulf"="palegreen2",
  "WesternAfrica"="green4",
  
  "NorthAmerica"="orchid3",
  "SouthAmericanTropicalForests"="deeppink1",
  "CentralAmericaCentralAndes"="darkviolet"
  )


putSpace<-function(string){
  for(let in LETTERS){
    string=str_replace_all(string,let,paste("\n",let,sep=""))
  }
  string=strsplit(string,split="")[[1]]
  return(paste(string[-1],collapse = ""))
}  
plot(c(1:length(c25)),pch=16,col=c25,cex=4)


colBlue=rgb(red=rep(0,1000),green=rep(0,1000),blue=rep(1,1000),alpha=seq(0,1,length=1000))
colRed=rgb(red=rep(1,1000),green=rep(0,1000),blue=rep(0,1000),alpha=seq(0,1,length=1000))

for(i in c("")){
  #"stage6lessITER_4surrogates","stage6moreITER_4surrogates","stage6moreITER_8surrogates")){
  
  setwd(paste("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage6/",i,sep=""))

  pdf("SummaryByPop.pdf")
  for(pop in c("PuertoMadryn","ACB","ASW","Colombia","Colombian","Maya","Mexican","Mixtec","Peru","PuertoRico")){
    a<-read.table(paste(pop,"/stage6.",pop,".SF.saveout",sep=""),stringsAsFactors = F,header=T)
    print(names(a)[! names(a) %in% names(c25)])
    a<-a[order(as.numeric(a$posterior.prob)),names(c25)]
    #a<-a[,names(c25)]
    #eq0=apply(a,1,check0)
    names(a)<-sapply(names(a), putSpace,USE.NAMES = F)
    aHead<-head(a,5)
    aTail<-tail(a,5)
    
    boxplot(a,las=2,cex.axis=0.3,main=paste(pop,"all iters"),col = c25)
    boxplot(aHead,las=2,cex.axis=0.3,main=paste(pop,"head"),col = c25)
    boxplot(aTail,las=2,cex.axis=0.3,main=paste(pop,"tail"),col = c25)
    outCorr=data.frame(matrix(1,length(a),length(a)))
    names(outCorr)<-row.names(outCorr)<-names(a)
    for(i in c(1:(length(a)-1))){
      for(j in c((i+1):(length(a)))){
        cor=cor.test(a[,i],a[,j])
        #plot(a[,i],a[,j],xlab=names(a)[i],ylab=names(a)[j],
        #     main=paste(names(a)[i],names(a)[j],
        #                "\nr2 = ",round(cor$estimate,digits=3),
        #                "( P = ",scientific(cor$p.value,digits=3),")"))
        outCorr[i,j]<-cor$estimate
        outCorr[j,i]<-ifelse(cor$estimate>0,
                           colBlue[round(cor$estimate*1000,digits=0)],
                           colRed[abs(round(cor$estimate*1000,digits=0))])
      
      }
    }
    plot(0,0,"n",xlim=c(0,(length(a)-1))+0.5,
       ylim=c(length(a),1)+0.5,axes=F,ann=F)
    axis(1,at=c(1:(length(a)-1)),labels = names(a)[-length(a)],cex.axis=0.3,las=2)
    axis(2,at=c(2:length(a)),labels = names(a)[-1],cex.axis=0.3,las=2)
    title(main=pop)
    for(i in c(1:(length(a)-1))){
      for(j in c((i+1):(length(a)))){
        rect(xleft=i-0.5,xright=i+0.5,
            ybottom = j-0.5,ytop = j+0.5,col=outCorr[j,i],border="black")
      }
    }
    rect(xleft=rep(length(a)-1,2001),xright=rep(length(a)-2,2001),
       ybottom=seq(1,10,length=2001),seq(1+9/2001,10+9/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
    text(x=rep(length(a)-0.9,9),y=seq(1,10,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.4)
    
    
  }
  dev.off()
  for(pop in c("PuertoMadryn","ACB","ASW","Colombia","Colombian","Maya","Mexican","Mixtec","Peru","PuertoRico")){
    all<-read.table(paste(pop,"/stage6.",pop,".SFbyIND.saveout",sep=""),stringsAsFactors = F,header=T)
    pdf(paste("SummaryByInd",pop,".pdf",sep=""))
    #a<-head(a[order(a$posterior.prob),names(c25)])
    for(ind in unique(all$target[all$target!="target"])){
      a<-all[all$target==ind,]
      a<-a[order(as.numeric(a$posterior.prob)),names(c25)]
      a<-data.frame(apply(a,2,as.numeric))
    
      #eq0=apply(a,1,check0)
      names(a)<-sapply(names(c25), putSpace,USE.NAMES = F)
      aHead<-head(a,5)
      aTail<-tail(a,5)
      
      boxplot(a,las=2,cex.axis=0.3,main=paste(ind,"all iters"),col = c25)
      boxplot(aHead,las=2,cex.axis=0.3,main=paste(ind,"head iters"),col = c25)
      boxplot(aTail,las=2,cex.axis=0.3,main=paste(ind,"tail iters"),col = c25)
      outCorr=data.frame(matrix(1,length(a),length(a)))
      names(outCorr)<-row.names(outCorr)<-names(a)
      for(i in c(1:(length(a)-1))){
        for(j in c((i+1):(length(a)))){
          cor=cor.test(a[,i],a[,j])
          #plot(a[,i],a[,j],xlab=names(a)[i],ylab=names(a)[j],
          #     main=paste(names(a)[i],names(a)[j],
          #                "\nr2 = ",round(cor$estimate,digits=3),
          #                "( P = ",scientific(cor$p.value,digits=3),")"))
          outCorr[i,j]<-cor$estimate
          outCorr[j,i]<-ifelse(cor$estimate>0,
                            colBlue[round(cor$estimate*1000,digits=0)],
                            colRed[abs(round(cor$estimate*1000,digits=0))])
      
        }
      }
      plot(0,0,"n",xlim=c(0,(length(a)-1))+0.5,
       ylim=c(length(a),1)+0.5,axes=F,ann=F)
      axis(1,at=c(1:(length(a)-1)),labels = names(a)[-length(a)],cex.axis=0.3,las=2)
      axis(2,at=c(2:length(a)),labels = names(a)[-1],cex.axis=0.3,las=2)
      title(main=ind)
      for(i in c(1:(length(a)-1))){
        for(j in c((i+1):(length(a)))){
          rect(xleft=i-0.5,xright=i+0.5,
            ybottom = j-0.5,ytop = j+0.5,col=outCorr[j,i],border="black")
        }
      }
      rect(xleft=rep(length(a)-1,2001),xright=rep(length(a)-2,2001),
       ybottom=seq(1,10,length=2001),seq(1+9/2001,10+9/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
      text(x=rep(length(a)-0.9,9),y=seq(1,10,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.4)
    }  
    dev.off()
  }
}
