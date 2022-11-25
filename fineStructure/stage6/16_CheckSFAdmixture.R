#!/bin/Rscript


require(stringr)
require(scales)

#####How we combine Admnixture: give for each component group the liost of population representing it
ListCombAdmix=list("SubsaharanAfrica"="Yoruba",
                   "CentralSouthAsia"="Kalash",
                   "MiddleEastNorthAfricaCaucasus"="Bedouin",
                   "Europe"="Finnish",
                   "NativeAmerican"="Karitiana")

####How we combine SF: 
###Final choice  of grouping: All "Oriental" clusters grouped together
listCombSF<-list(
  "SubsaharanAfrica"=c("SouthernAfrica","EasternAfrica","GuineanGulf","WesternAfrica","WesternAfrica2"),
  "EasternMediterranean"=c("NorthAfricaLevant","CentralSouthAsia",
                           "NorthAfricaCentralSouthAsia",
                           "Caucasus","Levant"),
  "Europe"=c("NorthEurope","Basque","WestEuropeCentralEurope","Sardinia","Italy","Spain"),
  "NativeAmerican"=c("CentralAmericaCentralAndes","SouthAmericanTropicalForests","NorthAmerica")
)





colorLISTS<-list()
colorLISTS[[paste("c",19)]] <- c(
  "SouthernAfrica"="springgreen",
  "EasternAfrica"="palegreen",
  "WesternAfrica2"="darkseagreen",
  "GuineanGulf"="forestgreen",
  "WesternAfrica"="darkolivegreen",
  
  "CentralSouthAsia"="brown",
  
  "Caucasus"="goldenrod2",
  "Levant"="goldenrod4",
  "NorthAfricaLevant"="lightgoldenrod",
  
  "Spain"="violetred1",
  "Italy"="violetred2",
  "Sardinia"="violetred3",
  "NorthEurope"="darkorange1",
  "WestEuropeCentralEurope"="darkorange2",
  "Basque"="darkorange3",
  
  "NorthAfricaCentralSouthAsia"="rosybrown",
  
  "NorthAmerica"="skyblue",
  "SouthAmericanTropicalForests"="cadetblue",
  "CentralAmericaCentralAndes"="steelblue"
)

colorLISTS[[paste("c",4)]]<-c("Europe"="darkorange",
      "EasternMediterranean"="brown",
      "NativeAmerican"="cadetblue",
      "SubsaharanAfrica"="seagreen")

colorLISTS[[paste("c",5)]]<-c(
  "SubsaharanAfrica"="seagreen",
  "MiddleEastNorthAfricaCaucasus"="goldenrod",
  "Europe"="darkorange",
  "CentralSouthAsia"="brown",
  "NativeAmerican"="cadetblue"
)



putSpace<-function(string){
  for(let in LETTERS){
    string=str_replace_all(string,let,paste("\n",let,sep=""))
  }
  string=strsplit(string,split="")[[1]]
  return(paste(string[-1],collapse = ""))
}  

putSpace2<-function(string){
  for(let in LETTERS){
    string=str_replace_all(string,let,paste(" ",let,sep=""))
  }
  
  string=strsplit(string,split="")[[1]]
  
  which<-which( string==" ")
  if(length(which)>1){
    which<-which[seq(1,length(which),2)]
    string[which]="\n"
  }
  return(paste(string[-1],collapse = ""))
}  


setwd(paste("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects//RAICES/fineStructure/Outputs/stage6/",sep=""))

SF<-read.table("ProportionsPerIndividual.tsv",header = T,stringsAsFactors =  F,sep="\t")
names(SF)[names(SF)=="NorthWestEurope"]="WestEuropeCentralEurope"
SF<-SF[,c("Target",names(colorLISTS[["c 19"]]))]
Admixture<-read.table("../../../Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned.5.AncestryComponentByIndividual.txt",stringsAsFactors = F,header=T,sep="\t")
Admixture$Ind<-paste(Admixture$Population,Admixture$Ind,sep="___")
if(sum(grepl("___",SF$Target))!= sum(SF$Target %in% Admixture$Ind)){
  stop("pb SF and Admixture Ind columns")
}

SF<-SF[grepl("___",SF$Target),]
Admixture<-Admixture[,! names(Admixture) %in% c("Region","latitude","longitude","set","VCFid","cex","Point","Color")]
names(Admixture)[c(1:2)]<-c("Ind","Population")
for(i in c(3:length(Admixture))){
  names(Admixture)[[i]]<-names(colorLISTS[[paste("c",length(Admixture)-2)]])[ which(colorLISTS[[paste("c",length(Admixture)-2)]]==names(Admixture)[i])] 
}
Admixture<-Admixture[,c("Ind","Population",names(colorLISTS[[paste("c",length(Admixture)-2)]]))]



pdf("CorrelationsAdmixture_RawSourceFinder.pdf")
#######comparing all CP/FS/SF clusters to all K in admixture...
colBlue=rgb(red=rep(0,1000),green=rep(0,1000),blue=rep(1,1000),alpha=seq(0,1,length=1000))
colRed=rgb(red=rep(1,1000),green=rep(0,1000),blue=rep(0,1000),alpha=seq(0,1,length=1000))

#for(pop in c(SF$Target[! grepl("___",SF$Target)],"All admixed")){
for(pop in c("PuertoMadryn","All admixed")){
  out<-data.frame(matrix(NA,ncol(Admixture)-2,ncol(SF)-1))
  names(out)<-names(SF)[-1]
  row.names(out)<-names(Admixture)[-c(1,2)]

  ###SF-Admixture
  plot(0,0,"n",xlab="Admixture components",ylab="SourceFinder Cluster",main=pop,
       xlim=c(-1.5,length(Admixture)-1),
       ylim=c(0,length(SF)+2),
       axes=F)
  y=0
  for(cluster in names(SF)[-1]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col=colorLISTS[[paste("c",19)]][cluster])
    text(x=x-1,y=y,labels = putSpace2(cluster),cex=0.3)
    for(K in names(Admixture)[-c(1,2)]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(SF[,c("Target",cluster)],Admixture[,c("Ind",K)],by.x="Target",by.y="Ind")
      }else{
        tmp<-merge(SF[grepl(paste(pop,"___",sep=""),SF$Target),c("Target",cluster)],Admixture[Admixture$Population==pop,c("Ind",K)],by.x="Target",by.y="Ind")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      #text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
      #                             ifelse(corP<1e-3,"3",
      #                                    ifelse(corP<1e-2,"2",
      #                                           ifelse(corP<1e-1,"1",""))))
      #)
      text(x=x,y=y,labels = scientific(corP,digits=3),cex=0.5)
      
      if(cluster==names(SF)[2]){
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(SF)+0.5,ytop=length(SF)+1.5,col=colorLISTS[[paste("c",5)]][K])
        text(x=x,y=length(SF)+1,labels = putSpace2(K),cex=0.3)
      }
    }
  }
  rect(xleft=rep(-1.5,2001),xright=rep(-0.5,2001),
       ybottom=seq(length(SF)+0.5,length(SF)+1.5,length=2001),
       ytop=seq(length(SF)+0.5+1/2001,length(SF)+1.5+1/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
  text(x=rep(-0.4,9),y=seq(length(SF)+0.5,length(SF)+1.5,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.4)
  
  
  
  #####Admixture-Admixture
  out<-data.frame(matrix(NA,ncol(Admixture)-2,ncol(Admixture)-2))
  names(out)<-names(Admixture)[-c(1,2)]
  row.names(out)<-names(Admixture)[-c(1,2)]
  plot(0,0,"n",xlab="Admixture components",ylab="Admixture components",main=paste("Correlations among components (Admixture)\n",pop,sep=""),
       xlim=c(-1.5,length(Admixture)-1),
       ylim=c(0,length(Admixture)+1),
       axes=F)
  y=0
  for(cluster in names(Admixture)[-c(1:2)]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col=colorLISTS[[paste("c",length(Admixture)-2)]][cluster])
    text(x=x-1,y=y,labels = putSpace2(cluster),cex=0.3)
    for(K in names(Admixture)[-c(1,2)]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(Admixture[,c("Ind",cluster)],Admixture[,c("Ind",K)],by="Ind")
      }else{
        tmp<-merge(Admixture[Admixture$Population==pop,c("Ind",cluster)],Admixture[Admixture$Population==pop,c("Ind",K)],by="Ind")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      #text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
      #                             ifelse(corP<1e-3,"3",
      #                                    ifelse(corP<1e-2,"2",
      #                                           ifelse(corP<1e-1,"1",""))))
      #)
      text(x=x,y=y,labels = scientific(corP,digits=3),cex=0.5)
      if(cluster==names(Admixture)[3]){
        print("ho")
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(Admixture)-0.5,ytop=length(Admixture)+0.5,col=colorLISTS[[paste("c",length(Admixture)-2)]][K])
        text(x=x,y=length(Admixture),labels = putSpace2(K),cex=0.3)
        
        
      }
    }
  }
  rect(xleft=rep(-1.5,2001),xright=rep(-0.5,2001),
       ybottom=seq(length(Admixture)-0.5,length(Admixture)+0.5,length=2001),
       ytop=seq(length(Admixture)-0.5+1/2001,length(Admixture)+0.5+1/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
  text(x=rep(-0.4,9),y=seq(length(Admixture)-0.5,length(Admixture)+0.5,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.2)
  
  #####SF-SF
  out<-data.frame(matrix(NA,ncol(SF)-1,ncol(SF)-1))
  names(out)<-names(SF)[-1]
  row.names(out)<-names(SF)[-1]
  plot(0,0,"n",xlab="SourceFinder Clusters",ylab="SourceFinder Clusters",main=paste("Correlations among components (SourceFind)\n",pop,sep=""),
       xlim=c(-1.5,length(SF)),
       ylim=c(0,length(SF)+2),
       axes=F)
  y=0
  for(cluster in names(SF)[-1]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col=colorLISTS[[paste("c",length(SF)-1)]][cluster])
    text(x=x-1,y=y,labels = putSpace2(cluster),cex=0.2)
    for(K in names(SF)[-1]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(SF[,c("Target",cluster)],SF[,c("Target",K)],by="Target")
      }else{
        tmp<-merge(SF[grepl(paste(pop,"___",sep=""),SF$Target),c("Target",cluster)],SF[grepl(paste(pop,"___",sep=""),SF$Target),c("Target",K)],by="Target")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      #text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
      #                             ifelse(corP<1e-3,"3",
      #                                    ifelse(corP<1e-2,"2",
      #                                           ifelse(corP<1e-1,"1",""))))
      #)
      text(x=x,y=y,labels = scientific(corP,digits=3),cex=0.25)
      if(cluster==names(SF)[2]){
        print("ho")
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(SF)+0.5,ytop=length(SF)+1.5,col=colorLISTS[[paste("c",length(SF)-1)]][K])
        text(x=x,y=length(SF)+1,labels = putSpace2(K),cex=0.2)
        
        
      }
    }
  }
  rect(xleft=rep(-1.5,2001),xright=rep(-0.5,2001),
       ybottom=seq(length(SF)+0.5,length(SF)+1.5,length=2001),
       ytop=seq(length(SF)+0.5+1/2001,length(SF)+1.5+1/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
  text(x=rep(-0.4,9),y=seq(length(SF)+0.5,length(SF)+1.5,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.2)
  
  
}

dev.off()



###Combining Components to better see the correlations between Admixture and SF
AdmixtureGROUP<-read.table("../../../Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned.5.MeanByGroup.txt",stringsAsFactors = F,header=T,sep="\t")
for(i in c(1:length(AdmixtureGROUP))){
  names(AdmixtureGROUP)[i]<-names(colorLISTS[[paste("c",length(AdmixtureGROUP))]])[ which(colorLISTS[[paste("c",length(AdmixtureGROUP))]]==names(AdmixtureGROUP)[i])] 
}



AdmixtureComb<-Admixture[,c("Ind","Population")]
for(i in names(ListCombAdmix)){
  comp=which(AdmixtureGROUP[ ListCombAdmix[[i]],]==max(AdmixtureGROUP[ ListCombAdmix[[i]],]))
  AdmixtureComb<-cbind(AdmixtureComb,Admixture[,names(AdmixtureGROUP)[ comp ]])
  names(AdmixtureComb)[length(AdmixtureComb)]<-i
}






SFcomb<-data.frame(matrix(NA,nrow(SF),1))
names(SFcomb)<-c("Target")
SFcomb$Target<-SF$Target

for(i in names(listCombSF)){
  print(i)
  if(length(listCombSF[[i]])>1){
    SFcomb<-cbind(SFcomb,apply(SF[,listCombSF[[i]] ],1,sum))
    
  }else{
    SFcomb<-cbind(SFcomb,SF[,listCombSF[[i]] ])
  }
  names(SFcomb)[length(SFcomb)]<-i
  
}



pdf("CorrelationsAdmixture_CombinatedSourceFinder.pdf")
for(pop in c("PuertoMadryn","All admixed")){
  #####SF-Admixture
  out<-data.frame(matrix(NA,ncol(AdmixtureComb)-2,ncol(SFcomb)-1))
  names(out)<-names(SFcomb)[-1]
  row.names(out)<-names(AdmixtureComb)[-c(1,2)]
  plot(0,0,"n",xlab="Admixture components (Combined)",ylab="SourceFind Clusters (Combined)",main=paste("Comparisons SourceFind Vs Admixture\n",pop,sep=""),
       xlim=c(-1.5,length(AdmixtureComb)-1),
       ylim=c(0,length(SFcomb)+2),
       axes=F)
  y=0
  for(cluster in names(SFcomb)[-1]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col= colorLISTS[[paste("c",length(SFcomb)-1)]][cluster])
    text(x=x-1,y=y,labels = putSpace2(cluster),cex=0.3)
    for(K in names(AdmixtureComb)[-c(1,2)]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(SFcomb[,c("Target",cluster)],AdmixtureComb[,c("Ind",K)],by.x="Target",by.y="Ind")
      }else{
        tmp<-merge(SFcomb[grepl(paste(pop,"___",sep=""),SFcomb$Target),c("Target",cluster)],AdmixtureComb[AdmixtureComb$Population==pop,c("Ind",K)],by.x="Target",by.y="Ind")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      #text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
      #                             ifelse(corP<1e-3,"3",
      #                                    ifelse(corP<1e-2,"2",
      #                                           ifelse(corP<1e-1,"1",""))))
      #)
      text(x=x,y=y,labels = scientific(corP,digits=3),cex=0.5)
      if(cluster==names(SFcomb)[2]){
        print("ho")
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(SFcomb)+0.5,ytop=length(SFcomb)+1.5,col=colorLISTS[[paste("c",length(AdmixtureComb)-2)]][K])
        text(x=x,y=length(SFcomb)+1,labels = putSpace2(K),cex=0.3)
        
        
      }
    }
  }
  rect(xleft=rep(-1.5,2001),xright=rep(-0.5,2001),
       ybottom=seq(length(SFcomb)+0.5,length(SFcomb)+1.5,length=2001),
       ytop=seq(length(SFcomb)+0.5+1/2001,length(SFcomb)+1.5+1/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
  text(x=rep(-0.4,9),y=seq(length(SFcomb)+0.5,length(SFcomb)+1.5,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.4)
  
  #####Admixture-Admixture
  out<-data.frame(matrix(NA,ncol(AdmixtureComb)-2,ncol(AdmixtureComb)-2))
  names(out)<-names(AdmixtureComb)[-c(1,2)]
  row.names(out)<-names(AdmixtureComb)[-c(1,2)]
  plot(0,0,"n",xlab="Admixture components (Combined)",ylab="Admixture components (Combined)",main=paste("Correlations among components (Admixture)\n",pop,sep=""),
       xlim=c(-1.5,length(AdmixtureComb)-1),
       ylim=c(0,length(AdmixtureComb)+1),
       axes=F)
  y=0
  for(cluster in names(AdmixtureComb)[-c(1:2)]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col=colorLISTS[[paste("c",length(AdmixtureComb)-2)]][cluster])
    text(x=x-1,y=y,labels = putSpace2(cluster),cex=0.3)
    for(K in names(AdmixtureComb)[-c(1,2)]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(AdmixtureComb[,c("Ind",cluster)],AdmixtureComb[,c("Ind",K)],by="Ind")
      }else{
        tmp<-merge(AdmixtureComb[AdmixtureComb$Population==pop,c("Ind",cluster)],AdmixtureComb[AdmixtureComb$Population==pop,c("Ind",K)],by="Ind")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      #text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
      #                             ifelse(corP<1e-3,"3",
      #                                    ifelse(corP<1e-2,"2",
      #                                           ifelse(corP<1e-1,"1",""))))
      #)
      text(x=x,y=y,labels = scientific(corP,digits=3),cex=0.5)
      if(cluster==names(SFcomb)[2]){
        print("ho")
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(AdmixtureComb)-0.5,ytop=length(AdmixtureComb)+0.5,col=colorLISTS[[paste("c",length(AdmixtureComb)-2)]][K])
        text(x=x,y=length(AdmixtureComb),labels = putSpace2(K),cex=0.3)
        
        
      }
    }
  }
  rect(xleft=rep(-1.5,2001),xright=rep(-0.5,2001),
       ybottom=seq(length(AdmixtureComb)-0.5,length(AdmixtureComb)+0.5,length=2001),
       ytop=seq(length(AdmixtureComb)-0.5+1/2001,length(AdmixtureComb)+0.5+1/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
  text(x=rep(-0.4,9),y=seq(length(AdmixtureComb)-0.5,length(AdmixtureComb)+0.5,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.4)

  #####SF-SF
  out<-data.frame(matrix(NA,ncol(SFcomb)-1,ncol(SFcomb)-1))
  names(out)<-names(SFcomb)[-1]
  row.names(out)<-names(SFcomb)[-1]
  plot(0,0,"n",xlab="SourceFinder Clusters (Combined)",ylab="SourceFinder Clusters (Combined)",main=paste("Correlations among components (SourceFind)\n",pop,sep=""),
       xlim=c(-1.5,length(SFcomb)),
       ylim=c(0,length(SFcomb)+2),
       axes=F)
  y=0
  for(cluster in names(SFcomb)[-1]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col=colorLISTS[[paste("c",length(SFcomb)-1)]][cluster])
    text(x=x-1,y=y,labels = putSpace2(cluster),cex=0.3)
    for(K in names(SFcomb)[-1]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(SFcomb[,c("Target",cluster)],SFcomb[,c("Target",K)],by="Target")
      }else{
        tmp<-merge(SFcomb[grepl(paste(pop,"___",sep=""),SFcomb$Target),c("Target",cluster)],SFcomb[grepl(paste(pop,"___",sep=""),SFcomb$Target),c("Target",K)],by="Target")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      #text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
      #                             ifelse(corP<1e-3,"3",
      #                                    ifelse(corP<1e-2,"2",
      #                                           ifelse(corP<1e-1,"1",""))))
      #)
      text(x=x,y=y,labels = scientific(corP,digits=3),cex=0.5)
      if(cluster==names(SFcomb)[2]){
        print("ho")
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(SFcomb)+0.5,ytop=length(SFcomb)+1.5,col=colorLISTS[[paste("c",length(SFcomb)-1)]][K])
        text(x=x,y=length(SFcomb)+1,labels = putSpace2(K),cex=0.3)
        
        
      }
    }
  }
  rect(xleft=rep(-1.5,2001),xright=rep(-0.5,2001),
       ybottom=seq(length(SFcomb)+0.5,length(SFcomb)+1.5,length=2001),
       ytop=seq(length(SFcomb)+0.5+1/2001,length(SFcomb)+1.5+1/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
  text(x=rep(-0.4,9),y=seq(length(SFcomb)+0.5,length(SFcomb)+1.5,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.4)
  
  
    
}

dev.off()

names(Admixture)[-c(1:2)]<-paste("Adm",names(Admixture)[-c(1:2)],sep=".")
names(SF)[-1]<-paste("SF",names(SF)[-1],sep=".")

names(AdmixtureComb)[-c(1:2)]<-paste("AdmComb",names(AdmixtureComb)[-c(1:2)],sep=".")
names(SFcomb)[-1]<-paste("SFComb",names(SFcomb)[-1],sep=".")
out<-merge(Admixture,AdmixtureComb,by=c("Ind","Population"))
out<-merge(out,SF,by.x="Ind",by.y="Target")
out<-merge(out,SFcomb,by.x="Ind",by.y="Target")

write.table(out,"../../../SummaryAncestryComponents_Admixture_SF.ByInd.tsv",col.names = T,row.names = F,sep="\t",quote=F)
outByPop<-out[0,c(1,1:ncol(out))]
names(outByPop)[c(1,2,3)]<-c("Population","N","Stat")

for(pop in unique(out$Population)){
  tmp<-out[ out$Population==pop,]
  print(paste("Pop",pop," N =",nrow(tmp)))
  tmp2<-data.frame(apply(tmp[,-c(1,2)],2,summary))
  outByPop<-rbind(outByPop,cbind(Population=pop,N=nrow(tmp),Stat=row.names(tmp2),tmp2))
}
write.table(outByPop,"../../../SummaryAncestryComponents_Admixture_SF.ByPOP.tsv",col.names = T,row.names = F,sep="\t",quote=F)
