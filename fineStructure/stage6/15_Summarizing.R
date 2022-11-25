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

numberIter=5


for(i in c("")){
  outTable<-data.frame(matrix(NA,0,length(c25)+1))
  names(outTable)<-c("Target",names(c25))
     
  pdf("BoxplotPerPops.pdf")                  
  #"stage6lessITER_4surrogates","stage6moreITER_4surrogates","stage6moreITER_8surrogates")){
  
  setwd(paste("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects//RAICES/fineStructure/Outputs/stage6/",i,sep=""))
  for(pop in c("PuertoMadryn","ACB","ASW","Colombia","Colombian","Maya","Mexican","Mixtec","Peru","PuertoRico")){
    outTablePOP<-data.frame(matrix(NA,0,length(c25)+1))
    names(outTablePOP)<-c("Target",names(c25))
    
    a<-read.table(paste(pop,"/stage6.",pop,".SF.saveout",sep=""),stringsAsFactors = F,header=T)
    a<-a[order(as.numeric(a$posterior.prob)),names(c25)]
    a<-data.frame(apply(a,2,as.numeric))
    names(a)<-names(c25)
    aHead<-tail(a,numberIter)
    
    outTablePOP[nrow(outTablePOP)+1,]<-c(pop,apply(aHead,2,mean))
    
    all<-read.table(paste(pop,"/stage6.",pop,".SFbyIND.saveout",sep=""),stringsAsFactors = F,header=T)
    for(ind in unique(all$target[all$target!="target"])){
      a<-all[all$target==ind,]
      a<-a[order(as.numeric(a$posterior.prob)),names(c25)]
      a<-data.frame(apply(a,2,as.numeric))
      
      names(a)<-names(c25)
      aHead<-tail(a,numberIter)
      
      outTablePOP[nrow(outTablePOP)+1,]<-c(paste(pop,ind,sep="___"),apply(aHead,2,mean))
    }
    outTable<-rbind(outTable,outTablePOP)
   
    outTablePOP[,-1]<-apply(outTablePOP[,-1],2,as.numeric)
    names(outTablePOP)[-1]<-sapply(names(outTablePOP)[-1], putSpace,USE.NAMES = F)
    b<-boxplot(outTablePOP[-1,-1],las=2,col=c25,main=paste(pop,"\n N = ",nrow(outTablePOP)-1))     
    points(x=c(1:(length(outTablePOP)-1)),y=outTablePOP[1,-1],pch=4,lwd=4)  
    means<-apply(outTablePOP[-1,-1],2,mean)
    points(x=c(1:(length(outTablePOP)-1)),y=means,pch=5,lwd=4)  
  }
  dev.off()
  outTable[grepl("___",outTable$Target),]
  write.table(outTable,"ProportionsPerIndividual.tsv",col.names = T,row.names = F,sep="\t",quote=F)
}
