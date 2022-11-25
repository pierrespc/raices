#!/bin/Rscript


params<-commandArgs(trailingOnly=T)
prefPhase=params[1]

listCHROM<-params[-1]

print(listCHROM)
print(prefPhase)


Ne<-c()
mu<-c()
wei<-c()
for(chr in listCHROM){
  print(chr)
  em<-read.csv(paste("stage1.chr",chr,".EMprobs.out",sep=""),stringsAsFactors = F,header=F,sep=" ")
  wh<-which(em$V1=="EMPAR")
  listInds=em$V3[wh]
  if(chr==listCHROM[1]){
    listIndsCHECK=listInds
  }else{
    if(length(listInds)!=length(listIndsCHECK)){
	print(length(listInds))
	print(length(listIndsCHECK))
      stop("pb inds across chr")
    }
    if(sum(listInds==listIndsCHECK)!=length(listIndsCHECK)){
      stop("pb ind corres across chr")
    }  
  }
  numem<-wh[2]-wh[1]-1
  if(length(unique(em[wh+numem,1]))!=1){
    stop("pb em read")
  }
  Ne<-cbind(Ne,as.numeric(em[wh+numem,4]))
  mu<-cbind(mu,as.numeric(em[wh+numem,5]))
  system(paste("sed -n 2p ",prefPhase,chr,"_alignedRef_phased.phase > chr",chr,".N",sep=""))
  wei<-c(wei,read.table(paste("chr",chr,".N",sep=""),header=F,stringsAsFactors = F)[1,1])
}


NEperind<-apply(Ne,2,mean)
MUperind<-apply(mu,2,mean)
NEmean<-weighted.mean(NEperind,wei)
MUmean<-weighted.mean(MUperind,wei)

write.table(data.frame(cbind(NEmean,MUmean),stringsAsFactors = F),"stage1.Combined",col.names = T,row.names = F,sep="\t",quote=F)


