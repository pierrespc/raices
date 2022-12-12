#!/bin/Rscript
require(stringr)
setwd(paste("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects//RAICES/",sep=""))

legend<-read.table("StartFilteredData/BASEFinalColors_ALL.tsv",stringsAsFactors = F,header=T,sep="\t")
starting<-read.table("StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.fam",stringsAsFactors = F,header=F)
subsampled<-read.table("StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling.fam",stringsAsFactors = F,header=F)
admixture<-read.table("StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.fam",stringsAsFactors = F,header=F)

SFids<-read.table("fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids",stringsAsFactors = F,header=F)
donorIndList<-read.table("fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.initialDonorList",stringsAsFactors = F,header=F)

Receiptors<-SFids$V1[ ! SFids$V1 %in% donorIndList$V1]
legend<-legend[ legend$Ind %in% starting$V2,]
donorFiles<-list.files(path = "fineStructure/Outputs/stage4/FinalGroups/",pattern = ".DonorInds")
donorList<-list()
for(do in donorFiles){
  name<-str_remove(do,".DonorInds")
  name<-ifelse(name=="NorthAfricaLevant","MiddleEast",
               ifelse(name=="NorthAfricaCentralSouthAsia","NorthAfrica",name))
  
  donorList[[name]]<-read.table(paste("fineStructure/Outputs/stage4/FinalGroups/",do,sep=""),stringsAsFactors =F, header=F)$V1
}
donorList[["Receiptors"]]<-Receiptors

out<-data.frame(matrix(NA,0,6))
names(out)<-c("Population","Set","N_NoRelatedStarting","N_Subsampled","N_Admixture","N_Groups_SourceFind")
for(Population in unique(legend$Population)){
  N_NoRelatedStarting=sum(starting$V1==Population)
  if(N_NoRelatedStarting==0){
    stop(paste("1 --- ",Population))
  }
  Set=paste(unique(legend$set[ legend$Population==Population]),collapse = "/")
  if(Set=="CENPAT"){
    Set="Present Study"
  }
  N_Subsampled=sum(subsampled$V1==Population)
  N_Admixture=sum(admixture$V1==Population)
  
  N_Groups_SourceFind=NULL
  N_SourceFind=0
  for(group in names(donorList)){
    N=sum(legend$Ind %in% donorList[[group]] & legend$Population == Population)
    if(N==0){
      next
    }
    if(N_SourceFind==0){
      N_Groups_SourceFind=paste(N,group)
    }else{
      N_Groups_SourceFind=paste(N_Groups_SourceFind,"; ",N," ",group,sep="")
    }
    N_SourceFind=N_SourceFind+N
  }
  if(N_SourceFind!=N_Admixture){
    warning(paste("2 --- ",Population))
  }
  if(N_SourceFind==0){
    N_Groups_SourceFind=""
  }
  out<-rbind(out,cbind(Population,Set,N_NoRelatedStarting,N_Subsampled,N_Admixture,N_Groups_SourceFind))  
}

if(sum(as.numeric(out$N_NoRelatedStarting)) != nrow(starting)){
  stop("pb startring")
}
if(sum(as.numeric(out$N_Subsampled)) != nrow(subsampled)){
  stop("pb subsampled")
}
if(sum(as.numeric(out$N_Admixture)) != nrow(admixture)){
  stop("pb admixture")
}

write.table(out,"SummaryIndividualsIncluded.tsv",col.names=T,row.names=F,sep="\t",quote=F)

