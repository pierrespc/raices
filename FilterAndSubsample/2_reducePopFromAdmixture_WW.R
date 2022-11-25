#!/bin/bash

require(stringr)
root="/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/"
K=11

indTab<-read.table(paste(root,"Admixture_WW/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling.pruned.",K,".AncestryComponentByIndividual.txt",sep=""),
			header=T,sep="\t",stringsAsFactors=F,comment.char="@")     
ind=indTab[ indTab$Population=="PuertoMadryn",str_starts(names(indTab),"X.")]
group<-read.table(paste(root,"Admixture_WW/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling.pruned.",K,".MeanByGroup.txt",sep=""),
			header=T,sep="\t",stringsAsFactors=F,comment.char="@")

maxEurope<-names(group)[which(group["European",]==max(group["European",]))]

group<-group[ ! row.names(group) %in% indTab$Region,]
print(dim(ind))
maxInPM<-apply(ind,2,max)

print(maxInPM)
KinPM<-names(ind)[which(maxInPM>0.02)]
keepPops<-c()
maxPop<-c()
for(pop in row.names(group)){
	if(unique(indTab$Region[ indTab$Population==pop])!="European"){
		max<-max(group[pop,names(group)!=maxEurope])
	}else{
		max<-max(group[pop,])
	}
	maxPop<-c(maxPop,max)
}

for(K in KinPM){
	print(K)
	#print(summary(ind[K]))
	#print(row.names(group[ group[,K]>0.2,]))
	keepPops<-unique(c(keepPops,row.names(group[ group[,K]==maxPop,])))
}

print(KinPM)

keptList<-indTab[ indTab$Population %in% keepPops,c("Ind","Population","Region")]
print(unique(keptList[,c("Population","Region")]))
#keptList<-keptList[order(keptList$Region,kepList$Population)]

write.table(keptList[,c("Population","Ind")],paste(root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.KEEP",sep=""),
			col.names=F,row.names=F,quote=F)

system(paste("module load plink/1.90b6.16;
		plink 	--bfile ",root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling \\
			--keep ",root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.KEEP \\
			--make-bed \\
			--out ",root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED",sep=""))

system(paste("module load plink/1.90b6.16;
		plink 	--bfile ",root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED \\
			--indep-pairwise 50 5 0.5 \\
			--out ",root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruning",sep=""))

system(paste("module load plink/1.90b6.16;
		plink 	--bfile ",root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED \\
			--extract ",root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruning.prune.in \\
			--make-bed \\
			--out ",root,"/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned",sep=""))
