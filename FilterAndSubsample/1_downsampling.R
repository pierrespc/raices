#!/bin/bash

root="/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/StartFilteredData/"

IND<-read.table(paste(root,"FinalColors_ALL.tsv",sep=""),stringsAsFactors=F,header=T,sep="\t")
listPOPnoSub<-IND$Population[ IND$Region %in% c("Argentina","Admixed American","AfroAmerican") | grepl("Native",IND$Region)]
TH=25
prefIN=paste(root,"/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree",sep="")
famIN<-read.table(paste(prefIN,".fam",sep=""),stringsAsFactors=F,header=F)
if(sum(!famIN$V1 %in% IND$Population)){
	print(unique(famIN$V1[!famIN$V1 %in% IND$Population ]))
	stop()
}
keep<-c()
remove<-c()
for(pop in unique(famIN$V1)){
	tmp=famIN[ famIN$V1==pop,]
	if(pop %in% listPOPnoSub | nrow(tmp)<=TH){
		print(paste(pop," kept all",nrow(tmp)))
		tmp2<-tmp
	}else{
		print(paste("downsampling ",pop,nrow(tmp)))
		tmp2=tmp[ sample(nrow(tmp),TH,replace=F),]
		remove<-rbind(remove,tmp[ ! tmp$V2 %in% tmp2$V2,])
	}
	keep<-rbind(keep,tmp2)
}

print(paste(nrow(remove),"removed +",nrow(keep),"kept =",nrow(famIN)))

write.table(keep,paste(root,"Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling.KEEP",sep=""),
		col.names=F,row.names=F,sep="\t",quote=F)
write.table(remove,paste(root,"Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling.REMOVE",sep=""),
		col.names=F,row.names=F,sep="\t",quote=F)

system(paste("module load plink/1.90b6.16;
		plink 	--bfile ",root,"Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree \\
			--remove ",root,"Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling.REMOVE \\
			--make-bed \\
			--out ",root,"Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling",sep=""))
