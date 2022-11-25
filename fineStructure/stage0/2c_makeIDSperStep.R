#!/bin/Rscript
#

usage="<ALL ids> <AdmixtureMeanByGroups> <AdmixurePerInd> <prefOUT> \n with \n
<ALL ids>: the .ids for fineStructure with all individuals in fineStructure files
<AdmixtureMeanByGroups>: the file average Admixture estimations per pop
<AdmixurePerInd>: the file with the estimation per individual
<prefOUT>: the prefize of the output file\n";

params=commandArgs(trailingOnly=T)
if(length(params)!=4){
	stop(usage)
}

allIDS=params[1]
MeanByGroups=params[2]
PerInd=params[3]
prefOut=params[4]


###get components for Native American
MEAN<-read.table(MeanByGroups,sep="\t",header=T,stringsAsFactors=F)

native<-MEAN[grepl("Native",row.names(MEAN)),]
#print(native)
condition<-function(vec,th=0.5){
	return(sum(vec>th)>0)
}

native<-names(native)[apply(native,2,condition)]
#print(native)


###get list of individuals as donors : 
#	-all no american individuals 
#	-all Native american with Native American specific ancestry >0.95

###get list of recipients :
#	-all individuals from admixed populations
#	- -all Native american with Native American specific ancestry <=0.95

IND<-read.table(PerInd,sep="\t",header=T,stringsAsFactors=F)
if(length(native)==1){
	IND$Native<-IND[,native]
}else{
	IND$Native=apply(IND[,native],1,sum)
}



ids<-read.table(allIDS,stringsAsFactors=F,header=F)
IND=IND[ IND$Ind %in% ids$V1,]

recipients<-IND$Ind[ IND$Region %in% c("Argentina","Admixed American","AfroAmerican")]
recipients<-c(recipients, IND$Ind[(IND$Native<=0.95 & grepl("Native",IND$Region))])
donors<-IND$Ind[ ! IND$Region %in% c("Argentina","Admixed American","AfroAmerican") & !grepl("Native",IND$Region)]
donors=c(donors,IND$Ind[ (IND$Native>0.95 & grepl("Native",IND$Region))])

print(paste("regions of recipients (N=",length(recipients),")",sep=""))
print(table(IND$Region[IND$Ind %in% recipients]))
print(paste("regions of donors (N=",length(donors),")",sep=""))
print(table(IND$Region[IND$Ind %in% donors]))
if(sum(! IND$Ind %in% c(recipients,donors))>0){
         print("regions of unsassigned inds")
         print(table(IND$Region[! IND$Ind %in% c(recipients,donors)]))
         stop()
}       
 
if(sum(IND$Ind %in% recipients & IND$Ind %in% donors)>0){
         print("individuals as donor and recipiens")
         print(IND[ IND$Ind %in% recipients & IND$Ind %in% donors,])
         stop()
}       

ids$V3=NA
ids$V3[ ids$V1 %in% recipients]<-0
ids$V3[ ids$V1 %in% donors]<-1

###checks
if(sum(ids$V3==0)!=length(recipients) | sum(ids$V3==1)!=length(donors) | sum(is.na(ids$V3))>0){
	print(paste(sum(ids$V3==0),"vs",length(recipients)))
	print(paste(sum(ids$V3==1),"vs",length(donors)))
	print(sum(is.na(ids$V3)))
	stop("issues")
}

write.table(ids,paste(prefOut,"_STEP1.ids",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
write.table(ids[ids$V3==1,1],paste(prefOut,".initialDonorList",sep=""),col.names=F,row.names=F,sep="\t",quote=F)

table<-table(ids$V2[ids$V3==1])
propsTab<-data.frame(cbind(names(table),table/sum(table)))
donorsSubSet=sample(donors,600)
ids$V3<-0
ids$V3[ ids$V1 %in% donorsSubSet]<-1
write.table(ids,paste(prefOut,"_STEP1.subSet.ids",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
write.table(ids[ids$V3==1,1],paste(prefOut,".subSet.initialDonorList",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
table<-table(ids$V2[ids$V3==1])
temp<-data.frame(cbind(names(table),table/sum(table)))
propsTab<-merge(propsTab,temp,by=names(propsTab)[1],all=T)
propsTab[is.na(propsTab[,3]),3]<-0
print("correlation of sample size in full vs subset donor list")
cor=cor.test(as.numeric(propsTab[,2]),as.numeric(propsTab[,3]),method="spearman")
print(cor)
pdf(paste(prefOut,".subSet.initialDonorList.pdf",sep=""))
plot(as.numeric(propsTab[,2]),as.numeric(propsTab[,3]),
	xlab="Proportion of donors per population. Full set",
	ylab="Proportion of donors per population. Sub set",
	main=paste("Rho = ",round(cor$estimate,digits=4)))
dev.off()

