#/aplic/R/bin/R

usage = "Rscript putGenMap_inMapFile.R  <MANY_1> <mapFileIn> <suff> <bim> <MANY_2> <GenMapIn> <mapFileOut>
With
<MANY_1> T if you splitted the plink map files for each chrom or F if all in one unique file
<mapFileIn> the plink map file (if <MANY_1> is T) or prefixe (if <MANY_1> is F and then your map files must be <prefixe>1.map <prefixe>2.map etc...)
<suff>: suffix for map files
<bim> if map file is bim or not
<MANY_2> T if you splitted the genetic map files for each chrom or F if all in one unique file
<GenMapIn> the genetic  map file. as provided by Hapmap (if <MANY_2> is T) or prefixe (if <MANY_2> is F and then your map files must be <prefixe>1.txt <prefixe>2.txt etc...)
<colInGenetMap> colomn with genetic map in genetmap file
<mapFileOut> the plink map file to write (if <MANY_1> is T) or prefixe (if <MANY_1> is F and then your map files must be <prefixe>1.map <prefixe>2.map etc...)\n"


parameters <- commandArgs(trailingOnly=T)

print(paste(parameters,sep="\n"))
if (length(parameters) == 8){
	MANY_1 <- parameters[1]
	if(MANY_1 == "F"){
		MANY_1=FALSE
	}else{
		if(MANY_1 == "T"){
			MANY_1=TRUE
		}else{
			stop("<MANY_1> must be T or F check usage\n")
		}
	}	
	
	mapFileIn<-parameters[2]
    suff<-parameters[3]
    BIM<-parameters[4]
    if(BIM == "F"){
        BIM=FALSE
        suffMap="map"
    }else{
        if(BIM == "T"){
            BIM=TRUE
            suffMap="bim"
        }else{
            stop("<bim> must be T or F check usage\n")
        }
    }
	MANY_2 <- parameters[5]
	if(MANY_2 == "F"){
		MANY_2=FALSE
	}else{
		if(MANY_2 == "T"){
			MANY_2=TRUE
		}else{
			stop("<MANY_2> must be T or F check usage\n")
		}
	}	
	GenMapIn<-parameters[6]
    column=as.numeric(parameters[7])
    if(is.na(column)){
            stop("<colInGentMap> must be numeric")
    }
	mapFileOut<-parameters[8]

}else{
	stop(usage)
}




#####
if(!MANY_1){
	MapGW<-read.table(mapFileIn,stringsAsFactors=F,header=F)
	print(paste(mapFileIn,"has been read"))
	outGW=MapGW[0,]
    listCHR=unique(MapGW[,1])
}else{
    listCHR=c(1:22)
}
if(!MANY_2){
	GenMapGW<-read.table(GenMapIn,stringsAsFactors=F,header=T)[,c(1,column)]
	names(GenMapGW)<-c("pos","cM")
	print(paste(GenMapIn,"has been read"))
}	



for(chr in listCHR){
	print(paste("now processing chr",chr))	
	if(!MANY_1){
		MapCHR<-MapGW[MapGW$V1==chr,]
	}else{
		#MapCHR<-read.table(paste(mapFileIn,chr,".phased",suffMap,sep=""),stringsAsFactors=F,header=F)
		MapCHR<-read.table(paste(mapFileIn,chr,suff,".",suffMap,sep=""),stringsAsFactors=F,header=F)
		print(paste(mapFileIn,chr,suff,".",suffMap," has been read",sep=""))
	}
	if(dim(MapCHR)[1]==0){
		next
	}
    
    if(chr>22){
        if(!MANY_1){
            outGW<-rbind(outGW,MapCHR)
        }else{
            #write.table(MapCHR,paste(mapFileOut,chr,".phased",suffMap,sep=""),quote=F,row.names=F,col.names=F)
	write.table(MapCHR,paste(mapFileOut,chr,suff,".",suffMap,sep=""),quote=F,row.names=F,col.names=F)
            print(paste(mapFileOut,chr,suff,".",suffMap," has been written",sep=""))
        }
        next
    }
	if(!MANY_2){
		#GenMapCHR<-GenMapGW[GenMapGW$chr==chr,]
		GenMapCHR<-GenMapGW
	}else{

		GenMapCHR<-read.table(paste(GenMapIn,chr,"_combined_b37.txt",sep=""),stringsAsFactors=F,header=T)[,c(1,column)]
		names(GenMapCHR)<-c("pos","cM")
		print(paste(GenMapIn,chr,".txt has been read",sep=""))
	}
	
	print(head(MapCHR))
	print(head(GenMapCHR))		
	MapCHR$V3<-approx(x=GenMapCHR$pos,y=GenMapCHR$cM,xout=MapCHR$V4,rule=2)$y
	
	
	if(!MANY_1){	
		outGW<-rbind(outGW,MapCHR)
	}else{
		#write.table(MapCHR,paste(mapFileOut,chr,"phased",suffMap,sep=""),quote=F,row.names=F,col.names=F)
		write.table(MapCHR,paste(mapFileOut,chr,suff,".",suffMap,sep=""),quote=F,row.names=F,col.names=F)
		print(paste(mapFileOut,chr,suff,".",suffMap," has been written",sep=""))
	}
}

if(!MANY_1){
	write.table(outGW,mapFileOut,quote=F,row.names=F,col.names=F)	
	print(paste(mapFileOut,"has been written"))
}




