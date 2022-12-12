#!/usr/bin/Rscript

require(colorspace)
library(RColorBrewer)
source("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/PLOTS/turboPalette.R")

folder="/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/"
pref="Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling.pruned"
listKNUM=11
OrderIND="None"
setwd(paste(folder,"/Admixture_WW/",sep=""))

orderRegions<-c("Subsaharan African","Middle Eastern - North African","Caucasus","European","East Asian","Central South Asian","Central Asian","Oceanian","Groenland","North Native American","Central Native American","South Native American","AfroAmerican","Admixed American","Argentina")
FocusReg<-"Argentina"
#FocusReg<-c("Ancient Central Chile","Central Chile","Patagonia","Ancient Patagonia","Tierra del Fuego","Ancient Tierra del Fuego")
fam<-read.table(paste(folder,"/StartFilteredData/",pref,".fam",sep=""),stringsAsFactors = F,header=F)
regions<-read.table(paste(folder,"/StartFilteredData/FinalColors_ALL.tsv",sep=""),stringsAsFactors=F,header=T,comment.char = "@",sep="\t")
regions<-regions[ regions$Ind %in% fam$V2,]

if(sum(!regions$Region %in% orderRegions)){
	print(paste(unique(regions$Region[!regions$Region %in% orderRegions]),collapse=" --- "))
	stop()
}
orderRegions<-orderRegions[ orderRegions %in% regions$Region]
adjustTextRight=300

ColorPalette<-turbo_colormap_data_HEX[round(seq(1,length(turbo_colormap_data_HEX),length.out=20))]
for(KNUM in listKNUM){
  #a<-read.table(paste("chr1.1KG.PopArg.pruned.",KNUM,".Q",sep=""),stringsAsFactors=F,header=F)
	a<-read.table(paste("K",KNUM,"/rep1/",pref,".rep1.K",KNUM,".Q",sep=""),stringsAsFactors=F,header=F)
#	outfile<-paste("chr1.1KG.PopArg.pruned.K",KNUM,".pdf",sep="")
	outfile<-paste(pref,".",KNUM,sep="")

  MyColors=turbo_colormap_data_HEX[round(seq(1,length(turbo_colormap_data_HEX),length.out=KNUM))]
  if(KNUM==11){
    MyColors=c(Yakut="grey40",EastAsia="grey30",
             Kalash="brown",Europe="orangered",
             Mozabite="goldenrod2",
             Oceania1="grey80",CentralAsia="grey60",
             Mbuti="grey10",Africakept="darkgreen",
             NatAm2="skyblue",
             NatAm="cadetblue")
  }
  numberK=dim(a)[2]
  numberInd=dim(a)[1]


  if(length(unique(MyColors))!=KNUM){
    print(MyColors[duplicated(MyColors)])
    stop("pb num col")
  }
  a$Ind<-fam$V2
  a$Population<-fam$V1
  
  a<-merge(a,regions,by=c("Ind","Population"))

  #a$Region[ is.na(a$Region)]<-a$

  numPops<-max(table(unique(a[,c("Population","Region")])$Region))
  numRegions<-1
  numInds<-max(table(unique(a[,c("Ind","Region")])$Region))
  if(dim(a)[1] != numberInd){
	  print(fam[ ! fam$V2 %in% a$Ind,])
	  stop(paste("your regions file and admixture output do not coincide: do not have the same number of Pops: a with",nrow(a),"and expected",numberInd))
  }

  out<-c()
  for(reg in orderRegions){
    print(reg)
    temp<-a[ a$Region == reg,]
    if(dim(temp)[1]==0){
      print(paste(reg,"skipped"))
      
      next
    }	
	  meanOverRegion<-vector(length=numberK)
	  names(meanOverRegion)<-paste("V",c(1:numberK),sep="")
	
	
	  meanOverPop<-data.frame(matrix(NA,length(unique(temp$Population)),numberK))
	  names(meanOverPop)<-paste("V",c(1:numberK),sep="")
	  rownames(meanOverPop)<-unique(temp$Population)
	
	  for(K in c(1:numberK)){
		  meanOverRegion[paste("V",K,sep="")]<-mean(temp[,paste("V",K,sep="")])
		  for(pop in unique(temp$Population)){
			  meanOverPop[pop,paste("V",K,sep="")]<-mean(temp[temp$Population==pop,paste("V",K,sep="")])
		  }
	  }
	
	  meanOverRegion<-meanOverRegion[order(as.numeric(meanOverRegion),decreasing=T)]
	  meanOverPop<-meanOverPop[,names(meanOverRegion)]
	  meanOverPop<-meanOverPop[order(meanOverPop[,1],decreasing=T),]
	  print(head(meanOverPop))
	  for(pop in rownames(meanOverPop)){
		  temp2<-a[ a$Population == pop,]
		  temp2<-temp2[order(temp2[,names(meanOverRegion)[1]],decreasing=T),]
		  out<-rbind(out,temp2)
	  }
  }

  if(OrderIND !="None"){
    orderInd<-read.table(OrderIND,stringsAsFactors=F,header=F)
    rownames(out)<-out$Ind
    out<-out[ orderInd$V1,]
    
  }
  

  #separ<-ceiling(numberInd/500)
  #separPop<-ceiling(numberInd/1000)
  separPop=1
  separ=1
  #pdf(paste(outfile,".pdf",sep=""), width=180,height=20)
  #par(mar=c(1, 2, 40, 2) + 0.1)

  ###plot Africa
  pdf(paste(outfile,"-A.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ITER = 0
  for(reg in "Subsaharan African"){
	  temp<-out[ out$Region == reg,]
	  Population=temp$Population[1]
	  ITER=ITER+1
	
	  Population2=Population
	  #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
	  axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
	  #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
	  for(ind in unique(temp$Ind)){
		  ybottom=0
		  if(temp$Population[ temp$Ind==ind] != Population){
			  Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
			  xleft=xleft+separPop
			  ITER=ITER+1
			
	  		#axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
			  Population2=Population
			  axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
		  	#axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
			  dimPrevRegion=dimPrevRegion+separPop
		  }
		  for(k in c(1:numberK)){
			  ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
			  rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
			  ybottom=ytop
		  }
		  xleft=xleft+1
	  }
	  xleft=xleft+separ
	  atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
	  dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }

  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line=4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "Subsaharan Africans",cex=10)
  dev.off()
  
  ###plot Mena 
  pdf(paste(outfile,"-B.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ITER = 0
  for(reg in c("Middle Eastern - North African","Caucasus")){
    temp<-out[ out$Region == reg,]
    Population=temp$Population[1]
    ITER=ITER+1
    
    Population2=Population
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    for(ind in unique(temp$Ind)){
      ybottom=0
      if(temp$Population[ temp$Ind==ind] != Population){
        Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
        xleft=xleft+separPop
        ITER=ITER+1
        
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        Population2=Population
        axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        dimPrevRegion=dimPrevRegion+separPop
      }
      for(k in c(1:numberK)){
        ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
        rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
        ybottom=ytop
      }
      xleft=xleft+1
    }
    xleft=xleft+separ
    atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
    dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }
  
  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line=4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "Middle Easterns\nNorth Africans\nand Caucasians",cex=10)
  dev.off()
  ###plot East Asia and Central Asia
  
  pdf(paste(outfile,"-C.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ITER = 0
  for(reg in c("East Asian","Central Asian")){
          temp<-out[ out$Region == reg,]
          Population=temp$Population[1]
          ITER=ITER+1

          Population2=Population
          #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
          axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
          #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
          for(ind in unique(temp$Ind)){
                  ybottom=0
                  if(temp$Population[ temp$Ind==ind] != Population){
                          Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
                          xleft=xleft+separPop
                          ITER=ITER+1

                        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
                          Population2=Population
                          axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
                        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
                          dimPrevRegion=dimPrevRegion+separPop
                  }
                  for(k in c(1:numberK)){
                          ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
                          rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
                          ybottom=ytop
                  }
                  xleft=xleft+1
          }
          xleft=xleft+separ
          atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
          dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }

  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line=4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "East Asians\nCentral Asians",cex=10)
  dev.off()

  ###plot Central Asians
  pdf(paste(outfile,"-D.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ITER = 0
  for(reg in c("Central South Asian")){
    temp<-out[ out$Region == reg,]
    Population=temp$Population[1]
    ITER=ITER+1
    
    Population2=Population
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    for(ind in unique(temp$Ind)){
      ybottom=0
      if(temp$Population[ temp$Ind==ind] != Population){
        Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
        xleft=xleft+separPop
        ITER=ITER+1
        
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        Population2=Population
        axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        dimPrevRegion=dimPrevRegion+separPop
      }
      for(k in c(1:numberK)){
        ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
        rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
        ybottom=ytop
      }
      xleft=xleft+1
    }
    xleft=xleft+separ
    atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
    dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }
  
  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line=4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "Central Asians\nCentral South Asians",cex=10)
  dev.off()

  ###plot Europe
  pdf(paste(outfile,"-E.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ITER = 0
  for(reg in "European"){
    temp<-out[ out$Region == reg,]
    Population=temp$Population[1]
    ITER=ITER+1
    
    Population2=Population
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    for(ind in unique(temp$Ind)){
      ybottom=0
      if(temp$Population[ temp$Ind==ind] != Population){
        Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
        xleft=xleft+separPop
        ITER=ITER+1
        
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        Population2=Population
        axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        dimPrevRegion=dimPrevRegion+separPop
      }
      for(k in c(1:numberK)){
        ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
        rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
        ybottom=ytop
      }
      xleft=xleft+1
    }
    xleft=xleft+separ
    atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
    dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }
  
  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line=4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "Europeans",cex=10)
  
  dev.off()
  
  
  ###plot Oceanian Native Americans
  pdf(paste(outfile,"-F.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ITER = 0
  for(reg in c("Oceanian","Groenland","North Native American","Central Native American","South Native American")){
    temp<-out[ out$Region == reg,]
    Population=temp$Population[1]
    ITER=ITER+1
    
    Population2=Population
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    for(ind in unique(temp$Ind)){
      ybottom=0
      if(temp$Population[ temp$Ind==ind] != Population){
        Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
        xleft=xleft+separPop
        ITER=ITER+1
        
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        Population2=Population
        axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        dimPrevRegion=dimPrevRegion+separPop
      }
      for(k in c(1:numberK)){
        ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
        rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
        ybottom=ytop
      }
      xleft=xleft+1
    }
    xleft=xleft+separ
    atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
    dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }
  
  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line=4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "Oceanians\nGroenland\nNative Americans",cex=10)
  dev.off()

  ###plot Afro 
  pdf(paste(outfile,"-G.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ITER = 0
  for(reg in c("AfroAmerican")){
          temp<-out[ out$Region == reg,]
          Population=temp$Population[1]
          ITER=ITER+1

          Population2=Population
          #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
          axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
          #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
          for(ind in unique(temp$Ind)){
                  ybottom=0
                  if(temp$Population[ temp$Ind==ind] != Population){
                          Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
                          xleft=xleft+separPop
                          ITER=ITER+1

                        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
                          Population2=Population
                          axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
                        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
                          dimPrevRegion=dimPrevRegion+separPop
                  }
                  for(k in c(1:numberK)){
                          ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
                          rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
                          ybottom=ytop
                  }
                  xleft=xleft+1
          }
          xleft=xleft+separ
          atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
          dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }

  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line=4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "Afro\nAmericans",cex=10)
  dev.off()
  ###plot Admixed 
  pdf(paste(outfile,"-H.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ITER = 0
  for(reg in c("Admixed American")){
    temp<-out[ out$Region == reg,]
    Population=temp$Population[1]
    ITER=ITER+1
    
    Population2=Population
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    for(ind in unique(temp$Ind)){
      ybottom=0
      if(temp$Population[ temp$Ind==ind] != Population){
        Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
        xleft=xleft+separPop
        ITER=ITER+1
        
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        Population2=Population
        axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=Population2,cex.axis=6,las=2,tick=T)
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population))),label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        dimPrevRegion=dimPrevRegion+separPop
      }
      for(k in c(1:numberK)){
        ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
        rect(xleft=xleft,xright=xleft+1,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
        ybottom=ytop
      }
      xleft=xleft+1
    }
    xleft=xleft+separ
    atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))+dimPrevRegion
    dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))+separ
  }
  
  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line=4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "Admixed\nAmericans",cex=10)
  dev.off()
  
  
  ####plot Puerto Madryn
  pdf(paste(outfile,"-I.pdf",sep=""), width=180,height=20)
  par(mar=c(1, 2, 40, 2) + 0.1)
  
  plot(0,0,"n",ylim=c(0,1),xlim=c(0,(numInds+separPop*(numPops-numRegions)+numRegions*(separ)))*1.2,main=paste("K =", numberK),ann=F,axes=F)
  title(paste("K =",numberK,sep=""))
  xleft=0
  atPop<-c()
  dimPrevRegion<- 0
  atRegion<-c()
  ZoomFac=2
  separ<-ceiling(numberInd/250)
  #separPop<-ceiling(numberInd/1000)
  separPop=1
  
  ITER = 0
  for(reg in orderRegions[ orderRegions %in% FocusReg]){
    temp<-out[ out$Region == reg,]
    Population=temp$Population[1]
    ITER=ITER+1
    
    Population2=Population
    axis(3,at=xleft+mean(c(0,sum(temp$Population==Population)))*ZoomFac,label=Population2,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population)))*ZoomFac,label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
    for(ind in unique(temp$Ind)){
      ybottom=0
      if(temp$Population[ temp$Ind==ind] != Population){
        Population=temp$Population[ temp$Ind==ind]
        #rect(xleft=xleft,xright=xleft+separPop,ybottom=0,ytop=1,col="black",border=NA)
        xleft=xleft+separPop*ZoomFac
        ITER=ITER+1
        Population2=Population
        axis(3,at=xleft+mean(c(0,sum(temp$Population==Population)))*ZoomFac,label=Population2,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        #axis(3,at=xleft+mean(c(0,sum(temp$Population==Population)))*ZoomFac,label=ITER,cex.axis=6,las=2,col.axis=unique(temp$Color[temp$Population==Population]),tick=T)
        dimPrevRegion=dimPrevRegion+separPop*ZoomFac
      }
      for(k in c(1:numberK)){
        ytop=ybottom+temp[ temp$Ind==ind,paste("V",k,sep="")]
        rect(xleft=xleft,xright=xleft+1*ZoomFac,ybottom=ybottom,ytop=ytop,col=MyColors[k],border=NA)
        ybottom=ytop
      }
      xleft=xleft+1*ZoomFac
    }
    xleft=xleft+separ*ZoomFac
    atRegion[reg]<-mean(c(0,sum(temp$Region==reg)))*ZoomFac+dimPrevRegion
    dimPrevRegion=atRegion[reg]+mean(c(0,sum(temp$Region==reg)))*ZoomFac+separ*ZoomFac
  }
  
  #axis(1,at=atRegion,label=names(atRegion),cex.axis=6,tick=F,line = 4)
  text(x=(numInds+separPop*(numPops-numRegions)+numRegions*(separ))*1.15,y=0.5,labels = "Present Study",cex=10)
  
  
  forLeg<-unique(out[,c("Population","Color")])
  #plot(0,0,"n",pch=0,ann=F,axes=F)
  #legend("topleft",legend=paste(c(1:dim(forLeg)[1]),": ",forLeg$Population,sep=""),ncol=10,cex=8,text.col=forLeg$Color)

  
  dev.off()
  
  meanByPop=data.frame(matrix(NA,0,KNUM))
  names(meanByPop)=paste(c(1:KNUM),sep="")
  iterP=0
  for(region in unique(out$Region)){
    iterP=iterP+1
    meanByPop[iterP,]<-apply(out[out$Region==region,paste("V",c(1:KNUM),sep="")],2,mean)
    rownames(meanByPop)[iterP]<-region
    for(pop in unique(out$Population[out$Region==region])){
      iterP=iterP+1
      meanByPop[iterP,]<-apply(out[out$Population==pop,paste("V",c(1:KNUM),sep="")],2,mean)
      rownames(meanByPop)[iterP]<-pop
    }
  }
  #write.table(meanByPop,paste("chr1.1KG.PopArg.pruned.",KNUM,".MeanByGroup.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  names(meanByPop)<-MyColors
  write.table(meanByPop,paste(outfile,".MeanByGroup.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)
  #  write.table(out,paste("chr1.1KG.PopArg.pruned.",KNUM,".AncestryComponentByIndividual.txt",sep=""),col.names=T,row.names=T,sep="\t",quote=F)  
  names(out)[names(out) %in% paste("V",c(1:KNUM),sep="")]<-MyColors
  write.table(out,paste(outfile,".AncestryComponentByIndividual.txt",sep=""),col.names=T,row.names=F,sep="\t",quote=F)
}


pdf(paste(pref,".CrossValidation.pdf",sep=""))
a<-read.table("CV.BESTruns",stringsAsFactors = F,header=T)
a<-a[ a$K >2,]
plot(CV~K,data=a,type="l",main="Cross-Validation Score",axes=F)
axis(2,at=seq(round(min(a$CV),digits = 4),round(max(a$CV),digits = 4),by = 2e-4))
axis(1,at=a$K)
points(CV~K,data=a,pch=4)
points(a$K[ a$CV==min(a$CV)],a$CV[ a$CV==min(a$CV)],col="red",pch=1,cex=3,lwd=2)
dev.off()


