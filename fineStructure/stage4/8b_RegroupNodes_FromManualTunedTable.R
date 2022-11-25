
require(stringr)
require(ape)
#library(dendextend)



###this step  has not been run in the cluster so we can check the outputs
source("/pasteur/zeus/projets/p02/Hotpaleo/pierre//Scripts/Tools/FineStructure_Rplot/FinestructureRcode/FinestructureLibrary.R")
source("/pasteur/zeus/projets/p02/Hotpaleo/pierre//Scripts/Tools/FineStructure_Rplot/FinestructureRcode/FinestructureDendrogram.R")


###SOME FUNCTIONS
######################################################################
######################################################################
######################################################################


####the home-made TVDcv function following this:
###To quantify the strength of differences between the inferred clusters
###we perform the following analyses.
###As noted above we can summarize the copying profiles of all the samples in a given cluster X
###to produce a characteristic copying vector x = (x1, x2,...,xn)
###; the average (across individuals in cluster X) proportion of each individual in cluster X's closest
###ancestry that is attributed to individuals from each of the clusters,
###Y = (Y1, Y2,...,Yn), where n is the number of inferred clusters.
TVD<-function(vec1,vec2){
  return(0.5*sum(abs(vec1-vec2)))
}

TVDmin<-function(vec1,vec2){
  return(0.5*min(abs(vec1-vec2)))
}

###function to get a subset of a string 
change<-function(string,split,pos){
  return(strsplit(string,split = split)[[1]][[pos]])
}

###calculate proportion
prop<-function(vec){
  return(vec/sum(vec))
}

###function to get all pairs of leafs connected by a direct ancestor in the tree
FindLeafPair<-function(treePhy,PairToFind){
  numTips=length(treePhy$tip.label)
  ###get Nodes leading to at least one leaf
  NodeToLeafs=as.data.frame(treePhy$edge[treePhy$edge[,2]<=numTips,])
  names(NodeToLeafs)<-c("Node","Leaf")
  if(max(table(NodeToLeafs$Node))>2){
    stop("ambiguous bifurcation")
  }
  ###get Nodes leading to exactly two leaf
  NodeToLeafs<-NodeToLeafs[ NodeToLeafs$Node %in% NodeToLeafs$Node[duplicated(NodeToLeafs$Node)],]
  pairFound=F
  for(nn in unique(NodeToLeafs$Node)){
    print(paste(treePhy$tip.label[NodeToLeafs$Leaf[NodeToLeafs$Node==nn]],collapse = " "))
    if(paste(treePhy$tip.label[NodeToLeafs$Leaf[NodeToLeafs$Node==nn]],collapse = " ")==PairToFind){
      if(pairFound){
        stop("found twice???")
      }
      pairFound=T
    }
  }
  return(pairFound)
}

### for fancier labels in plot
getLabelFUSION<-function(lab){
  tmp<-strsplit(lab,split = "\\)")[[1]][1]
  tmp<-strsplit(tmp,split = "\\(")[[1]][2]
  lab1<-strsplit(tmp,"and")[[1]][1]
  lab2<-strsplit(tmp,"and")[[1]][2]
  lab1<-str_replace_all(lab1," ","")
  lab2<-str_replace_all(lab2," ","")
  return(c(lab1,lab2))
}


######################################################################
######################################################################
######################################################################



setwd("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES//fineStructure/Outputs/stage4/")
treefile<-"../stage3/stage3.tree.Other.xml" ## finestructure tree file
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttreeSTART<-extractTree(treexml) ## extract the tree into ape's phylo format

tableTuned<-read.table("Choice_AllIter_OtherOrder.txt",stringsAsFactors = F,header=T,sep="\t")

###get starting clusters 
Colors<-read.table("../../../StartFilteredData/FinalColors_ALL.tsv",stringsAsFactors = F,header=T,sep="\t")
ttree<-ttreeSTART
startingClusters<-list()
n=0
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
keep<-c()
while(mapstate!=""){
  n=n+1
  cluster<-strsplit(mapstate, split=")",fixed=T)[[1]][1]
  cluster<-str_remove(string = cluster,pattern = "\\(")
  cluster<-strsplit(cluster,",")[[1]]
  pops<-Colors$Population[ Colors$Ind %in% cluster]
  pops<-pops[order(pops)]
  pops<-as.data.frame(table(pops))
  pops<-paste(pops$Freq,pops$pops,sep=" ")
  pops<-paste(pops,collapse = ";")
  ###if cluster have more than 1 tip, we merge the labels and remove all other tips
  if(length(cluster)>1){
      ttree<-drop.tip(ttree,tip=cluster[-1])
  }
  ttree$tip.label[ ttree$tip.label==cluster[1]]<-paste("X",n,sep="")
  
  ##rename Tip
  cluster<-paste(cluster,collapse = ",")
  startingClusters[[paste("X",n,sep="")]]<-list(pops=pops,ind=cluster,label=paste("X",n,sep=""))
  mapstate<-str_remove(mapstate,paste("\\(",cluster,"\\)",sep=""))
  
}





###BACKUP
startingClusters_BU<-startingClusters
ttree_BU<-ttree



#####THIS IS THE PROCESS FROM LEAFS TO ROOTS
startingClusters<-startingClusters_BU
ttree<-ttree_BU
#tdend<-myapetodend(ttree,factor=1)

keepOn=T
TH=1





for(i in c(1:nrow(tableTuned))){
  cluster1<-tableTuned$cluster1[i]
  cluster2<-tableTuned$cluster2[i]
  TVDcv<-tableTuned$TVDcv[i]
  Merge<-ifelse(tableTuned$NewLabel[i]=="NO","n","y")
  if(Merge=="y"){
    PairMerge<-paste(cluster1,cluster2)
    ##get the closest pairs
    PairFound<-FindLeafPair(ttree,PairMerge)
    if(!PairFound){
      stop(paste(tableTuned[i,],collapse = " "))
    }
  
  
    ###save the info for that iteration
    pop1=startingClusters[[cluster1]][["pops"]]
    pop2=startingClusters[[cluster2]][["pops"]]
    label<-paste("comparing",pop1," with ",pop2,"(",cluster1,"and",cluster2,")")
    print(label)
    listInd1<-strsplit(startingClusters[[cluster1]][["ind"]],",")[[1]]
    listInd2<-strsplit(startingClusters[[cluster2]][["ind"]],",")[[1]]
    region1<-paste(unique(Colors$Region[Colors$Ind%in%listInd1]),collapse=";")
    region2<-paste(unique(Colors$Region[Colors$Ind%in%listInd2]),collapse=";")
    label=paste(label,paste(region1,region2,sep=" vs "),sep="\n")
    
    label=paste(label,"\nTVD = ",TVDcv,sep="")
    ##"we merge both leafs and save it at cluster1")
    startingClusters[[cluster1]][["ind"]]<-paste(startingClusters[[cluster1]][["ind"]],startingClusters[[cluster2]][["ind"]],sep=",")
    listInd<-strsplit(startingClusters[[cluster1]][["ind"]],",")[[1]]
    pops<-Colors$Population[ Colors$Ind %in% listInd]
    pops<-pops[order(pops)]
    pops<-as.data.frame(table(pops))
    pops<-paste(pops$Freq,pops$pops,sep=" ")
    pops<-paste(pops,collapse = ";")
    startingClusters[[cluster1]][["pops"]]<-pops
    startingClusters[[cluster1]][["label"]]<-tableTuned$NewLabel[i]
    print(pops)
    ###"we delete cluster 2 from tree
    startingClusters<-startingClusters[ - which(names(startingClusters) == cluster2)]
    ttree<-drop.tip(ttree,cluster2)
    ttree$node.label<-NULL
   
  }
  
}

###get Merging Path for each cluster
system("mkdir TVDsequences")
system("mkdir FinalGroups")
TVDpath<-list()
for( tip in ttree$tip.label){
    ###retireve the last merging with that cluster
    MergingPath<-tableTuned[ (tableTuned$cluster1 ==tip |
                              tableTuned$cluster2==tip  ) 
                             & tableTuned$NewLabel!="NO",]
    #print(tip)
    #print(MergingPath)
    ##get the first STOP line
    stopLine<-tableTuned[ (tableTuned$cluster1 ==tip |
                             tableTuned$cluster2==tip  ) & tableTuned$NewLabel=="NO",][1,]
    stopLine$Color="red"
    if(tip=="X80"){
      nameLabel="Basque"
    }else{
      if(tip=="X58"){
        nameLabel="WesternAfrica2"
      }else{
        if(tip=="X91"){
          nameLabel="Sardinia"
        }else{
          MergingPath$Color="black"
          nameLabel=MergingPath$NewLabel[nrow(MergingPath)]
          if(nameLabel != startingClusters[[tip]]$label){
            stop(paste("pb label",tip))
          }
        }
      }
    }
    
    MergingPath<-rbind(MergingPath,stopLine)
    TVDpath[[ nameLabel ]]<-MergingPath
    pdf(paste("TVDsequences/",nameLabel,".pdf",sep=""))
    plot(MergingPath$TVDcv,col=MergingPath$Color,ylab="TVDcv",xlab="Iter",main=nameLabel)
    dev.off()
    write.table(MergingPath,paste("TVDsequences/",nameLabel,".tsv",sep=""),
                col.names=T,row.names = F,sep="\t",quote=F)
    write.table(strsplit(startingClusters[[tip]]$ind,split=",")[[1]],
                paste("FinalGroups/",nameLabel,".DonorInds",sep=""),
                col.names=F,row.names = F,sep="\t",quote=F)
    write.table(strsplit(startingClusters[[tip]]$pops,split=",")[[1]],
                paste("FinalGroups/",nameLabel,".DonorPops",sep=""),
                col.names=F,row.names = F,sep="\t",quote=F)
    
    
    
}




final<-ttree

print(final$tip.label)

numPerLine=5
for(lab in ttree$tip.label){
  tmp<-strsplit(startingClusters[[lab]]$pops,";")[[1]]
  out<-startingClusters[[lab]]$label
  for(i in seq(1,length(tmp),numPerLine)){
    start=i
    end=min(c(i+(numPerLine-1)),length(tmp))
    out<-paste(out,paste(tmp[c(start:end)],collapse = ";"),sep = "\n")
  }
  out<-str_replace(out,"X80\n","Basque\n")
  out<-str_replace(out,"X58\n","WesternAfrica2\n")
  out<-str_replace(out,"X91\n","Sardinia\n")
  final$tip.label[ final$tip.label==lab]<-out
  
}

print(final)

pdf("potentialFINAL_ManulTuning.pdf")
plot(final,use.edge.length = FALSE, node.depth = 2,cex=0.4)
dev.off()

write.tree(final,file="potentialFINAL_ManulTuning.nwk",append = F,tree.names = F)
###write the 
