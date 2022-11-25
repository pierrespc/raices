
require(stringr)
require(ape)
library(dendextend)

source("/pasteur/zeus/projets/p02/Hotpaleo/pierre//Scripts/Tools/FineStructure_Rplot/FinestructureRcode/FinestructureLibrary.R")
source("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Tools/FineStructure_Rplot/FinestructureRcode/FinestructureDendrogram.R")


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
GetLeafPairs<-function(treePhy,listPairsNotToTest){
  numTips=length(treePhy$tip.label)
  ###get Nodes leading to at least one leaf
  NodeToLeafs=as.data.frame(treePhy$edge[treePhy$edge[,2]<=numTips,])
  names(NodeToLeafs)<-c("Node","Leaf")
  if(max(table(NodeToLeafs$Node))>2){
    stop("ambiguous bifurcation")
  }
  ###get Nodes leading to exactly two leaf
  NodeToLeafs<-NodeToLeafs[ NodeToLeafs$Node %in% NodeToLeafs$Node[duplicated(NodeToLeafs$Node)],]
  pairs=list()
  for(nn in unique(NodeToLeafs$Node)){
    if(! paste(treePhy$tip.label[NodeToLeafs$Leaf[NodeToLeafs$Node==nn]],collapse = " ") %in% listPairsNotToTest){
      pairs[[length(pairs)+1]]<-treePhy$tip.label[NodeToLeafs$Leaf[NodeToLeafs$Node==nn]]
    }
  }
  return(pairs)
}


###function to get the Leaf pair with minimum TVD
GetLowestTVD<-function(LeafPairs,ChunckMatrixCluster){
  tvdPairs<-c()
  for(i in c(1:length(LeafPairs))){
    tvdPairs[i]<-TVD(ChunckMatrixCluster[LeafPairs[[i]][1],],ChunckMatrixCluster[LeafPairs[[i]][2],])
  }
  whMin<-which(tvdPairs==min(tvdPairs))
  if(length(whMin)>1){
    warning("different minimum for tvd")
    whMin=sample(whMin,1)
  }
  return(c(LeafPairs[[whMin]],tvdPairs[whMin]))
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



setwd("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs//")
chunkfile<-"stage2/output.chunkcounts.out" ## chromopainter chunkcounts file
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence
treefile<-"stage3/stage3.tree.Other.xml" ## finestructure tree file
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttreeSTART<-extractTree(treexml) ## extract the tree into ape's phylo format



###get starting clusters 
Colors<-read.table("../../StartFilteredData/FinalColors_ALL.tsv",stringsAsFactors = F,header=T,sep="\t")
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
  startingClusters[[paste("X",n,sep="")]]<-list(pops=pops,ind=cluster)
  mapstate<-str_remove(mapstate,paste("\\(",cluster,"\\)",sep=""))
  
}



###get Sum of chunck count per donor CLuster
chunkcountInd<-data.frame(matrix(NA,nrow(dataraw),length(startingClusters)))
row.names(chunkcountInd)<-row.names(dataraw)
for(ind in row.names(chunkcountInd)){
  for(cluster in names(startingClusters)){
    which<-which(row.names(dataraw) %in% strsplit(startingClusters[[cluster]][["ind"]],split=",")[[1]])
    chunkcountInd[ind,cluster]<-sum(dataraw[ind,which])    
  }
}


###BACKUP
chunkcountInd_BU<-chunkcountInd
startingClusters_BU<-startingClusters
ttree_BU<-ttree







#####THIS IS THE PROCESS FROM LEAFS TO ROOTS
chunkcountInd<-chunkcountInd_BU
startingClusters<-startingClusters_BU
ttree<-ttree_BU
#tdend<-myapetodend(ttree,factor=1)

keepOn=T
TH=1

listTREE<-list()
outITER<-c()
iter=0

listNotMerged<-c()

while(keepOn & length(startingClusters)>2){
  iter=iter+1
  print(paste("iter",iter))
  ###now calculate the proportion per individual of closest ancestry to other cluster
  chunkcountIndProp<-t(apply(chunkcountInd,1,prop))
  
  ###now get the avergae across proportion
  chunkcountCluster<-data.frame(matrix(NA,length(startingClusters),length(startingClusters)))
  names(chunkcountCluster)<-row.names(chunkcountCluster)<-names(startingClusters)
  for(cluster1 in names(startingClusters)){
    which1<-which(row.names(dataraw) %in% strsplit(startingClusters[[cluster1]][["ind"]],split=",")[[1]])
    for(cluster2 in names(startingClusters)){
      which2<-which(row.names(dataraw) %in% strsplit(startingClusters[[cluster2]][["ind"]],split=",")[[1]])
      chunkcountCluster[cluster1,cluster2]<-mean(chunkcountIndProp[which1,cluster2])
    }
  }
  
  ###now get pair of leafs to combine
  
  ##get the closest pairs
  ClosestPair<-GetLowestTVD(GetLeafPairs(ttree,listNotMerged),chunkcountCluster)
  cluster1<-ClosestPair[1]
  cluster2<-ClosestPair[2]
  TVDcv<-as.numeric(ClosestPair[3])
  
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
  diff=str_detect(region1,";") | str_detect(region2,";") | region1!=region2
  #if(! diff){
  #  Merge="y"
  #}else{
    Merge="y"
    while(! Merge %in% c("y","n")){
      Merge<-readline(label)
    }
  #}
  outITER<-rbind(outITER,cbind("mergeAccepted"=Merge,"diff"=diff,"cluster1"=cluster1,"cluster2"=cluster2,"pops1"=pop1,"pops2"=pop2,"Regions1"=region1,"Regions2"=region2,"TVDcv"=TVDcv))
  #print(TVDmin(chunkcountCluster[cluster1,],chunkcountCluster[cluster2,]))
  #if(TVDcv < TH){
  if(Merge=="y"){
    ##"we merge both leafs and save it at cluster1")
    startingClusters[[cluster1]][["ind"]]<-paste(startingClusters[[cluster1]][["ind"]],startingClusters[[cluster2]][["ind"]],sep=",")
    listInd<-strsplit(startingClusters[[cluster1]][["ind"]],",")[[1]]
    pops<-Colors$Population[ Colors$Ind %in% listInd]
    pops<-pops[order(pops)]
    pops<-as.data.frame(table(pops))
    pops<-paste(pops$Freq,pops$pops,sep=" ")
    pops<-paste(pops,collapse = ";")
    startingClusters[[cluster1]][["pops"]]<-pops
    print(pops)
    ###"we delete cluster 2 from tree
    startingClusters<-startingClusters[ - which(names(startingClusters) == cluster2)]
    ttree<-drop.tip(ttree,cluster2)
    ttree$node.label<-NULL
   
    ##"we actualize chunckCountPerInd for that cluster"
    chunkcountInd<-chunkcountInd[,names(chunkcountInd)!=cluster2]
    for(ind in row.names(chunkcountInd)){
        which<-which(row.names(dataraw) %in% listInd)
        chunkcountInd[ind,cluster1]<-sum(dataraw[ind,which])    
    }
  }else{
    #keepOn=F   
    listNotMerged<-c(listNotMerged,paste(cluster1,cluster2))
    print(listNotMerged)
  }
  listTREE[[iter]]<-list("TVDcv"=TVDcv,"tree"=ttree,"clusters"=startingClusters,"label"=label)
  
}
system("mkdir stage4")
write.table(outITER,"stage4/AllIter_OtherOrder.tsv",sep="\t",col.names = T,row.names = T,quote = F)


pdf("stage4/AllIter_OtherOrder.pdf")
ttree<-ttree_BU
clustersWILL<-getLabelFUSION(listTREE[[1]]$label)
plot(ttree,use.edge.length = FALSE, node.depth = 2,cex=0.4,
     tip.color = ifelse(ttree$tip.label %in% clustersWILL,"black","lightgrey"),
     main=listTREE[[1]]$label)

for(i in c(1:length(listTREE))){
  if(i!=length(listTREE)){
    clustersWILL<-getLabelFUSION(listTREE[[i+1]]$label)
  }else{
    clustersWILL=character(length = 0)
  }
  ttree<-listTREE[[i]]$tree
  plot(ttree,use.edge.length = FALSE, node.depth = 2,cex=0.4,
       tip.color = ifelse(ttree$tip.label %in% clustersWILL,"black","lightgrey"),
       main=ifelse(i!=length(listTREE),listTREE[[i+1]]$label,"Root Reached"))
}

dev.off()

final<-listTREE[[118]]


ttree<-final$tree
numPerLine=5
for(lab in ttree$tip.label){
  tmp<-strsplit(final$clusters[[lab]]$pops,";")[[1]]
  for(i in seq(1,length(tmp),numPerLine)){
    start=i
    end=min(c(i+(numPerLine-1)),length(tmp))
    if(i==1){
      out<-paste(tmp[c(start:end)],collapse = ";")
    }else{
      out<-paste(out,paste(tmp[c(start:end)],collapse = ";"),sep = "\n")
    }
    #print(out)
  }
  ttree$tip.label[ ttree$tip.label==lab]<-out
}


tdend<-myapetodend(ttree,factor=1)
pdf("stage4/potentialFINAL_OtherOrder.pdf")
plot(ttree,use.edge.length = FALSE, node.depth = 2,cex=0.4)
dev.off()
