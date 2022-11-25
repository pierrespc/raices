#!/bin/Rscript
require(ape)
require(dendextend)
require(stringr)
require(XML)

####function copied from FS scripts 
myapetodend<-function(ttree,simplify=FALSE,tol=0.1,factor=0.25){
  nodelab<-ttree$node.label
  ttree$node.label<-NULL
  htree<-my.as.hclust.phylo(ttree,tol)
  dend<-as.dendrogram(htree)
  #       if(!is.null(nodelab))dend<-setNodeLabels(dend,nodelab)
  dend<-dendrapply(dend,positiveheights)
  if(simplify) dend<-simpledend(dend)
  dend<-dendrapply(dend,flattenheights,factor=factor)
  if(!is.null(nodelab))dend<-setNodeLabels(dend,nodelab)
  dend
}
############
c19<-list(
  "EasternAfrica"="palegreen",
  "SouthernAfrica"="springgreen",
  "WesternAfrica2"="darkseagreen",
  "GuineanGulf"="forestgreen",
  "WesternAfrica"="darkolivegreen",
  
  "CentralSouthAsia"="brown",
  
  "Caucasus"="goldenrod2",
  "Levant"="goldenrod4",
  "MiddleEast"="lightgoldenrod",
  
  "NorthAfrica"="rosybrown",
  
  "Spain"="violetred1",
  "Italy"="violetred2",
  "Sardinia"="violetred3",
  "NorthEurope"="darkorange1",
  "WestEuropeCentralEurope"="darkorange3",
  "Basque"="orangered",
  
  "NorthAmerica"="skyblue",
  "SouthAmericanTropicalForests"="cadetblue",
  "CentralAmericaCentralAndes"="steelblue"
)
#################


setwd("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage4")

final<-read.tree("potentialFINAL_ManulTuning.nwk")
final$tip.label<-str_replace(final$tip.label,"_",":\n")

for(i in c(1:length(final$tip.label))){
  final$tip.label[i]<-str_replace_all(final$tip.label[i],"_"," ")
  final$tip.label[i]<-str_replace_all(final$tip.label[i],"-","; ")
  tmp<-strsplit(final$tip.label[i],"; ")[[1]]
  for(j in c(1:length(tmp))){
      tmp[j]<-paste(tmp[j],ifelse(j%%5==0,";\n","; "),sep="")
  }
  final$tip.label[i]<-paste(tmp,collapse = "")
}


final$tip.label<-str_replace(final$tip.label,"NorthWestEurope","WestEuropeCentralEurope")
final$tip.label<-str_replace(final$tip.label,"NorthAfricaLevant","MiddleEast")
final$tip.label<-str_replace(final$tip.label,"NorthAfricaCentralSouthAsia","NorthAfrica")
final$tip.label<-str_replace(final$tip.label,"; Ossetian","-Ossetian")
final$tip.label<-str_replace(final$tip.label,"; Jew","-Jew")
#pdf("potentialFINAL_ManulTuning.pdf")
tdend<-myapetodend(final,factor=1)

orderStuff<-c("EasternAfrica","WesternAfrica","WesternAfrica2","GuineanGulf","SouthernAfrica",
              "Sardinia","Italy","Spain","Basque",
              "WestEuropeCentralEurope","NorthEurope",
                      "CentralSouthAsia",
                      "NorthAfrica","MiddleEast","Levant","Caucasus",
                      
                      "NorthAmerica","CentralAmericaCentralAndes","SouthAmericanTropicalForests"
                      )
orderNum<-c()
for(p in orderStuff){
  w<-which(str_starts(final$tip.label,paste(p,":\n",sep="")))
  if(length(w)!=1){
    stop(paste(c(i,w),collapse = " "))
  }
  orderNum<-c(orderNum,w)
}
orderNum<-orderNum[c(length(orderNum):1)]
tdend<-rotate(tdend,orderNum)
tr <- as.phylo(tdend)

colList<-c()
for(i in tr$tip.label){
  tmp<-strsplit(i,split=":")[[1]][1]
  col=c19[[tmp]]
  if(length(col)!=1){
    stop(paste(c(tmp,col),collapse = " "))
  }
  colList<-c(colList,col)
}
pdf("FinalTree.pdf")
plot(tr,use.edge.length = FALSE, node.depth = 2,cex=0.4,tip.color =  colList)
dev.off()
#plot.dendrogram(tdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=0.6,las=2),edgePar=list(p.lwd=0,t.srt=0,t.off=0.3),yaxt="n",height=0.5,dLeaf=0.2)
#dev.off()
