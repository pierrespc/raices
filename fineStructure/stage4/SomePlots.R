#!/bin/Rscript

setwd("~/Documents/OtherCollaborations/RAICES/outputsFineStructure/")
source("~/Documents/PostDocPasteur/Scripts/Tools/FineStructure_Rplot/FinestructureRcode/FinestructureDendrogram.R")
source("~/Documents/PostDocPasteur/Scripts/Tools/FineStructure_Rplot/FinestructureRcode/FinestructureLibrary.R")

Colors<-read.table("../FinalColors.tsv",stringsAsFactors = F,header=T,sep="\t")

chunkfile<-"stage2/output.chunkcounts.out" ## chromopainter chunkcounts file
mcmcfile<-"stage3/stage3.mcmc.longer.xml" ## finestructure mcmc file
#treefile<-"stage3/stage3.tree.Other.xml.FSformat.xml" ## finestructure tree file
treefile<-"stage3/stage3.tree.Other.xml" ## finestructure tree file
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame




treexml<-xmlTreeParse(treefile) ## read the tree as xml format


ttree<-extractTree(treexml) ## extract the tree into ape's phylo format
## If you dont want to plot internal node labels (i.e. MCMC posterior assignment probabilities)
## now is a good time to remove them via:
ttree$node.label<-NULL
tdend<-myapetodend(ttree,factor=1) # convert to dendrogram format



## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering

mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations

popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only


popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointsComplete(popdend) # needed for obscure dendrogram reasons

popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend
popdendclear<-fixMidpointsComplete(popdendclear) # needed for obscure dendrogram reasons

pdf(file="stage3.tree.Other.SimplePopLabelDend.pdf",height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
plot.dendrogram(popdendclear,horiz=FALSE,nodePar=list(cex=0,lab.cex=0.6,las=2),edgePar=list(p.lwd=0,t.srt=0,t.off=0.3),yaxt="n",height=0.5,dLeaf=0.2)
dev.off()

pdf(file="stage3.tree.Other.RAWdend.pdf",height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
plot.dendrogram(tdend,horiz=FALSE,nodePar=list(cex=0,lab.cex=0.6,las=2,col=),edgePar=list(p.lwd=0,t.srt=0,t.off=0.3),yaxt="n",height=0.5,dLeaf=0.2)
dev.off()


row.names(Colors)<-Colors$Ind
Colors<-Colors[ttree$tip.label,]
ttreeOtherLab<-ttree
ttreeOtherLab$tip.label<-Colors$Population
  
pdf(file="stage3.tree.Other.RAWtreePOP.pdf",height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
plot(ttreeOtherLab,cex=0.1,direction = "upwards",tip.color = Colors$Color,use.edge.length = F)
dev.off()
pdf(file="stage3.tree.Other.RAWtreeIND.pdf",height=6,width=12)
par(mar=c(8,2,2,2),cex=0.8)
plot(ttree,cex=0.1,direction = "upwards",tip.color = Colors$Color,use.edge.length = F)
dev.off()

########################
## PAIRWISE COINCIDENCES

fullorder<-labels(tdend) # the order according to the tree
mcmcmatrixraw<-as.matrix(read.csv(meancoincidencefile,row.names=1)) # read in the pairwise coincidence file we created earlier
mcmcmatrix<-mcmcmatrixraw[fullorder,fullorder] 
mapstatematrix<-groupingAsMatrix(mapstatelist)[fullorder,fullorder] # map state for reference



