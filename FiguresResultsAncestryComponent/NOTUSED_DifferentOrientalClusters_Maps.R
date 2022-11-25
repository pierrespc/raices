#!/bin/Rscript

#########
new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

#' Convenient functions
new_scale_fill <- function() {
  new_scale("fill")
}

new_scale_color <- function() {
  new_scale("colour")
}

new_scale_colour <- function() {
  new_scale("colour")
}
#############


system("mkdir /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/FiguresResultsAncestryComponent/")
setwd("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/FiguresResultsAncestryComponent/")
require(ggplot2)
require(scatterpie)
require(maps)
require(stringr)

world <- map_data('world')

infoREF<-read.table("../StartFilteredData/FinalColors_ALL.tsv",stringsAsFactors = F,header=T,sep="\t",comment.char = "@")
infoDUP<-unique(infoREF[,! names(infoREF) %in%  c("Ind","VCFid","Region","set","cex","Point")])

info<-infoDUP[0,]
for(i in unique(infoDUP$Population)){
  tmp<-infoDUP[infoDUP$Population==i,]
  if(nrow(tmp)==1){
    info<-rbind(info,tmp)
  }else{
    tmp$latitude[1]<-mean(tmp$latitude)
    tmp$longitude[1]<-mean(tmp$longitude)
    info<-rbind(info,tmp[1,])
  }
  
}



KperGroup<-read.table("../Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned.5.MeanByGroup.txt",
                      stringsAsFactors = F,header=T,sep="\t")

c5<-list("SubsaharanAfrica"="seagreen",
      "MiddleEastNorthAfricaCaucasus"="goldenrod",
      "Europe"="darkorange",
      "CentralSouthAsia"="brown",
      "NativeAmerican"="cadetblue")

c19<-list(
  "SouthernAfrica"="springgreen",
  "EasternAfrica"="palegreen",
  "WesternAfrica2"="darkseagreen",
  "GuineanGulf"="forestgreen",
  "WesternAfrica"="darkolivegreen",
  
  "CentralSouthAsia"="brown",
  
  "Caucasus"="goldenrod2",
  "Levant"="goldenrod4",
  "NorthAfricaLevant"="lightgoldenrod",
  
  "Spain"="violetred1",
  "Italy"="violetred2",
  "Sardinia"="violetred3",
  "NorthEurope"="darkorange1",
  "WestEuropeCentralEurope"="darkorange3",
  "Basque"="orangered",
  
  "NorthAfricaCentralSouthAsia"="rosybrown",
  
  "NorthAmerica"="skyblue",
  "SouthAmericanTropicalForests"="cadetblue",
  "CentralAmericaCentralAndes"="steelblue"
)
for(i in names(c5)){
  names(KperGroup)[names(KperGroup)==c5[[i]]]<-i
}
colors<-unlist(c5)
KperGroup$Population<-row.names(KperGroup)

KperGroup<-merge(KperGroup,info,by="Population")

scatters<-c()

for(pop in unique(KperGroup$Population)){
  scatters<-rbind(scatters,cbind(
    Population=pop,
    x=KperGroup$longitude[KperGroup$Population==pop],
    y=KperGroup$latitude[KperGroup$Population==pop],
    radius=sum(infoREF$Population==pop),
    KperGroup[KperGroup$Population==pop,names(c5)]))
}


xlim=c(-160,70)
ylim=c(-60,90)
ratYX=(ylim[2]-ylim[1])/(xlim[2]-xlim[1])

scatters$radius<-ifelse(scatters$Population =="PuertoMadryn",10,3)

pdf("PieChartsAdmixture.pdf",height = 10*ratYX,width = 10)
p <- ggplot(world, aes(long, lat)) +
  geom_scatterpie(aes(x=x, y=y,r=radius),
      data=scatters, cols=names(c5), alpha=1)+
  geom_map(map=world, aes(map_id=region), fill=NA, color="black")+
  scale_fill_manual(values = colors)+
  coord_cartesian(xlim=c(-160,70),ylim=c(-60,90))+
  #labs(title = "Genetic ancestry proportions (Admixture K= 8)\n") +
  theme_classic()

print(p)
dev.off()




###########NOW FOR FINESTRUCTURE
groups<-list(
  "SouthernAfrica"="SubsaharanAfrica",
  "EasternAfrica"="SubsaharanAfrica",
  "GuineanGulf"="SubsaharanAfrica",
  "WesternAfrica"="SubsaharanAfrica",
  "WesternAfrica2"="SubsaharanAfrica",
  
  "NorthAfricaLevant"="CentralSouthAsiaNorthAfricaLevant",
  "CentralSouthAsia"="CentralSouthAsiaNorthAfricaLevant",
  
  "NorthAfricaCentralSouthAsia"="CentralSouthAsiaNorthAfrica",
  
  "Caucasus"="CaucasusLevant",
  "Levant"="CaucasusLevant",
  
  "NorthEurope"="Europe",
  "Basque"="Europe",
  "WestEuropeCentralEurope"="Europe",
  "Sardinia"="Europe",
  "Italy"="Europe",
  "Spain"="Europe",
  
  "CentralAmericaCentralAndes"="NativeAmerican",
  "SouthAmericanTropicalForests"="NativeAmerican",
  "NorthAmerica"="NativeAmerican")

scattersFS<-scatters[,c(1:4)]
for(i in names(groups)){
  scattersFS<-cbind(scattersFS,0)
  names(scattersFS)[length(scattersFS)]<-i
}

###reading the group population
filesPops<-list.files(path="../fineStructure/Outputs//stage4/FinalGroups/",pattern = ".DonorPops")


sumIndsPerPop=list()
for(file in filesPops){
  tmp<-strsplit(readLines(paste("../fineStructure/Outputs//stage4/FinalGroups/",file,sep=""))[1],split = ";")[[1]]
  cluster<-str_replace(file,".DonorPops","")
  if(cluster=="NorthWestEurope"){
    cluster="WestEuropeCentralEurope"
  }
  for(cc in tmp){
    namePop<-strsplit(cc," ")[[1]][2]
    n<-as.numeric(strsplit(cc," ")[[1]][1])
    if(is.null(sumIndsPerPop[[namePop]])){
      sumIndsPerPop[[namePop]]=n
    }else{
      sumIndsPerPop[[namePop]]=n+sumIndsPerPop[[namePop]]
    }
    scattersFS[ scattersFS$Population==namePop,cluster]<-scattersFS[ scattersFS$Population==namePop,cluster]+n
    if(cluster=="NativeAmerican"){
      print(file)
      print(cc)
      print(n)
    }
  }  
}

for(pop in names(sumIndsPerPop)){
  scattersFS[ scattersFS$Population==pop,-c(1:4)]<-scattersFS[ scattersFS$Population==pop,-c(1:4)]/sumIndsPerPop[[pop]]
}

admixed<-read.table("../fineStructure/Outputs//stage6/ProportionsPerIndividual.tsv",stringsAsFactors =F,header=T) 

###check all donor read
print(scattersFS$Population[ ! scattersFS$Population %in% names(sumIndsPerPop)])
print(scattersFS$Population[ ! scattersFS$Population %in% c(names(sumIndsPerPop),admixed$Target)])

scattersFS<-scattersFS[ scattersFS$Population %in% c(names(sumIndsPerPop),"PuertoMadryn"),]

namePop="PuertoMadryn"
for(file in names(admixed)[-1]){
  if(file=="NorthWestEurope"){
    cluster="WestEuropeCentralEurope"
  }else{
    cluster=file
  }

  scattersFS[ scattersFS$Population==namePop,cluster]<-scattersFS[ scattersFS$Population==namePop,cluster]+
    mean(admixed[ grepl(paste(namePop,"___",sep=""),admixed$Target),file])
  
}

colors<-unlist(c19)
pdf("PieChartsSourceFinder.pdf",height = 10*ratYX,width = 10)
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="black")+
  geom_scatterpie(aes(x=x, y=y,r=radius),
                  data=scattersFS, cols=names(c19), alpha=1)+
  
  scale_fill_manual(values = colors)+
  coord_cartesian(xlim=c(-160,70),ylim=c(-60,90))+
  #$labs(title = "Proportion of Donor individuals assigned to each cluster\nSourcefind genetic ancestry proportions in Puerto Madryn") +
  theme_classic()

print(p)
dev.off()




#######################FINAL PLOT

####make data frame for pie charts Puerto Madryn
scatter19GROUP_PM<-scattersFS[ scattersFS$Population=="PuertoMadryn",]
scatter5GROUP_PM<-scatter19GROUP_PM[,c(1:4)]
scatter5GROUP_PMtrick<-scatter19GROUP_PM[,c(1:4)]
for(gn in unique(groups)){
    
    sum<-sum(scatter19GROUP_PM[,names(groups)[groups==gn] ])
    scatter5GROUP_PM<-cbind(scatter5GROUP_PM,sum)
    names(scatter5GROUP_PM)[length(scatter5GROUP_PM)]<-gn
    scatter5GROUP_PMtrick<-cbind(scatter5GROUP_PMtrick,sum)
    names(scatter5GROUP_PMtrick)[length(scatter5GROUP_PMtrick)]<-names(groups)[groups == gn][1]
}
for(gn in names(colors)){
  if(! gn %in% names(scatter5GROUP_PMtrick)){
    scatter5GROUP_PMtrick<-cbind(scatter5GROUP_PMtrick,0)
    names(scatter5GROUP_PMtrick)[length(scatter5GROUP_PMtrick)]<-gn
  }
}


scattersFS<-scattersFS[ scattersFS$Population!="PuertoMadryn",]


LONG=scatter19GROUP_PM$x[1]
LAT=scatter19GROUP_PM$y[1]


size1=3
size2=2.5


####prepare the annotation bar
colorsLEG<-data.frame(matrix(0,0,3))
names(colorsLEG)=c("pop","col","size")
lastGroup<-"NONE"
for(i in names(groups)){
  if(groups[[i]] != lastGroup){
    lastGroup=groups[[i]]
    colorsLEG<-rbind(colorsLEG,cbind(pop=paste(lastGroup," (",round(scatter5GROUP_PM[,lastGroup],digits=4)*100,"%)",sep=""),
                                     col=c19[[i]],size=size1))
  }
  colorsLEG<-rbind(colorsLEG,cbind(pop=paste(i," (",round(scatter19GROUP_PM[,i],digits=4)*100,"%)",sep=""),
                                   col=c19[[i]],size=size2))
}
colorsLEG$fill=ifelse(colorsLEG$size==3,NA,colorsLEG$col)
colorsLEG$border=ifelse(colorsLEG$size==3,NA,"black")
colorsLEG$col=ifelse(colorsLEG$size==3,colorsLEG$col,"black")
colorsLEG$hjust=ifelse(colorsLEG$size==3,"centre","inward")

putSpace<-function(string){
  for(let in LETTERS){
    string=str_replace_all(string,let,paste(" ",let,sep=""))
  }
  string=strsplit(string,split="")[[1]]
  return(paste(string[-1],collapse = ""))
}  
colorsLEG$pop<-sapply(colorsLEG$pop,putSpace,USE.NAMES = F)

xleft=-158
xright=-113
ytop=50
ybottom=-65
step=(ytop-ybottom)/(length(groups)+length(unique(groups))-1)

pdf("PieChartsSourceFinder_FINAL.pdf",height = 10*ratYX,width = 10)
p <- ggplot(world, aes(long, lat)) +
  geom_map(map=world, aes(map_id=region), fill=NA, color="black")+
  geom_scatterpie(aes(x=x, y=y,r=radius),
                  data=scattersFS, cols=names(colors), alpha=1)+
  geom_scatterpie(aes(x=x+40, y=y-10,r=radius),
                  data=scatter5GROUP_PMtrick, cols=names(colors), alpha=1)+
  geom_scatterpie(aes(x=x+40, y=y+15,r=radius),
                  data=scatter19GROUP_PM, cols=names(colors), alpha=1)+
  geom_point(aes(x=x,y=y),data=scatter5GROUP_PMtrick,shape=8,col="black",size=3,stroke=3)+
  scale_fill_manual(values = colors)+
  coord_cartesian(xlim=c(-160,70),ylim=c(-60,90))+
  labs(#title = paste("Proportion of \"donor\" individuals assigned to each cluster\t / \t Genetic ancestry proportion estimates in Puerto Madryn",sep=""),
       fill="19 Clusters") +
  annotate(geom="text",x=LONG+40,y=LAT+28,label="19 Clusters",size=5)+
  annotate(geom="text",x=LONG+40,y=LAT+3,label="6 Cluster Groups",size=5)+
  annotate(geom="rect",xmin=-190,xmax=xright,ymax=ytop+1*step,ymin=-90,fill='white',color="black")+
  annotate(geom="rect",xmin = xleft,xmax=xleft+step,ymin = seq(ytop,ybottom,-step),ymax = seq(ytop,ybottom,-step)+step*0.9,fill=colorsLEG$fill,color=colorsLEG$border)+
  annotate(geom="text",x=xleft+step*1.1,y=seq(ytop,ybottom,-step)+step/2,label=colorsLEG$pop,color=colorsLEG$col,size=as.numeric(colorsLEG$size),hjust = colorsLEG$hjust,vjust = "centre")+
  
  theme_classic()+
  theme(legend.position = "none")
  

print(p)

dev.off()
