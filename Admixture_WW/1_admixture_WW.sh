#!/bin/bash

root=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/
outFolder=$root/Admixture_WW/
inFolder=$root/StartFilteredData/
mkdir $outFolder


pref=Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling
module load plink/1.90b6.16

#plink --bfile $inFolder/$pref --indep-pairwise 50 5 0.5 --out $inFolder/$pref.pruning
#plink --bfile $inFolder/$pref --extract $inFolder/$pref.pruning.prune.in --make-bed --out $inFolder/$pref.pruned

/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Admixture/1_runAdmixture.sh $outFolder/ $inFolder/$pref.pruned metapaleo 12 1
