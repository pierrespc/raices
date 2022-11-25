#!/bin/bash

root=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

folderOUT=$root/Admixture_WW/
prefIN=$root/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling.pruned
numKmin=2
numKmax=5
numrep=1
RefInfo=$root/StartFilteredData/FinalColors_ALL.tsv
LabInfo=none
PopOrder=$root/StartFilteredData/PONG.ORDERPOPS_ALL



/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Admixture/6_plotPONG.sh $folderOUT \
		$prefIN \
		$numKmin \
		$numKmax \
		$numrep \
		$RefInfo \
		$LabInfo \
		$PopOrder
