#!/bin/bash


root=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES
prefIN=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES_1stSubmission/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered


mkdir $root/StartFilteredData/
perl /pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Tools/King/getListIndToRemove_from_RelatedPairs.pl \
	$prefIN \
	$prefIN.fam \
	0.088 \
	0.0000001 \
	$root/StartFilteredData//Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered

	

module load plink/1.90b6.16
plink --bfile $prefIN \
	--remove $root/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.RemoveKin0.088 \
	--make-bed --out $root/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree

