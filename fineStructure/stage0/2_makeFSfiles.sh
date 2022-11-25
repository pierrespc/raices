#!/bin/bash


SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/stage0/

commandConvert="$SinguExec $fsImage impute2chromopainter.pl"
#commandConvert="perl $folderScripts/2b_impute2chromopainter.pl"
commandFs="$SinguExec $fsImage fs"
module load singularity
module load java


rootRef=/pasteur/zeus/projets/p02/Hotpaleo/common_data/db/ref_datasets/Human/1000GP_Phase3/


folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/
cd $folder 

mkdir fineStructure
mkdir fineStructure/Inputs/


a=$(ls  shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr{1..22}_alignedRef.bim | wc -l)
echo $a
if [ $a -lt 22 ]
then
	echo "not 22 chr with alignedRef"
	answer="?"
	while [ $answer != "y" ] && [ $answer != "n" ]
	do 
		echo "copy $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr<CHR>.bim to $folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr<CHR>_alignedRef.bim????"
		echo "type y/n"
		read answer
	done
	if [ $answer == "n" ]
	then
		echo "well fix this then"
		exit
	else
		suff="_alignedRef"
		for chr in {1..22}
		do
			if [ ! -e shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_alignedRef.bim" ]
			then
				echo coyping for chr$chr
				cp shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr.bim shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_alignedRef.bim"
			else
				echo already there for chr$chr
			fi
		done
	fi
else
	suff="_alignedRef"
fi
echo "we have fixed bim files to :$folder/shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr<CHR>$suff.bim"
a=$(ls  shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr{1..22}$suff.bim | wc -l)
if [ $a != 22 ]
then
	echo "still not enough *bim files in $folder/shapeITphased/"
	exit
fi
a=$(ls  shapeITphased/withCM.Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr{1..22}$suff.bim | wc -l)
if [ $a != 22 ]
then
	Rscript $folderScripts/2a_putGenMap_inMapFile.R \
		T \
		shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr \
		$suff \
		T \
		T \
		$rootRef/genetic_map_chr \
		3 \
		shapeITphased/withCM.Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr
else
	echo cM ok!
fi
for chr in {1..22}
do
	echo $chr
	echo "start.pos recom.rate.perbp" > fineStructure/Inputs/chr$chr.recomb
	awk 'BEGIN{prevCM=0;prevBP=0} {rate=($3-prevCM)/($4-prevBP); print $4,rate; prevCM=$3; prevBP=$4}' shapeITphased/withCM.Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr$suff.bim >> fineStructure/Inputs/chr$chr.recomb
	if [  ! -s fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_alignedRef_phased.phase" ]
        then
                $commandConvert \
			shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_alignedRef_phased.haps" \
			fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_alignedRef_phased"
        else
                echo fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_alignedRef_phased.phase" already generated
        fi


done


####we have two steps in the analyze:
##step 1 clustering donors (everything except Puerto Madryn)
##step 2 painting Puerto Madryn individuals



awk '{if(NR>2){print $2,$1,1}}' shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr22_alignedRef_phased.sample > fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids



Rscript $folderScripts/2c_makeIDSperStep.R fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned.5.MeanByGroup.txt Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned.5.AncestryComponentByIndividual.txt fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1

 ###poor trick because soem phased files removed one Puerto Madryn individual (no problem for first steps of chromopainter because only for donors.
###I'll see how I can handle this for subsequent steps

#for chr in {22..1}
#do
#	awk '{if(NR>2){print $2,$1,1}}' shapeITphased/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_alignedRef_phased.sample" > fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_ALL.ids"
#	Rscript $folderScripts/2c_makeIDSperStep.R fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr"_ALL.ids" Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned.5.MeanByGroup.txt Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned.5.AncestryComponentByIndividual.txt fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$chr
#done
