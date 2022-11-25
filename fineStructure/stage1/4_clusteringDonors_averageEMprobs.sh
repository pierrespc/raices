#!/bin/bash



folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/stage1/

commandConvert="$SinguExec $fsImage impute2chromopainter.pl"
#commandConvert="perl $folderScripts/2b_impute2chromopainter.pl"
commandFs="$SinguExec $fsImage fs"
commandCombine="$SinguExec $fsImage chromocombine"

module load singularity
module load java


###Following: Ongaro et al. 2019. First cluster donor pops and then analyze admixed individuals as recipients.


###stage1
#estimated the mutation/emission and the switch rate parameters
#with ten steps of the Expectation–Maximization (E–M) algorithm 
#on a subset of chromosomes {4, 10, 15, 22}
#Using any individual except admixed individuals (including Native amerfican with <0.95 of native american ancestry) both as “donor” and “recipients.”
# For computation reason we randomly selected a subset of 600 donors for stage1


emStage1=10

mkdir $folder/fineStructure/Outputs/
cd $folder/fineStructure/Outputs/


###for parallelization of stage1 and stage2
step=20
totalInds=$(awk '{if($3==1)print $0}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids | wc -l)


#nR=$(cat $folder/fineStructure//Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.initialDonorList | wc -l)
#echo $nR
mkdir stage1
cd stage1 
if [ ! -e stage1.Combined ]
then
	for i in 3 7 10 18 22
	do

		if [ ! -s stage1.chr$i.EMprobs.out ]
		then
			
                         ind1=1
                         ind2=$step
                         count=1
                         while [ $ind1 -le $totalInds ]
                         do
				if [ $ind2 -gt $totalInds ]
				then
					ind2=$totalInds
				fi
				#echo here
				if [ ! -s stage1.chr$i.paral/$ind1.$ind2.EMprobs.out ]
	                        then
					echo stage1.chr$i.paral/$ind1.$ind2.EMprobs.out does not exist
					rm stage1.chr$i.EMprobs.out
					exit 
				fi
				cat stage1.chr$i.paral/$ind1.$ind2.EMprobs.out >> stage1.chr$i.EMprobs.out
                                let ' ind1 = ind1 + step '
                                let ' ind2 = ind2 + step '
                                let ' count = count + 1 '
				
			done
			
		else
			echo stage1.chr$i already generated
		fi
	done
	###average Ne and m across chromosomes and indviduals
	Rscript $folderScripts/4a_averageStage1.R $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr {22,18,10,7,3}
else
	echo stage1.Combined already generated
fi

mu=$(awk '{if(NR==2)print $2}' stage1.Combined)
Ne=$(awk '{if(NR==2)print $1}' stage1.Combined)

echo $mu $Ne

