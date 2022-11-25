#!/bin/bash



folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/

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
# At the end, we did not did this: for computation reason we randomly selected a subset of 600 donors for stage1. We did on the whole set of individuals


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
count=0 
if [ ! -e stage1.Combined ]
then
	for i in 3 7 10 18
	#for i in 22
	do
		if [ ! -s stage1.chr$i.EMprobs.out ]
		then
			 mkdir stage1.chr$i.paral/
                         mkdir stage1.chr$i.paral/logs/
                         ind1=1
                         ind2=$step
                         while [ $ind1 -le $totalInds ]
                         do
				if [ $ind2 -gt $totalInds ]
				then
					ind2=$totalInds
				fi
				if [ ! -s stage1.chr$i.paral/$ind1.$ind2.EMprobs.out ]
	                        then
					echo stage1.chr$i.paral/logs/$ind1.$ind2
                                        jobS1paral=$(sbatch -J S1.$i.$ind1.$ind2 -o stage1.chr$i.paral/logs/$ind1.$ind2.o -e stage1.chr$i.paral/logs/$ind1.$ind2.e --mem=4G --time=12:00:00 --cpus-per-task=3 \
						--wrap "$commandFs cp \
                                                -g $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$i\"_alignedRef_phased.phase\" \
                                                -r $folder/fineStructure/Inputs/chr$i.recomb \
                                                -t $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids \
                                                -a $ind1 $ind2 \
						-i $emStage1 \
						-iM \
						-in \
						-o stage1.chr$i.paral/$ind1.$ind2")
				fi
                                let ' ind1 = ind1 + step '
                                let ' ind2 = ind2 + step '
                                let ' count = count + 1 '
				
			done
			
		else
			echo stage1.chr$i already generated
		fi
	done
	echo $count
else
	echo stage1.Combined already generated
fi

