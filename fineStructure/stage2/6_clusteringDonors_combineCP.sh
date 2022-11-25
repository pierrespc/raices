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
# For computation reason we randomly selected a subset of 600 donors for stage1


###for parallelization of stage1 and stage2
step=20
totalInds=$(awk '{if($3==1)print $0}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids | wc -l)

cd $folder/fineStructure/Outputs/stage2/

####now stage 2: (describe by Ongaro et al. 2019): 
#We reconstruct each individual’s chromosomes as a series of genomic fragments inherited (copied) from a set of donor individuals,
# using the information on the allelic state of recipient and donors at each available position. 
#Briefly, we ‘painted’ the genomic profile of each donor as the combination of fragments received from other donor individuals.




if [ ! -e output.chunkcounts.out ]
then
	for i in {22..1}
	#for i in {22..1}
	do

		if [ ! -s stage2.chr$i.regionsquaredchunkcounts.out ]
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
				
				if [ ! -s stage2.chr$i.paral/$ind1.$ind2.regionsquaredchunkcounts.out ]
				then
					echo  stage2.chr$i.paral/$ind1.$ind2.regionsquaredchunkcounts.out does not exit!
					exit 
				fi
				let ' ind1 = ind1 + step '
				let ' ind2 = ind2 + step '
				let ' count = count + 1 '
				
			done
			jobS2comb=$(sbatch -J S2comb.$i -o chr$i.comb.o -e chr$i.comb.e --mem=4G --time=01:00:00 \
                                         --wrap "$commandCombine -d stage2.chr$i.paral/ -o stage2.chr$i")
		else
			echo stage2.chr$i already generated
		fi
	done
	jobS2combALL=$(sbatch -J S2comb.ALL -o ALL.comb.o -e ALL.comb.e --mem=4G --time=01:00:00 \
                                          --wrap "$commandCombine -l -m stage2.chr{1..22}")
else
	echo output.chunkcounts.out already generated
fi

