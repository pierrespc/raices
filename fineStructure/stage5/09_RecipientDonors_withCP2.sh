#!/bin/bash



folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/stage5/

commandConvert="$SinguExec $fsImage impute2chromopainter.pl"
#commandConvert="perl $folderScripts/2b_impute2chromopainter.pl"
commandCP="$SinguExec $fsImage ChromoPainterv2 "
commandCombine="$SinguExec $fsImage chromocombine"

module load singularity
module load java


###stage5
#estimated the mutation/emission and the switch rate parameters


###for parallelization of stage1 and stage2
totalAdmInd=$(awk '{if($3==0)print $0}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids | wc -l)
totalDonorInd=$(awk '{if($3==1)print $0}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids | wc -l)
totalInd=$(awk '{if($3==1)print $0}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids | wc -l)
let ' diff = totalInd - totalDonorInd - totalAdmInd ' 
if [[ $diff != 0 ]]
then
	echo "pb num inds: $totalInd not $totalDonorInd + $totalAdmInd"
	exit
fi



mkdir $folder/fineStructure/Outputs/stage5/


cd $folder/fineStructure/Outputs/stage1/
mu=$(awk '{if(NR==2)print $2}' stage1.Combined)
Ne=$(awk '{if(NR==2)print $1}' stage1.Combined)

echo $mu $Ne

cd ../stage4/FinalGroups/
listCluster=$(ls *.DonorInds)

cd ../../stage5

rm stage5.*DonorList*
nTotalDonor=0
rm stage5.*ids*
iterClu=0
for clusterFile in $listCluster
do
	let ' iterClu = iterClu + 1 '
	cluster=$(echo $clusterFile | sed  s/.DonorInds//g)
	nCluster=$(wc -l ../stage4/FinalGroups/$clusterFile | awk '{print $1}')
	let ' nTotalDonor = nTotalDonor + nCluster '
	let ' nCluster = nCluster * 2 '
	echo -e $cluster"\tD"  >> stage5.DonorList
	echo -e pop$iterClu"\tD"  >> stage5.DonorList2
	awk -v c=$cluster -v OFS="\t" '{print $1,c,1}' ../stage4/FinalGroups/$clusterFile >> stage5.OnlyDonors.ids
	awk -v c=$iterClu -v OFS="\t" '{print $1,"pop"c,1}' ../stage4/FinalGroups/$clusterFile >> stage5.OnlyDonors.ids2
	
done

listAdmPop=$(awk -v OFS="\t" '{if($3==0)print $2}'  $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids | sort  | uniq)
#for popTest in $listAdmPop
for popTest in $listAdmPop
do
	mkdir $popTest
	cat stage5.OnlyDonors.ids > $popTest/stage5.$popTest.ids
	cat stage5.OnlyDonors.ids2 > $popTest/stage5.$popTest.ids2
	awk -v p=$popTest -v OFS="\t" '{if($2==p)print $1,$2,1}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids >> $popTest/stage5.$popTest.ids
	awk -v p=$popTest -v OFS="\t" '{if($2==p)print $1,$2,1}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids >> $popTest/stage5.$popTest.ids2
	cat  stage5.DonorList > $popTest/stage5.$popTest.DonorRecipList
	echo -e $popTest"\tR" >> $popTest/stage5.$popTest.DonorRecipList 
done
if [[ $nTotalDonor !=  $totalDonorInd ]]
then 
	echo "nTotalDonor found: $nTotalDonor vs $totalDonorInd expected"
	exit
fi
####now stage 2: (describe by Ongaro et al. 2019): 
#We reconstruct each individual’s chromosomes as a series of genomic fragments inherited (copied) from a set of donor individuals,
# using the information on the allelic state of recipient and donors at each available position. 
#Briefly, we ‘painted’ the genomic profile of each donor as the combination of fragments received from other donor individuals.

let ' totalInd = totalInd * 2 '
mkdir logs/
for popTest in $listAdmPop
#for popTest in PuertoRico
do
	if [[ ! -s $popTest/stage5.$popTest.regionsquaredchunkcounts.out ]]
	then
		#for i in 22
		for i in {21..1}
		do
			echo $popTest chr$i
			nL=$(wc -l $popTest/chr$i/stage5.$popTest.chr$i.regionsquaredchunkcounts.out | awk '{print $1}')
			let ' nL = nL - 1 '
			nIndTestPop=$(awk -v p=$popTest '{if($2==p)print $0}'  $popTest/stage5.$popTest.ids | wc -l)
			if [[ $nIndTestPop != $nL ]]
			#if [ ! -s $popTest/chr$i/stage5.$popTest.chr$i.regionsquaredchunkcounts.out ]
			then
				echo "$popTest/chr$i ($nL vs $nIndTestPop)"
				mkdir $popTest/chr$i
				numLinePhase=$(wc -l $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$i"_alignedRef_phased.phase" | awk '{print $1}')
				numLineNew=$(wc -l  $popTest/chr$i/stage5.$popTest.chr$i.phase | awk '{print $1}')
				numIndPhaseFile=$(wc -l $popTest/stage5.$popTest.ids | awk '{print $1}')
				let ' numLineExpect = ( numIndPhaseFile * 2 ) + 3 '
				if [[ $numLineNew == $numLineExpect ]]
				then
					echo  $popTest/chr$i/stage5.$popTest.chr$i.phase already generated
				else
					echo "generating  $popTest/chr$i/stage5.$popTest.chr$i.phase ($numLineNew != $numLineExpect)"
					 let ' nIndTestPop = totalDonorInd + nIndTestPop '
					###the four header lines NUMhaplotype Donors; NUM total ind; NUM snps; Positions
					if [[ $nIndTestPop != $numIndPhaseFile ]]
					then
						echo "issue pb ind in phase file $popTest/chr$i/stage5.$popTest.chr$i.phase ($nIndTestPop vs $numIndPhaseFile)"
						exit
					fi
					let ' nIndTestPop = nIndTestPop * 2 '
					
		                        echo $nIndTestPop > $popTest/chr$i/stage5.$popTest.chr$i.phase
			                sed -n 2,3p $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$i"_alignedRef_phased.phase" >> $popTest/chr$i/stage5.$popTest.chr$i.phase
					###prepare the reordered phase file
					while read line
					do
						a=($line)
						id=${a[0]}
						numFound=$(awk -v id=$id 'BEGIN{s=0} {if($1==id)s=s+1} END{print s}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids)
						if [[ $numFound != 1 ]]
						then
							echo "pb chr$i $popTest: $id found $numFound times"
							exit
						fi
						### get the line for phase
						line=$(awk -v id=$id '{if($1==id)print NR}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids)
						let ' line1 = ( line - 1 ) * 2 + 4 '
						let ' line2 = ( line - 1 ) * 2 + 5 '
						if [[ $line1 -lt 4 ]] || [[ $line1 -gt $numLinePhase ]] || [[ $line2 -lt 4 ]] || [[ $line2 -gt $numLinePhase ]]
						then
							echo "pb chr$i $popTest: $id line weird ($line-->$line1 $line2)"
							exit
						fi
						awk -v line1=$line1 -v line2=$line2 '{if(NR==line1 || NR==line2)print $0}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$i"_alignedRef_phased.phase" >> $popTest/chr$i/stage5.$popTest.chr$i.phase
				
					done < $popTest/stage5.$popTest.ids
				fi
				sbatch -J chr$i.$popTest -o logs/stage5.$popTest.chr$i.o -e logs/stage5.$popTest.chr$i.e --mem=4G --cpus-per-task=3 --time=72:00:00 -A metapaleo -p metapaleo \
					--wrap "$commandCP \
					-g $folder/fineStructure/Outputs/stage5/$popTest/chr$i/stage5.$popTest.chr$i.phase \
					-r $folder/fineStructure/Inputs/chr$i.recomb \
					-t $folder/fineStructure/Outputs/stage5/$popTest/stage5.$popTest.ids \
					-f $folder/fineStructure/Outputs/stage5/$popTest/stage5.$popTest.DonorRecipList 0 0 \
					-k 100 \
					-i 10 \
					-M $mu \
					-n $Ne \
					-ip \
					-o  $folder/fineStructure/Outputs/stage5/$popTest/chr$i/stage5.$popTest.chr$i"
			else
				echo  $popTest/chr$i/stage5.$popTest.chr$i already generated
			fi
		done
	else
		echo $popTest/stage5.$popTest already done
	fi
done

