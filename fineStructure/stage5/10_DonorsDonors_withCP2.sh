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

mkdir Donors/
for clusterFile in $listCluster
do
	let ' iterClu = iterClu + 1 '
        cluster=$(echo $clusterFile | sed  s/.DonorInds//g)
	mkdir Donors/$cluster/
        cat  stage5.DonorList > Donors/$cluster/stage5.$cluster.DonorRecipList
        echo -e $cluster"\tR" >> Donors/$cluster/stage5.$cluster.DonorRecipList 
done
####now stage 2: (describe by Ongaro et al. 2019): 
#We reconstruct each individual’s chromosomes as a series of genomic fragments inherited (copied) from a set of donor individuals,
# using the information on the allelic state of recipient and donors at each available position. 
#Briefly, we ‘painted’ the genomic profile of each donor as the combination of fragments received from other donor individuals.


###prepare the phase files
#for i in 22
for i in {21..1}
do
	mkdir Donors/chr$i/
	numLinePhase=$(wc -l $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$i"_alignedRef_phased.phase" | awk '{print $1}')
	numLineNew=$(wc -l  Donors/chr$i/stage5.Donors.chr$i.phase | awk '{print $1}')
	numIndPhaseFile=$(wc -l  stage5.OnlyDonors.ids | awk '{print $1}')
	let ' numLineExpect = (numIndPhaseFile * 2) + 3 '
	if [[ $numLineNew == $numLineExpect ]]
	then
		echo  "Donors/chr$i/stage5.Donors.chr$i.phase already generated ($numLineNew $numLineExpect)"
	else
		echo "generating  Donors/chr$i/stage5.Donors.chr$i.phase ($numLineNew != $numLineExpect)"
 	 	nIndTestPop=$nTotalDonor
		###the four header lines NUMhaplotype Donors; NUM total ind; NUM snps; Positions
		if [[ $nIndTestPop != $numIndPhaseFile ]]
		then
			echo "issue pb ind in phase file Donors/chr$i/stage5.Donors.chr$i.phase ($nIndTestPop vs $numIndPhaseFile)"
			exit
		fi
		let ' nIndTestPop = nIndTestPop * 2 '
	        echo $nIndTestPop > Donors/chr$i/stage5.Donors.chr$i.phase
	        sed -n 2,3p $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$i"_alignedRef_phased.phase" >> Donors/chr$i/stage5.Donors.chr$i.phase
		###prepare the reordered phase file
		while read line
		do
			a=($line)
			id=${a[0]}
			numFound=$(awk -v id=$id 'BEGIN{s=0} {if($1==id)s=s+1} END{print s}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids)
			if [[ $numFound != 1 ]]
			then
				echo "pb chr$i Donors: $id found $numFound times"
				exit
			fi
			### get the line for phase
			line=$(awk -v id=$id '{if($1==id)print NR}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids)
			let ' line1 = ( line - 1 ) * 2 + 4 '
			let ' line2 = ( line - 1 ) * 2 + 5 '
			if [[ $line1 -lt 4 ]] || [[ $line1 -gt $numLinePhase ]] || [[ $line2 -lt 4 ]] || [[ $line2 -gt $numLinePhase ]]
			then
				echo "pb chr$i Donors: $id line weird ($line-->$line1 $line2)"
				exit
			fi
			awk -v line1=$line1 -v line2=$line2 '{if(NR==line1 || NR==line2)print $0}' $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1.chr$i"_alignedRef_phased.phase" >> Donors/chr$i/stage5.Donors.chr$i.phase
				
		done < stage5.OnlyDonors.ids
	fi
	for clusterFile in $listCluster
	do	
		cluster=$(echo $clusterFile | sed  s/.DonorInds//g)
		nLineCounts1=$(wc -l Donors/$cluster/chr$i/stage5.$cluster.chr$i.chunkcounts.out | awk '{print $1}')
		nLineCounts2=$(wc -l Donors/$cluster/chr$i/stage5.$cluster.chr$i.chunklengths.out | awk '{print $1}')
		nLineCounts3=$(wc -l Donors/$cluster/chr$i/stage5.$cluster.chr$i.regionsquaredchunkcounts.out | awk '{print $1}')
		nLineExpected=$(awk -v p=$cluster '{if($2 == p)print $0}' $folder/fineStructure/Outputs/stage5/stage5.OnlyDonors.ids | wc -l ); 
		let ' nLineCounts1 = nLineCounts1 - 1 '
		let ' nLineCounts2 = nLineCounts2 - 1 '
		let ' nLineCounts3 = nLineCounts3 - 1 '
		if [[ $nLineCounts1 == $nLineExpected ]] && [[ $nLineCounts2 == $nLineExpected ]] && [[ $nLineCounts3 == $nLineExpected ]]
		then
			echo "Donors/$cluster/chr$i/stage5.$cluster.chr$i already made ($nLineCounts1 and $nLineCounts2 vs $nLineExpected)"
		else
			echo "making Donors/$cluster/chr$i/stage5.$cluster.chr$i ($nLineCounts1 and $nLineCounts2 and $nLineCounts3 vs $nLineExpected)"
			mkdir Donors/$cluster/chr$i
			#sbatch -J chr$i.$cluster -o logs/stage5.$cluster.chr$i.o -e logs/stage5.$cluster.chr$i.e --mem=4G --cpus-per-task=3 --time=100:00:00 -A metapaleo -p metapaleo \
			sbatch -J chr$i.$cluster -o logs/stage5.$cluster.chr$i.o -e logs/stage5.$cluster.chr$i.e --mem=10G --cpus-per-task=1 -p metapaleo -A metapaleo \
					--wrap "$commandCP \
					-g $folder/fineStructure/Outputs/stage5/Donors/chr$i/stage5.Donors.chr$i.phase \
					-r $folder/fineStructure/Inputs/chr$i.recomb \
					-t $folder/fineStructure/Outputs/stage5/stage5.OnlyDonors.ids \
					-f $folder/fineStructure/Outputs/stage5/Donors/$cluster/stage5.$cluster.DonorRecipList 0 0 \
					-k 100 \
					-ip \
					-M $mu \
					-n $Ne \
					-i 10 \
					-o  $folder/fineStructure/Outputs/stage5/Donors/$cluster/chr$i/stage5.$cluster.chr$i"
		fi
	done
done
