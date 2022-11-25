#!/bin/bash



folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/stage5/

commandConvert="$SinguExec $fsImage impute2chromopainter.pl"
#commandConvert="perl $folderScripts/2b_impute2chromopainter.pl"
formatCPmanual=T
if [[ $formatCPmanual == "F" ]]
then
	commandCP="$SinguExec $fsImage fs cp "
else
	commandCP="$SinguExec $fsImage chromopainter "
fi
commandCombine="$SinguExec $fsImage chromocombine"

module load singularity
module load java


###stage5
#estimated the mutation/emission and the switch rate parameters



cd $folder/fineStructure/Outputs/stage5/
mkdir logs/
cd ../stage4/FinalGroups/
listCluster=$(ls *.DonorInds)
cd ../../stage5/
listAdmPop=$(awk -v OFS="\t" '{if($3==0)print $2}'  $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids | sort  | uniq)

if [[ ! -s stage5.RecipientDonors.regionsquaredchunkcounts.out ]]
then
	listIN=""
	for clusterFile in $listCluster
	do
		echo $cluster
                cluster=$(echo $clusterFile | sed  s/.DonorInds//g)
        	tmp="$listIN $folder/fineStructure/Outputs/stage5/Donors/$cluster//stage5.$cluster "
		listIN=$tmp
	done
	for testPop in $listAdmPop
	do
		echo $testPop
		tmp="$listIN $folder/fineStructure/Outputs/stage5/$testPop/stage5.$testPop "
		listIN=$tmp
	done
        sbatch -J comb.ALL -o logs/stage5.ALL.comb.o -e logs/stage5.ALL.comb.e --mem=2G --cpus-per-task=1 --time=00:10:00 --qos fast -p common,dedicated \
               --wrap "$commandCombine -v -u \
               -o $folder/fineStructure/Outputs/stage5/stage5.RecipientDonors \
               $listIN"
fi


rm $folder/fineStructure/Outputs/stage5/stage5.RecipientDonors.ids
while read line
do
	a=($line)
	id=${a[0]}
	if [[ $id == "Recipient" ]] || [[ $id == "#Cfactor" ]]
	then
		continue
	fi
	don=$(grep -w $id stage5.OnlyDonors.ids | wc -l )
	echo $id
	if [[ $don == 1 ]]
	then
		grep -w $id stage5.OnlyDonors.ids >> $folder/fineStructure/Outputs/stage5/stage5.RecipientDonors.ids
	else
		if [[ $don == 0 ]]
		then
			rec=$(grep -w $id $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids | wc -l)
			if [[ $rec != 1 ]]
			then
				echo $id recipient found $rec
				exit
			else
				grep -w $id $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_ALL.ids >> $folder/fineStructure/Outputs/stage5/stage5.RecipientDonors.ids
			fi
		else
			echo $id donor found $don
		fi
	fi
done < $folder/fineStructure/Outputs/stage5/stage5.RecipientDonors.chunkcounts.out
