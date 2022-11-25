#!/bin/bash

folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/stage6/

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

commandSourceFind="R < /pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Tools/sourcefindV2/sourcefindv2_PierreCall.R "
module load singularity
module load java


###stage5
#estimated the mutation/emission and the switch rate parameters



mkdir $folder/fineStructure/Outputs/stage6/
cd $folder/fineStructure/Outputs/stage6/
cd ../stage4/FinalGroups/
listCluster=$(ls *.DonorInds)
listIN=""
for clusterFile in $listCluster
do
	cluster=$(echo $clusterFile | sed  s/.DonorInds//g)
	tmp="$listIN$cluster "
	listIN=$tmp
done
cd ../../stage6/
listAdmPop=$(awk -v OFS="\t" '{if($3==0)print $2}'  $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids | sort  | uniq)
mkdir logs/
for popTest in $listAdmPop
do
	if [[ ! -e /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage6/$popTest/stage6.$popTest.SF.saveout ]]
	then
		mkdir $popTest
		echo "self.copy.ind: 0
num.surrogates: 8
exp.num.surrogates: 4
input.file.ids: /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage5/stage5.RecipientDonors.ids
input.file.copyvectors: /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage5/stage5.RecipientDonors.chunklengths.out
save.file: /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage6/$popTest/stage6.$popTest.SF.saveout
copyvector.popnames: $listIN
surrogate.popnames: $listIN
target.popnames: $popTest
num.slots: 100
num.iterations: 400000
num.burnin: 100000
num.thin: 1000" > $popTest/stage6.SF.$popTest.params
		sbatch -J SF.$popTest -o logs/stage6.SF.$popTest.o -e logs/stage6.SF.$popTest.e --mem=6GB --time=10:00:00 \
			--wrap "Rscript /pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Tools/sourcefindV2/sourcefindv2_PierreCall.R $folder/fineStructure/Outputs/stage6/$popTest/stage6.SF.$popTest.params"
	else
		echo /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage6/$popTestn/stage6.$popTest.SF.saveout already done
	fi
	if [[ ! -e /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage6/$popTest/stage6.$popTest.SFbyIND.saveout ]]
	then
		indTest=$(awk -v p=$popTest '{if($2==p)print $1}' $folder/fineStructure/Outputs/stage5/stage5.RecipientDonors.ids)
		listINindTest=""
		for ind in $indTest
		do
			tmp="$listINindTest$ind "
			listINindTest=$tmp
		done
		awk -v OFS="\t" -v p=$popTest '{if($2==p){$2=$1}; print $0}' $folder/fineStructure/Outputs/stage5/stage5.RecipientDonors.ids > $popTest/stage6.$popTest.SFbyIND.ids
                echo "self.copy.ind: 0
num.surrogates: 8
exp.num.surrogates: 4
input.file.ids: /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage6/$popTest/stage6.$popTest.SFbyIND.ids
input.file.copyvectors: /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage5/stage5.RecipientDonors.chunklengths.out
save.file: /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage6/$popTest/stage6.$popTest.SFbyIND.saveout
copyvector.popnames: $listIN
surrogate.popnames: $listIN
target.popnames: $listINindTest
num.slots: 100
num.iterations: 400000
num.burnin: 100000
num.thin: 1000" > $popTest/stage6.SFbyIND.$popTest.params
		sbatch -J SFi.$popTest -o logs/stage6.SFbyIND.$popTest.o -e logs/stage6.SFbyIND.$popTest.e --mem=6GB --time=10:00:00 \
			--wrap "Rscript /pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Tools/sourcefindV2/sourcefindv2_PierreCall.R $folder/fineStructure/Outputs/stage6/$popTest/stage6.SFbyIND.$popTest.params"
	else
		echo /pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage6/$popTest/stage6.$popTest.SFbyIND.saveout already done
	fi
done

