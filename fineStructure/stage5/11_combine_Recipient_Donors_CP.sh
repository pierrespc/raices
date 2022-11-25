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
listAdmPop=$(awk -v OFS="\t" '{if($3==0)print $2}'  $folder/fineStructure/Inputs/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.MAF0.0000001.GENO1.MIND1_STEP1.ids | sort  | uniq)
mkdir logs/
for popTest in $listAdmPop
#for popTest in PuertoMadryn
do
	if [[ ! -s $popTest/stage5.$popTest.regionsquaredchunkcounts.out ]]
	then
		listIN=""
		for i in {1..22}
		do
			tmp=$listIN" "$folder/fineStructure/Outputs/stage5/$popTest/chr$i/stage5.$popTest.chr$i
			listIN=$tmp
		done
		sbatch -J comb.$popTest -o logs/stage5.$popTest.comb.o -e logs/stage5.$popTest.comb.e --mem=4G --cpus-per-task=3 --time=00:10:00 --qos fast -p common,dedicated \
			--wrap "$commandCombine -v \
			-o $folder/fineStructure/Outputs/stage5/$popTest/stage5.$popTest \
			$listIN"
	else
		echo $popTest/stage5.$popTest already done
	fi
done

cd ../stage4/FinalGroups/
listCluster=$(ls *.DonorInds)

cd ../../stage5
for clusterFile in $listCluster
do
	cluster=$(echo $clusterFile | sed  s/.DonorInds//g)
	if [[ ! -s Donors/$cluster/stage5.$cluster.regionsquaredchunkcounts.out ]]
	then
		listIN=""
		for i in {1..22}
		do
			tmp=$listIN" "$folder/fineStructure/Outputs/stage5/Donors/$cluster/chr$i/stage5.$cluster.chr$i
			listIN=$tmp
		done
		sbatch -J comb.$cluster -o logs/stage5.$cluster.comb.o -e logs/stage5.$cluster.comb.e --mem=4G --cpus-per-task=3 --time=00:10:00 --qos fast -p common,dedicated \
			 --wrap "$commandCombine -v \
			-o $folder/fineStructure/Outputs/stage5/Donors/$cluster/stage5.$cluster \
			$listIN"
	else
		echo Donors/$cluster/stage5.$cluster.regionsquaredchunkcounts.out already done
	fi
done

