#!/bin/bash

step=$1
qos=$2

Kmax=15
numrep=10
root=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/
rootScript=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/Admixture

folder=Admixture
mkdir $root/$folder
cd $root/$folder

dataFolder=$root//StartFilteredData/
pref=Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.pruned

if [[ $step == 1 ]]
then
	/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Admixture/1_runAdmixture.sh $(pwd)/ $dataFolder/$pref $qos $Kmax $numrep
else
	if [[ $step == 2 ]]
	then
		/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Admixture/2_GetBestAdmixtureRuns.sh  $(pwd)/ $pref $Kmax $numrep
	else
		if [[ $step == 3 ]]
		then
			/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/Admixture/3_runAdmixtureBestLllik_withCV.sh $(pwd)/ $dataFolder/$pref $qos $Kmax
		else
			if [[ $step == 4 ]]
			then
				echo "K CV" > $root/$folder/BestRUNperK/CV.BESTruns
				K=2
				while [[ $K -le $Kmax ]]
				do
					cv=$(grep CV $root/$folder/BestRUNperK/$pref.K$K.out | awk '{print $4}')
					echo $K $cv >> $root/$folder/BestRUNperK/CV.BESTruns
					let ' K = K + 1 '
				done
				ls  $root/$folder/BestRUNperK/CV.BESTruns
				rm $root/$folder/BestRUNperK/CrossValidation.pdf
				Rscript $rootScript/2_MakePlotAdmixture.R $root/$folder/ $pref $Kmax
				if [[ ! -e $root/$folder/BestRUNperK/CrossValidation.pdf ]]
				then
					echo pdf not generated correctly for $folder
					exit
				fi
			else
				echo step $step not recognized
				exit
			fi
		fi
	fi
fi

