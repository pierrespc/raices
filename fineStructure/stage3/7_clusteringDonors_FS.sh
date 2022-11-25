#!/bin/bash



folder=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/

SinguExec="singularity exec --home $HOME:/home/$USER --bind /pasteur "

folderImg=/pasteur/zeus/projets/p02/Hotpaleo/common_data/VMs/singularity/
fsImage=$folderImg/evolbioinfo-finestructure-v4.1.1.img

folderScripts=/pasteur/zeus/projets/p02/Hotpaleo/pierre/Scripts/RAICES/fineStructure/stage3/

commandConvert="$SinguExec $fsImage impute2chromopainter.pl"
#commandConvert="perl $folderScripts/2b_impute2chromopainter.pl"
commandFs="$SinguExec $fsImage fs"

commandGreedy="$SinguExec $fsImage finestructuregreedy.sh" 
commandCombine="$SinguExec $fsImage chromocombine"

module load singularity
module load java



###Following: Gnecchi Ruscone 2019


###stage3: fineStructure


cd $folder/fineStructure/Outputs/

mkdir stage3
cd stage3
if [ ! -e stage3.mcmc.xml ]
then
	###x: 1e6 buirn in iter
	###y: 2e6 number of itereations
	###z: thinning every 1e4 iterations
	jobS3_1=$(sbatch -J S3.1 -o S3.1.o -e S3.1.e --mem=8G --cpus-per-task=1 -A metapaleo -p metapaleo \
                       --wrap "$commandFs fs \
				-x 1000000 \
				-y 2000000 \
				-z 10000 \
				../stage2/output.chunkcounts.out \
				stage3.mcmc.xml")
else
	echo stage3.mcmc.xml already generated
fi

#continue the previous run for 1M additional steps, treating the original run as burnin.
if [ ! -e stage3.mcmc.longer.xml ]
then
	jobS3_2=$(sbatch -J S3.2 -o S3.2.o -e S3.2.e --mem=8G --cpus-per-task=1 -A metapaleo -p metapaleo \
		--wrap "$commandFs fs \
		-x 0 \
		-y 1000000 \
		-z 10000 \
		../stage2/output.chunkcounts.out \
		stage3.mcmc.xml \
		stage3.mcmc.longer.xml")                       
else
	echo stage3.mcmc.longer.xml already generated
fi
## Infers a tree, using the best state seen in out.mcmc.longer.xml as the initial state.
if [ ! -e stage3.tree.xml ]
then
	jobS3_3=$(sbatch -J S3.3 -o S3.3.o -e S3.3.e --mem=8G --cpus-per-task=1 --qos=normal --time=06:00:00 \
		--wrap "$commandFs fs \
		-m T \
		-x 0 \
		-t 100000000 \
		../stage2/output.chunkcounts.out \
		stage3.mcmc.longer.xml \
		stage3.tree.xml")
else
	echo stage3.tree.xml already generated
fi


#Infers a tree, using (-T 1) the maximum concordance state over out.mcmc.xml as the initial state. This is reported with full likelihood ordering (-k 2), useful for cutting at a given number of ppulations K (but may look bad in the GUI).
if [ ! -e stage3.tree.Other.xml ]
then
	jobS3_4=$(sbatch -J S3.4 -o S3.4.o -e S3.4.e --mem=8G --cpus-per-task=1 --qos=normal --time=06:00:00 \
               --wrap "$commandFs fs \
		-m T \
		-k 2 \
		-T 1 \
		-t 100000000 \
		../stage2/output.chunkcounts.out \
		stage3.mcmc.longer.xml \
		stage3.tree.Other.xml")
fi
