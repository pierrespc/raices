#!/bin/perl

use strict;
use warnings;
use Cwd;

my $usage = "runShapeIT_Prephase.pl <FolderBEDFIle> <PlinkBEDFile> <MissingInd> <MissingSNP> <MAF> <mapFile> <refHaps> <legendFile> <sampleFile> <FolderOut> <states> <WindowSize> <thread> <prune> <main> <effSize> <commandPlink> <commandShapeit> <folderScripts>


BE CAREFULL WORK FOR ONE UNIQUE CHR

FolderBEDFIle: folder where the bed file is
PlinkBEDFile: pref file for -B option MUST BE ON ONE SINGLE CHROMOSOME
MissingInd : threshold to remove indiviudals with high missing genotype rate (--mind in plink)
MissingSNP : threshold to remove snps with high missing genotype rate (--geno in plink)
MAF: threshold for maf filtering (--maf in plink)
mapFile: file for -M option
refHaps: file for -R option (1st)
legendFile: file for -R option (2nd)
sampleFile: file for -R option (3rd)
Folderout: where the outfile will be saved
States: number of states for --states option (--state option; default is 100)
WindowSize; the size of the window considered for phasing (in Mb) (--window option; default is 2)
thread: the number of thread (--thread option; default is 4)
burn: number of iterations to burn (--burn option; default is 7)
prune: number of iterations to prune (--prune option; default is 8)
main: number of main iterations (--main option; default is 7)
effSize: effective pop size (--effective-size; default is 15000)
commandPlink: a string with the command to run plink 
commandShapeit: a string with the command to ruin shapeit
folderScripts: folder with the Flipping script
\n";


my $FolderBEDFile = shift or die $usage."\n\nMISSING FolderBEDFile\n";
print "FolderBEDFile: :".$FolderBEDFile."\n";
my $PlinkBEDFile = shift or die $usage."\n\nMISSING PlinkBEDFile\n";
print "PlinkBEDFile: :".$PlinkBEDFile."\n";
my $MissingInd = shift or die $usage."\n\nMISSING MissingInd\n";
print "MissingInd: :".$MissingInd."\n";
my $MissingSNP = shift or die $usage."\n\nMISSING MissingSNP\n";
print "MissingSNP: :".$MissingSNP."\n";
my $MAF = shift or die $usage."\n\nMISSING MAF\n";
print "MAF: :".$MAF."\n";
my $mapFile = shift or die $usage."\n\nMISSING mapFile\n";
print "mapFile: :".$mapFile."\n";
my $refHaps = shift or die $usage."\n\nMISSING refHaps\n";
print "refHaps: :".$refHaps."\n";
my $legendFile = shift or die $usage."\n\nMISSING legendFile\n";
print "legendFile: :".$legendFile."\n";
my $sampleFile = shift or die $usage."\n\nMISSING sampleFile\n";;
print "sampleFile: :".$sampleFile."\n";
my $FolderOut = shift or die $usage."\n\nMISSING FolderOut\n";
print "FolderOut: :".$FolderOut."\n";
my $States = shift or die $usage."\n\nMISSING States\n";;
print "States: :".$States."\n";
my $WindowSize = shift or die $usage."\n\nMISSING WindowSizwe\n";;
print "WindowSize: :".$WindowSize."\n";
my $thread = shift or die $usage."\n\nMISSING Thread\n";;
print "thread: :".$thread."\n";
my $burn = shift or die $usage."\n\nMISSING burn\n";;
print "burn :".$burn."\n";
my $prune = shift or die $usage."\n\nMISSING prune\n";;
print "prune: :".$prune."\n";
my $main = shift or die $usage."\n\nMISSING main\n";;
print "main: :".$main."\n";
my $effSize=shift or die $usage."\n\nMISSING effSize\n";;
print "effSize: :".$effSize."\n";
my $commandPlink=shift or die $usage."\n\nMISSING commandPlink\n";;
print "commandPlink: :".$commandPlink."\n";
my $commandShapeit=shift or die $usage."\n\nMISSING commandShapeit\n";;
print "commandShapeit: :".$commandShapeit."\n";
my $folderScripts=shift or die $usage."\n\nMISSING folderScripts\n";;
print "folderScrips: :".$folderScripts."\n";

print "PROCESS CHECKING between reference ane studied set\n";
my $inPhase;
if(! -e $FolderOut."/".$PlinkBEDFile.".alignments.snp.strand"){
    system($commandShapeit."\\
     -check\\
     -B ".$FolderBEDFile."/".$PlinkBEDFile."\\
     -M ".$mapFile."\\
      --input-ref ".$refHaps."\\
     ".$legendFile."\\
      ".$sampleFile."\\
       --output-log ".$FolderOut."/".$PlinkBEDFile.".alignments");
}else{
    print "Already made\n";
}
if(! -e $FolderOut."/".$PlinkBEDFile.".alignments.snp.strand"){
	print $FolderOut."/".$PlinkBEDFile.".alignments.snp.strand Not generated\n\n
	===>>>> EVERYHING SEEMS OK!\n";
	$inPhase=$FolderBEDFile."/".$PlinkBEDFile;
	
}else{
	print "\nFLIP STRAND\n";
	system("perl ".$folderScripts."1a_Flip_accordingTolegend.pl ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand ".$FolderBEDFile."/".$PlinkBEDFile.".bim");
	#if( -e $FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.DuplicatedNames" ){
	#	print "check positions in ".$FolderOut."/".$PlinkBEDFile."\n\n It seems that you have duplicated name\n\n";
	#	print "press ENTER to keep on\n";
	#	<STDIN>;
	#}

	print "\nREMOVE SNP NOT IN LEGEND FILE\n";
	system("cat ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.NotPresent > ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude");
	system("cat ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.NotFlippable >> ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude");
	system("cat ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.DuplicatedNames >> ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude");
	system("sort ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude | uniq > ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude.tmp");
	system("mv ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude.tmp ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude");

	system("$commandPlink --bfile ".$FolderBEDFile."/".$PlinkBEDFile." --exclude ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude --maf ".$MAF." --geno ".$MissingSNP." --mind ".$MissingInd." --make-bed --out ".$FolderOut."/".$PlinkBEDFile."_alignedRef --noweb");

	#system("rm ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand");
	#system("rm ".$FolderOut."/".$PlinkBEDFile.".alignments.snp.strand.exclude");

	print "PROCESS 2nd CHECKING between reference ane studied set\n";
	system($commandShapeit." -check -B ".$FolderOut."/".$PlinkBEDFile."_alignedRef -M ".$mapFile." -R ".$refHaps." ".$legendFile." ".$sampleFile." --output-log ".$FolderOut."/".$PlinkBEDFile.".2ndAlignments");

	if( -s $FolderOut."/".$PlinkBEDFile.".2ndAlignments.snp.strand"){
	    print "your file seems to have a problem we couldn0t solve...We removed the corresponding SNPs. check the files:\n -".$FolderOut."/".$PlinkBEDFile.".2ndAlignments.snp.strand\n\n";

	    print "press ENTER to keep on\n";
	    <STDIN>;


	}else{
	    print "Alignment with reference OK!!!!\n";
	}
	$inPhase=$FolderOut."/".$PlinkBEDFile."_alignedRef";
}
print "RUN PREPHASING\n";

system($commandShapeit."\\
    -B ".$inPhase."\\
     -M ".$mapFile."\\
      --input-ref ".$refHaps."\\
     ".$legendFile."\\
      ".$sampleFile."\\
      -O ".$FolderOut."/".$PlinkBEDFile."_alignedRef_phased\\
      --window ".$WindowSize."\\
      --states ".$States."\\
      --thread ".$thread."\\
      --burn ".$burn." --prune ".$prune." --main ".$main." --effective-size ".$effSize);
#system("rm ".$FolderOut."/".$PlinkBEDFile."_alignedRef.bim");
#system("rm ".$FolderOut."/".$PlinkBEDFile."_alignedRef.bed");
#system("mv ".$FolderOut."/".$PlinkBEDFile."_alignedRef.fam ".$FolderOut."/".$PlinkBEDFile."_alignedRef_phased.fam");


