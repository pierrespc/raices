#!/bin/perl

use strict;
use warnings;


my $fam="/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/StartFilteredData/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP.Filtered.2ndDegree.SubSampling_REDUCED.fam";



chdir("/pasteur/zeus/projets/p02/Hotpaleo/pierre/Projects/RAICES/fineStructure/Outputs/stage3/")|| die("cannot change directory");



open(FAM,$fam)|| die("can't open $fam\n");

my %hash;
foreach my $line(<FAM>){
	chomp $line;
	my @split=split(/\s+/,$line);
	$hash{$split[1]}=$split[0];
}
close(FAM);

#foreach my $key (keys %hash){
#	print $key." ".$hash{$key}."\n";
#}

for my $treeXML ("stage3.tree.Other.xml","stage3.tree.xml"){
	open(XML,$treeXML) || die("can't open $treeXML\n");
	open(OUT,">$treeXML.Parsed.xml") || die("can't write $treeXML.Parsed.xml\n");

	foreach my $line(<XML>){
        	chomp $line;
		if($line =~ /<Tree>/){
			#print $line."\n";
			foreach my $key (keys %hash){
				$line=~s/\($key:/\($key-$hash{$key}:/g;
				$line=~s/,$key:/,$key-$hash{$key}:/g;
			}
			#print $line."\n";
		}
		if($line =~ /<Pop>/ ){
			foreach my $key (keys %hash){
				$line=~s/\($key\)/\($key-$hash{$key}\)/g;
				$line=~s/,$key,/,$key-$hash{$key},/g;
				$line=~s/\($key,/\($key-$hash{$key},/g;
				$line=~s/,$key\)/,$key-$hash{$key}\)/g;
				
                        }
		}
		print OUT $line."\n";
	}
	close(XML);
	close(OUT);
}
