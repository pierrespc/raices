#!/bin/perl


use strict;
use warnings;
my $usage="<pref>?\n";

my $pref=shof or die $usage;


foreach suf in 
folder=params[1]
chr=params[2]
step=as.numeric(params[3])
Ntot=as.numeric(params[4])


for(ind1 in seq(1,Ntot,step)){
	ind2=ind1+step
	if(ind1==1){
		out<-read.table(paste("stage2.chr",chr,".paral/",ind1,".",ind2,"
	

