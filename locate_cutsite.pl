#!/usr/bin/perl

# The goal of this script is to go through the rainbow trout genome and identify location of cut sites for different restriction enzymes used in GBS prep

# Usage: ./locate_cutsite.pl chr.fa
# Can glob, i.e.
# ./locate_cutsite.pl *.fa

use warnings;

@chromosomes = @ARGV;

# Recognition sequence for restriction enzymes of interest
$eco = GAATTC;
$pst = CTGCAG;
$sbf = CCTGCAGG;

#Fragment length
$frag = 150;

#Outfiles for individual enzyme cut site locations
open EcoOUT, "> eco_cut_pos.csv" or die "Could not open EcoOUT";
open PstOUT, "> pst_cut_pos.csv" or die "Could not open PstOUT";
open SbfOUT, "> sbf_cut_pos.csv" or die "Could not open SbfOUT";

# Print headers
print EcoOUT "Chromosome,Position \n";
print PstOUT "Chromosome,Position \n";
print SbfOUT "Chromosome,Position \n";


open NeighborOUT, "> neighbors_eco_pst.csv" or die "Could not open NeighborOUT";
print NeighborOUT "Chromosome,Eco_cutsite,Pst_cutsite,Distance", "\n"; #print header

open SUMMARY, "> summary_all_chr.csv" or die "Could not open SUMMARY \n";
print SUMMARY "Chromosome,nEco,nPst,nSbf,nOverlap_Eco_Pst,Length \n";

$neighbors_allchr = 0;
$all_eco = 0;
$all_pst = 0;
$all_sbf = 0;

foreach $chrfile (@chromosomes){
    
    $infile = $chrfile;

    $chrfile =~ m/chr(\d+)\.fasta/;
    $chr = $1;

    print "Working on chromosome $chr, $chrfile \n";

    # Open connection to genome file
    open IN, $infile or die "Failed to open $infile \n";

    <IN>; # burn the header line

    $total_length = 0;
    $ncut_eco = 0;
    $ncut_pst = 0;
    $ncut_sbf = 0;

## Maybe want hash so can use exists() function
## Never mind, loops over arrays are not so slow. 
    @cutposEco = ();
    @cutposPst = ();
    @cutposSbf = ();

    while(<IN>){
	
	chomp;     
	$total_length = $total_length + length $_;
	
	if(m/$eco/){
	    $ncut_eco++;
	    @halves = split $eco, $_;
	    #$cutpos{$eco} = ($total_length + length $halves[0]);
	    #print $cutpos{$eco}, "\n";
	    push @cutposEco, ($total_length + length $halves[0]);
	    print EcoOUT $chr,",",($total_length + length $halves[0]), "\n";
	    
	}
	
	if(m/$pst/){
	    
	    $ncut_pst++;
	    @halves = split $pst, $_;
	    push @cutposPst, ($total_length + length $halves[0]);
	    print PstOUT $chr,",",($total_length + length $halves[0]), "\n";
	    
	}
	
	if(m/$sbf/){
	
	    $ncut_sbf++;
	    @halves = split $sbf, $_;
	    push @cutposSbf, ($total_length + length $halves[0]);
	    print SbfOUT $chr,($total_length + length $halves[0]), "\n";
	
	}

    }

    $neighbors = 0;

    foreach $i (@cutposEco){
   
	foreach $j (@cutposPst){
	
	    if($j > ($i-$frag) && $j < ($i+$frag)){
		$neighbors++;
		print NeighborOUT $chr, ",", $i, ",", "$j", ",", abs($i-$j), "\n";
	    }

	}
    }

    $neighbors_allchr = $neighbors_allchr + $neighbors;
    $all_eco = $all_eco + $ncut_eco;
    $all_pst = $all_pst + $ncut_pst;
    $all_sbf = $all_sbf + $ncut_sbf;
    
    #print SUMMARY "Chromosome,nEco,nPst,nSbf,nOverlap_Eco_Pst,Length \n";
    print SUMMARY "$chr,$ncut_eco,$ncut_pst,$ncut_sbf,$neighbors,$total_length \n";
   
    print " $ncut_eco EcoRI restriction sites \n";
    print " $ncut_pst PstI restriction sites \n";
    print " $ncut_sbf SbfI restriction sites \n";
    print "There are $neighbors sites within $frag bp of each other", "\n";
    print "Done with chromosome $chr \n";
    
}

print "Done with all chromosomes. \n $ncut_eco EcoRI cut sites \n $ncut_pst PstI cut sites \n $ncut_sbf SbfI cut sites \n $neighbors_allchr neighboring sites overall for Eco and Pst.\n";
