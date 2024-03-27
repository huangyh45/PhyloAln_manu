#!/usr/bin/perl

use strict;
my $indir=$ARGV[0];
my $outdir=$ARGV[1];
my $taxa=$outdir;
$taxa=~s/$indir\///;
open(F,"<$indir/orthograph.conf") or die "\nError: $indir/orthograph.conf can't open!\n";
open(FF,">$outdir/$taxa.conf") or die "\nError: $outdir/$taxa.conf can't open!\n";
while(<F>) {
    if($_=~/input-file         = /) {print FF "input-file         = $outdir/$taxa.fas\n";}
    elsif($_=~/species-name       = /) {print FF "species-name       = $taxa\n";}
    elsif($_=~/output-directory   = /) {print FF "output-directory   = $outdir\n";}
    else {print FF $_;}
}
close(FF);
