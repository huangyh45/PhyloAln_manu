#!/usr/bin/env perl

use Getopt::Std;
use Bio::Tools::CodonTable;
use strict;
our $gffreaddir='/public/lihaosen/anaconda3/bin';
our $orthofinderdir='/public/lihaosen/anaconda3/envs/orthofinder/bin';
our $busco='/public/lihaosen/anaconda3/envs/busco/bin/busco';
our $buscodir='/public/lihaosen/anaconda3/envs/busco/bin';
our $buscolineage='insecta';
our $protein_annotate='';
our $kinfindir='/public/lihs/annotation/kinfin:/public/lihs/annotation/kinfin/scripts';
our $pythondir='/public/lihaosen/anaconda3/envs/python2.7/bin'; #keep empty if it is default PATH

# the input parameters
my %opt=('o'=>'ortho','p'=>8,'x'=>'.fna','y'=>'.gff','z'=>'.faa','w'=>'.cds','v'=>'.annotation.tsv','u'=>'.faa.tsv','l'=>$buscolineage);
getopts('o:p:g:f:a:c:n:i:x:y:z:w:v:u:l:srbetkh',\%opt);
usage() if $opt{h};
my $outdir=$opt{o};
my $cpu=$opt{p};
my $genomedir=$opt{g};
my $genomesuffix=$opt{x};
my $gffdir=$opt{f};
my $gffsuffix=$opt{y};
my $aadir=$opt{a};
my $aasuffix=$opt{z};
my $cdsdir=$opt{c};
my $cdssuffix=$opt{w};
my $annodir=$opt{n};
my $annosuffix=$opt{v};
my $iprdir=$opt{i};
my $iprsuffix=$opt{u};
my $nofilter=$opt{s};
my $norename=$opt{r};
my $nobusco=$opt{b};
my $lineage=$opt{l};
my $nomanual=$opt{e};
my $noannotate=$opt{t};
my $nokinfin=$opt{k};

sub usage {

die "
perl $0
Group the proteins and generate the orthogroups by OrthoFinder and annotate the orthogroups by KinFin

Usage: 
-o   the output directory(default='ortho')
-p   the maximum CPUs to be used(default=8)
-g   the genome FASTA files
-x   the suffix of the genome FASTA files(default='.fna')
-f   the directory containing the prediction GFF3 files
-y   the suffix of the prediction GFF3 files(default='.gff')
-a   the directory containing the protein FASTA files
-z   the suffix of the protein FASTA files(default='.faa')
-c   the directory containing the CDS FASTA files
-w   the suffix of the CDS FASTA files(default='.cds')
-n   the directory containing the annotation TSV files
-v   the suffix of the annotation TSV files(default='.annotation.tsv')
-i   the directory containing the InterProScan annotation TSV files
-u   the suffix of the InterProScan annotation TSV files(default='.faa.tsv')
-s   skip the step to filter out the invalid protein sequences with '*', '.' or '?', the termination codons will be excluded automatically ignoring this parameter
-r   skip the step to rename the sequences
-b   skip the step to run BUSCO to examine the completeness of the protein sets
-l   the lineage of BUSCO analysis(default='$buscolineage')
-e   continue automatically after preparing the data and running BUSCO analysis without manual check
-t   skip the step to annotate the proteins when '-n' is not provided
-k   skip the step to annotate the orthogroups
-h   this help message

Example:
Group the proteins from the genomes and prediction GFF3 files and generate the orthogroups with BUSCO analysis and annotation using 35 CPUs:
$0 -g genomedir -f gffdir -o outputdir -p 35
Group the proteins(some species from the genomes and prediction GFF3 files and others from AA and CDS sequence files, all with annotations in the format of '.txt') and generate the orthogroups without renaming the sequences:
$0 -g genomedir -f gffdir -a aadir -c cdsdir -o outputdir -p 35 -n annotationdir -v .txt -r

";
}

die "\nError: no input or output directory was set!\n" unless $outdir&&($genomedir&&$gffdir||$aadir);
unless($outdir=~/^\//) {$outdir=$ENV{'PWD'}.'/'.$outdir;}
unless(-d $outdir) {mkdir($outdir) or die "\nError: fail to create the output directory '$outdir'!\n";}
else {print "\nWarning: the output directory '$outdir' has existed, please ensure if you want to continue the previous run on the basis of this directory!\n";}
our $logdir="$outdir/log";
our $okdir="$outdir/ok";
our $tempdir="$outdir/temp";
unless(-d $okdir) {mkdir($okdir) or die "\nError: fail to create the directory '$okdir'!\n";}
unless(-d $logdir) {mkdir($logdir) or die "\nError: fail to create the directory '$logdir'!\n";}
if(`ls $tempdir`) {runcmd("rm -r $tempdir/*","prepare");}
elsif(!(-d $tempdir)) {mkdir($tempdir) or die "\nError: fail to create the directory '$tempdir'!\n";}
our $path=$ENV{'PATH'};
my (%info,%taxon);

if($genomedir&&$gffdir) {
	printf "Extracting the protein/CDS sequences and selecting the longest isoforms...\n";
	unless(-d "$outdir/gene2longest_isoform") {mkdir("$outdir/gene2longest_isoform") or die "\nError: fail to create the directory '$outdir/gene2longest_isoform'!\n";}
	opendir(DIR,$genomedir) or die "\nError: fail to open the genome directory '$genomedir'!\n";
	my @files=readdir(DIR);
	closedir(DIR);
	my $suffix=quotemeta $genomesuffix;
	foreach my $file (@files) {
		unless($file=~s/$suffix$//) {next;}
		$taxon{$file}='gff';
		if(-e "$okdir/prepare.ok") {next;}
		my $genomefile="$genomedir/$file".$genomesuffix;
		my $gfffile="$gffdir/$file".$gffsuffix;
		my (%map,$seqid,%tempinfo,%t1);
		open(F,"<$gfffile") or die "\nError: $gfffile can't open!\n";
		while(<F>) {
                if(/^#/) {next;}
                chomp;
                my @arr=split("\t",$_);
                if(@arr!=9) {next;}
                if($arr[2]=~/RNA$/i) {
                        my ($gene,$trans);
                        if($arr[8]=~/ID=([^;]+)/i) {$trans=$1;}
                        if($arr[8]=~/parent=([^;]+)/i) {$gene=$1;}
                        $map{$trans}=$gene;
                }
        }
        close(F);
		runcmd("gffread -g $genomefile -y $tempdir/$file.pep.fasta $gfffile","prepare",$gffreaddir);
		runcmd("gffread -g $genomefile -x $tempdir/$file.cds.fasta $gfffile","prepare",$gffreaddir);
		open(P,"<$tempdir/$file.pep.fasta") or die "\nError: fail to open '$tempdir/$file.pep.fasta'!\n";
		while(<P>) {
			chomp;
			if($_=~/^\>(\S+)/) {
				$seqid=$1;
				unless($map{$seqid}) {$map{$seqid}=$seqid;}
			}
			else {
				$tempinfo{$seqid}{'seq'}.="$_\n";
				$tempinfo{$seqid}{'length'}+=length($_);
			}
		}
		close(P);
		$info{$file}{'invalid'}=0;
		foreach $seqid (keys %tempinfo) {
			unless($nofilter) {
				if($tempinfo{$seqid}{'seq'}=~/[\*\.\?]/) {$info{$file}{'invalid'}+=1;next;}
			}
			my $thislen=$tempinfo{$seqid}{'length'};
			unless($t1{$map{$seqid}}) {$t1{$map{$seqid}}=$seqid;}
			else {
				my $genelen=$tempinfo{$t1{$map{$seqid}}}{'length'};
				if($thislen>$genelen) {$t1{$map{$seqid}}=$seqid;}
				elsif($thislen==$genelen) {
					my ($num1,$num2);
					if($seqid=~/(\d+)$/) {$num1=$1;}
					if($t1{$map{$seqid}}=~/(\d+)$/) {$num2=$1;}
					if($num1<$num2) {$t1{$map{$seqid}}=$seqid;}
				}
			}
		}
		open(L,">$outdir/gene2longest_isoform/$file.tsv") or die "\nError: fail to open '$outdir/gene2longest_isoform/$file.tsv'!\n";
		foreach my $geneid (keys %t1) {
			my $seqid=$t1{$geneid};
			print L "$geneid\t$seqid\n";
			$info{$file}{$seqid}{'aa'}=$tempinfo{$seqid}{'seq'};
		}
		close(L);
		open(C,"<$tempdir/$file.cds.fasta") or die "\nError: fail to open '$tempdir/$file.cds.fasta'!\n";
		while(<C>) {
			chomp;
			if($_=~/^\>(\S+)/) {
				$seqid=$1;
			}
			elsif($info{$file}{$seqid}) {
				$info{$file}{$seqid}{'cds'}.="$_\n";
			}
		}
		close(C);
	}
	unless(-e "$okdir/prepare.ok") {runcmd("rm $tempdir/*.cds.fasta $tempdir/*.pep.fasta","prepare");}
}

if($aadir) {
	printf "Reading the protein sequences...\n";
	opendir(DIR,$aadir) or die "\nError: fail to open the protein sequence directory '$aadir'!\n";
	my @files=readdir(DIR);
	closedir(DIR);
	my $suffix=quotemeta $aasuffix;
	foreach my $file (@files) {
		unless($file=~s/$suffix$//) {next;}
		if($taxon{$file}) {printf "\nWarning: '$file' is found both in the genome and protein sequence directory, use the protein sequences in '$aadir'!\n";} 
		$taxon{$file}='aa';
		if(-e "$okdir/prepare.ok") {next;}
		my $aafile="$aadir/$file".$aasuffix;
		my $seqid;
		open(A,"<$aafile") or die "\nError: fail to open '$aafile'!\n";
		while(<A>) {
			chomp;
			if($_=~/^\>(\S+)/) {
				$seqid=$1;
			}
			else {
				$info{$file}{$seqid}{'aa'}.="$_\n";
			}
		}
		close(A);
	}
}

unless(-e "$okdir/prepare.ok") {
if($cdsdir) {
	printf "Reading the CDS sequences...\n";
	opendir(DIR,$cdsdir) or die "\nError: fail to open the CDS sequence directory '$cdsdir'!\n";
	my @files=readdir(DIR);
	closedir(DIR);
	my $suffix=quotemeta $cdssuffix;
	foreach my $file (@files) {
		unless($file=~s/$suffix$//) {next;}
		unless($taxon{$file}) {
			printf "\nWarning: '$file' is found in the CDS sequence directory without protein sequences, skip it!\n";
			next;
		}
		elsif($taxon{$file} eq 'gff') {printf "\nWarning: '$file' is found both in the genome and CDS sequence directory, use the CDS sequences in '$cdsdir'!\n";}
		my $cdsfile="$cdsdir/$file".$cdssuffix;
		my $seqid;
		open(CD,"<$cdsfile") or die "\nError: fail to open '$cdsfile'!\n";
		while(<CD>) {
			chomp;
			if($_=~/^\>(\S+)/) {
				$seqid=$1;
				if($info{$file}{$seqid}{'cds'}) {$info{$file}{$seqid}{'cds'}='';}
			}
			elsif($info{$file}{$seqid}) {
				$info{$file}{$seqid}{'cds'}.="$_\n";
			}
		}
		close(CD);
	}
}

printf "Checking and managing the sequences...\n";
my $codontable=Bio::Tools::CodonTable->new(-id=>1);
foreach my $taxa (keys %taxon) {
	$info{$taxa}{'valid'}=0;
	unless($info{$taxa}{'invalid'}) {$info{$taxa}{'invalid'}=0;}
	foreach my $seqid (keys %{$info{$taxa}}) {
		if($seqid eq 'invalid'||$seqid eq 'valid') {next;}
		unless($info{$taxa}{$seqid}{'aa'}) {die "\nError: fail to find the protein sequence of '$seqid' in '$taxa'!\n";}
                unless($info{$taxa}{$seqid}{'cds'}) {die "\nError: fail to find the CDS sequence of '$seqid' in '$taxa'!\n";}
		if($taxon{$taxa} eq 'aa') {
			if($info{$taxa}{$seqid}{'aa'}=~s/\*\n$/\n/||$info{$taxa}{$seqid}{'aa'}=~s/\.\n$/\n/) {}
			unless($nofilter) {
				if($info{$taxa}{$seqid}{'aa'}=~/[\*\.\?]/) {
					$info{$taxa}{'invalid'}+=1;
					delete $info{$taxa}{$seqid};
					next;
				}
			}
		}
		my $aaseq=$info{$taxa}{$seqid}{'aa'};
		$aaseq=~s/\n//g;
		my $cdsseq=$info{$taxa}{$seqid}{'cds'};
		$cdsseq=~s/\n//g;
		if(length($cdsseq)%3!=0) {
			my $cdslen=length($cdsseq);
			if($cdslen-$cdslen%3==3*length($aaseq)) {
				my $tail=$cdslen%3;
				printf "\nWarning: the length of the CDS sequence of '$seqid' in '$taxa'(=$cdslen) cannot be devided by 3, remove the incomplete codon($tail bases) at the end of the CDS sequence!\n";
				$info{$taxa}{$seqid}{'cds'}=~s/\n?\S\n$/\n/i;
                                $cdsseq=~s/.$//;
				if($tail>1) {
					$info{$taxa}{$seqid}{'cds'}=~s/\n?\S\n$/\n/i;
					$cdsseq=~s/.$//;
				}
			}
			else {
				die "\nError: the length of the CDS sequence of '$seqid' in '$taxa'(=$cdslen) cannot be devided by 3!\n".$info{$taxa}{$seqid}{'cds'}.$info{$taxa}{$seqid}{'aa'};
			}
		}
		if($info{$taxa}{$seqid}{'cds'}=~s/\n?T\n?A\n?[AG]\n$/\n/i||$info{$taxa}{$seqid}{'cds'}=~s/\n?T\n?G\n?A\n$/\n/i) {$cdsseq=~s/T..$//i;}
		for(my $i=0;$i+2<length($cdsseq);$i=$i+3) {
			my $codon=substr($cdsseq,$i,3);
			my $base=substr($aaseq,$i/3,1);
			if(uc($codontable->translate($codon)) ne uc($base)) {
				if(uc($base) eq 'X') {
                                        my $tbase=uc($codontable->translate($codon));
                                        if($tbase eq 'B'||$tbase eq 'Z'||$tbase eq 'J') {
                                                my ($aastart,$k);
                                                my $taaseq=$info{$taxa}{$seqid}{'aa'};
                                                for(my $j=0;$j<length($taaseq);$j++) {
                                                        $base=substr($taaseq,$j,1);
                                                        if($base ne "\n") {
                                                                if($k==$i/3) {
                                                                        if($base eq 'X') {printf "\nWarning: $k-$base-$codon of the protein sequence of '$seqid' in '$taxa' will be changed into $tbase!\n";}
									else {die "\nError: fail to change $k-$base-$codon of the protein sequence of '$seqid' in '$taxa' into $tbase!\n";}
                                                                        last;
                                                                }
                                                                $k++;
                                                        }
                                                        $aastart.=$base;
                                                }
                                                $aastart=quotemeta $aastart;
                                                $info{$taxa}{$seqid}{'aa'}=~s/^($aastart)X/$1$tbase/;
                                                next;
                                        }
                                }
				my $j=$i/3;
				die "\nError: the CDS sequence of '$seqid' in '$taxa' does not match its protein sequence: $j-$base-$codon!\n";
			}
		}
		if(length($cdsseq)/3!=length($aaseq)) {
			my $cdslen=length($cdsseq);
			my $aalen=length($aaseq);
			die "\nError: the length of the CDS sequence of '$seqid' in '$taxa'(=$cdslen) does not match its protein sequence(=$aalen)!\n";
		}
		$info{$taxa}{'valid'}+=1;
	}
	printf "Taxa: $taxa, valid proteins: $info{$taxa}{'valid'}, invalid proteins: $info{$taxa}{'invalid'}\n";
	delete $info{$taxa}{'valid'};
	delete $info{$taxa}{'invalid'};
}

if($annodir) {
	printf "Reading the annotations...\n";
	unless(-d "$outdir/annotation") {mkdir("$outdir/annotation") or die "\nError: fail to create the directory '$outdir/annotation'!\n";}
	foreach my $taxa (keys %taxon) {
		my $annofile="$annodir/$taxa".$annosuffix;
		my $seqid;
		open(ANN,"<$annofile") or die "\nError: fail to open the annotation file '$annofile'!\n";
		while(my $line=<ANN>) {
			unless($info{$taxa}{'annoheader'}) {$info{$taxa}{'annoheader'}=$line;}
			elsif($line=~s/^([^\t]+)\t//) {
				$seqid=$1;
				if($info{$taxa}{$seqid}) {$info{$taxa}{$seqid}{'anno'}=$line;}
			}
		}
		close(ANN);
	}
}

if($iprdir) {
	printf "Reading the InterProScan annotations...\n";
	unless(-d "$outdir/annotation") {mkdir("$outdir/annotation") or die "\nError: fail to create the directory '$outdir/annotation'!\n";}
	foreach my $taxa (keys %taxon) {
		my $iprfile="$iprdir/$taxa".$iprsuffix;
		my $seqid;
		open(IP,"<$iprfile") or die "\nError: fail to open the InterProScan annotation file '$iprfile'!\n";
		while(my $line=<IP>) {
			if($line=~s/^([^\t]+)\t//) {
				$seqid=$1;
				if($info{$taxa}{$seqid}) {push(@{$info{$taxa}{$seqid}{'ipr'}},$line);}
			}
		}
		close(IP);
	}
}

printf "Renaming and/or preparing the data...\n";
unless(-d "$outdir/cds") {mkdir("$outdir/cds") or die "\nError: fail to create the directory '$outdir/cds'!\n";}
unless($norename) {unless(-d "$outdir/rename_map") {mkdir("$outdir/rename_map") or die "\nError: fail to create the directory '$outdir/rename_map'!\n";}}
if($iprdir) {open(IPR,">$outdir/annotation/all.interproscan.tsv") or die "\nError: fail to open '$outdir/annotation/all.interproscan.tsv'!\n";}
foreach my $taxa (keys %taxon) {
	open(AA,">$outdir/$taxa.faa") or die "\nError: fail to open '$outdir/$taxa.faa'!\n";
	open(CDS,">$outdir/cds/$taxa.cds") or die "\nError: fail to open '$outdir/cds/$taxa.cds'!\n";
	if($annodir) {
		open(ANNO,">$outdir/annotation/$taxa.annotation.tsv") or die "\nError: fail to open '$outdir/annotation/$taxa.annotation.tsv'!\n";
		print ANNO $info{$taxa}{'annoheader'};
		delete $info{$taxa}{'annoheader'};
	}
	unless($norename) {open(MAP,">$outdir/rename_map/$taxa.tsv") or die "\nError: fail to open '$outdir/rename_map/$taxa.tsv'!\n";}
	my $n=0;
	foreach my $seqid (keys %{$info{$taxa}}) {
		$n+=1;
		my $seqname=$seqid;
		unless($norename) {
			$seqname="$taxa.$n";
			print MAP "$seqid\t$seqname\n";
		}
		print AA ">$seqname\n$info{$taxa}{$seqid}{'aa'}";
		print CDS ">$seqname\n$info{$taxa}{$seqid}{'cds'}";
		if($annodir) {
			unless($info{$taxa}{$seqid}{'anno'}) {die "\nError: no annotation line of '$seqid' in '$taxa'!\n";}
			print ANNO "$seqname\t$info{$taxa}{$seqid}{'anno'}";
		}
		if($iprdir&&$info{$taxa}{$seqid}{'ipr'}) {
			my @iprs=@{$info{$taxa}{$seqid}{'ipr'}};
			for(my $i=0;$i<@iprs;$i++) {
				print IPR "$seqname\t$iprs[$i]";
			}
		}
	}
	close(AA);
	close(CDS);
	if($annodir) {close(ANNO);}
	unless($norename) {close(MAP);}
}
if($iprdir) {close(IPR);}
%info=();
open(OK,">$okdir/prepare.ok") or die "\nError: fail to open '$okdir/prepare.ok'!\n";
close(OK);
}

unless($nobusco||-e "$okdir/busco.ok") {
	printf "Examining the completeness of the protein sets by BUSCO...\n";
	unless(-d "$outdir/busco") {mkdir("$outdir/busco") or die "\nError: fail to create the directory '$outdir/busco'!\n";}
	chdir("$outdir/busco") or die "\nError: fail to enter the directory '$outdir/busco'!\n";
	foreach my $taxa (keys %taxon) {
		runcmd("$busco -i $outdir/$taxa.faa -o $taxa -m prot -l $lineage -c $cpu --offline","busco",$buscodir);
	}
	open(SUM,">$outdir/busco/summary.tsv") or die "\nError: fail to open '$outdir/busco/summary.tsv'!\n";
	print SUM "Species\tComplete BUSCOs (C)\tComplete and single-copy BUSCOs (S)\tComplete and duplicated BUSCOs (D)\tFragmented BUSCOs (F)\tMissing BUSCOs (M)\tTotal\tSummary\n";
	while(<*\/short_summary*.txt>) {
		my ($taxa,$C,$S,$D,$F,$M,$total,$summary);
		open(BUSCO,"<$_") or die "\nError: fail to open '$_'!\n";
		while(<BUSCO>) {
			if($_=~/\/([^\/]+).faa/) {$taxa=$1;}
			elsif($_=~/^\s*(C:\d.+,n:\d+)\s*$/) {$summary=$1;}
			elsif($_=~/^\s*(\d+)\s*Complete BUSCOs \(C\)\s*$/) {$C=$1;}
			elsif($_=~/^\s*(\d+)\s*Complete and single-copy BUSCOs \(S\)\s*$/) {$S=$1;}
			elsif($_=~/^\s*(\d+)\s*Complete and duplicated BUSCOs \(D\)\s*$/) {$D=$1;}
			elsif($_=~/^\s*(\d+)\s*Fragmented BUSCOs \(F\)\s*$/) {$F=$1;}
			elsif($_=~/^\s*(\d+)\s*Missing BUSCOs \(M\)\s*$/) {$M=$1;}
			elsif($_=~/^\s*(\d+)\s*Total BUSCO groups searched\s*$/) {$total=$1;}
		}
		close(BUSCO);
		print SUM "$taxa\t$C\t$S\t$D\t$F\t$M\t$total\t$summary\n";
	}
        close(SUM);
	open(OK,">$okdir/busco.ok") or die "\nError: fail to open '$okdir/busco.ok'!\n";
	close(OK);
}

my $taxnum=int(keys %taxon);
printf "\n$taxnum proteomes are ready.\n";
unless($nomanual||-e "$okdir/manual.ok") {
	printf "Please check your data and BUSCO results.\n";
	printf "If you ensure to continue the analysis, just run the command again.\n";
	open(OK,">$okdir/manual.ok") or die "\nError: fail to open '$okdir/manual.ok'!\n";
        close(OK);
	exit;
}

unless(-e "$okdir/orthofinder.ok") {
	printf "Grouping the proteins and generating the orthogroups by OrthoFinder...\n";
	my $logtext=runcmd("orthofinder -f $outdir -M msa -t $cpu -X","orthofinder",$orthofinderdir);
	my ($ofddir,$logstr);
	if($logtext=~/Species-by-species orthologues directory:\s+(\S+)\/Orthologues\//) {$ofddir=$1;}
        elsif($logtext=~/Results:\s+(\S+)\//) {$ofddir=$1;}
	else {die "\nError: fail to find the OrthoFinder output directory, please check $logdir/orthofinder.log!\n";}
	runcmd("mv $ofddir/* $outdir/OrthoFinder","orthofinder");
	rmdir($ofddir) or die "\nError: fail to delete the directory '$ofddir'!\n";
	$ofddir=quotemeta $ofddir;
	open(LIN,"<$outdir/OrthoFinder/Log.txt") or die "\nError: fail to open '$outdir/OrthoFinder/Log.txt'!\n";
	while(<LIN>) {$logstr.=$_;}
	close(LIN);
	$logstr=~s/$ofddir/$outdir\/OrthoFinder/g;
	open(LOUT,">$outdir/OrthoFinder/Log.txt") or die "\nError: fail to open '$outdir/OrthoFinder/Log.txt'!\n";
        print LOUT $logstr;
        close(LOUT);
	open(OK,">$okdir/orthofinder.ok") or die "\nError: fail to open '$okdir/orthofinder.ok'!\n";
	close(OK);
}

unless($annodir||$noannotate) {
	printf "Annotating the proteins by protein_annotate.sh...\n";
	unless(-d "$outdir/annotation") {mkdir("$outdir/annotation") or die "\nError: fail to create the directory '$outdir/annotation'!\n";}
	unless(-d "$outdir/annotation/misc") {mkdir("$outdir/annotation/misc") or die "\nError: fail to create the directory '$outdir/annotation/misc'!\n";}
	chdir("$outdir/annotation/misc") or die "\nError: fail to enter the directory '$outdir/annotation/misc'!\n";
	foreach my $taxa (keys %taxon) {
		if(-e "$okdir/annotate_$taxa.ok") {next;}
		runcmd("cp $outdir/$taxa.faa $outdir/annotation/misc","protein_annotate");
		runcmd("$protein_annotate $taxa.faa","protein_annotate");
		runcmd("rm $taxa.faa","protein_annotate");
		runcmd("mv $taxa.annotation.tsv $outdir/annotation","protein_annotate");
		open(OK,">$okdir/annotate_$taxa.ok") or die "\nError: fail to open '$okdir/annotate_$taxa.ok'!\n";
		close(OK);
	}
	runcmd("cat *.faa.tsv > $outdir/annotation/all.interproscan.tsv","protein_annotate");
}

unless($nokinfin||-e "$okdir/kinfin.ok") {
	printf "Annotating the orthogroups by KinFin...\n";
	unless(-d "$outdir/kinfin") {mkdir("$outdir/kinfin") or die "\nError: fail to create the directory '$outdir/kinfin'!\n";}
	chdir("$outdir/kinfin") or die "\nError: fail to enter the directory '$outdir/kinfin'!\n";
	open(CONF,">config.csv") or die "\nError: fail to open 'config.csv'!\n";
	print CONF "#IDX,TAXON\n";
	open(SP,"<$outdir/OrthoFinder/WorkingDirectory/SpeciesIDs.txt") or die "\nError: fail to open '$outdir/OrthoFinder/WorkingDirectory/SpeciesIDs.txt'!\n";
	while(<SP>) {
		chomp;
		s/\.faa$//;
		s/: /,/;
		print CONF "$_\n";
	}
	close(SP);
	close(CONF);
	if($pythondir) {$path="$pythondir:$ENV{'PATH'}";$ENV{'PATH'}="$pythondir:$ENV{'PATH'}";}
	runcmd("iprs2table.py -i $outdir/annotation/all.interproscan.tsv --domain_sources Pfam","kinfin",$kinfindir);
	runcmd("kinfin --cluster_file $outdir/OrthoFinder/Orthogroups/Orthogroups.txt --config_file config.csv --sequence_ids_file $outdir/OrthoFinder/WorkingDirectory/SequenceIDs.txt --species_ids_file $outdir/OrthoFinder/WorkingDirectory/SpeciesIDs.txt --fasta_dir $outdir --functional_annotation functional_annotation.txt","kinfin",$kinfindir);
	runcmd("functional_annotation_of_clusters.py all --cluster_domain_annotation kinfin_results/cluster_domain_annotation.Pfam.txt --cluster_counts_by_taxon kinfin_results/cluster_counts_by_taxon.txt --outprefix Pfam","kinfin",$kinfindir);
	runcmd("functional_annotation_of_clusters.py all --cluster_domain_annotation kinfin_results/cluster_domain_annotation.IPR.txt --cluster_counts_by_taxon kinfin_results/cluster_counts_by_taxon.txt --outprefix IPR","kinfin",$kinfindir);
	runcmd("functional_annotation_of_clusters.py all --cluster_domain_annotation kinfin_results/cluster_domain_annotation.GO.txt --cluster_counts_by_taxon kinfin_results/cluster_counts_by_taxon.txt --outprefix GO","kinfin",$kinfindir);
	runcmd("plot_cluster_size_distribution.py --infile kinfin_results/cluster_counts_by_taxon.txt","kinfin",$kinfindir);
	runcmd("generate_network.py --cluster_summary kinfin_results/TAXON/TAXON.cluster_summary.txt --config_file config.csv","kinfin",$kinfindir);
	open(OK,">$okdir/kinfin.ok") or die "\nError: fail to open '$okdir/kinfin.ok'!\n";
	close(OK);
}

sub runcmd {
open(LOG,">>$logdir/$_[1].log") or die "\nError: fail to open '$logdir/$_[1].log'!\n";
if($_[2]) {
        printf "Setting PATH: +$_[2]\n";
        print LOG "Setting PATH: +$_[2]\n";
        $ENV{'PATH'}="$_[2]:$ENV{'PATH'}";
}
printf "$_[0]\n";
print LOG "$_[0]\n";
my ($iferr,$ifwarn,$templog);
$templog=`$_[0]`;
$iferr=$?;
printf $templog;
print LOG $templog;
close(LOG);
if($templog=~/error/i || $templog=~/fault/i || $templog=~/warn/i || $templog=~/sh: line/) {$ifwarn=1;}
if($iferr) {die "\nError in '$_[0]', please check $logdir/$_[1].log!\n";}
if($ifwarn) {printf "\nWarning: there may be some errors or warnings in '$_[0]', please check $logdir/$_[1].log!\n";}
$ENV{'PATH'}=$path;
return $templog;
}
