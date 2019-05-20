#!/usr/bin/perl
# 2016/12/20 GC含量を調べる

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV != 6 or $ARGV[0] eq '--help'){
	die "Usage : $0 -f [fasta file] -r [resolution] -o [output file]\n";
}
my %opt;
getopts("f:o:r:", \%opt);

my $FILE_fasta = $opt{f};
my $RESOLUTION = $opt{r};
my $FILE_out = $opt{o};


my $chr = '';
my $seq = '';
my $fh_fasta = IO::File->new($FILE_fasta) or die "cannot open $FILE_fasta: $!";
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";

while($_ = $fh_fasta->getline()){
	s/\r?\n//;
	if(m/^>(\w+)$/){
		if($seq ne ""){
			&getGC($chr, $seq);
		}
		$chr = $1;
		$seq = "";
		next;
	}
	$seq .= $_;
}
&getGC($chr, $seq);
$fh_fasta->close();
$fh_out->close();

sub getGC{
	my ($chr, $seq) = @_;
	$seq = uc($seq);
	for(my $i=0; $i < length($seq); $i+= $RESOLUTION){
		my $start = $i;
		my $end = $i + $RESOLUTION - 1;
		my $len = $end - $start + 1;
		my $target_seq = substr($seq, $start, $RESOLUTION);
		my $GC = $target_seq =~ tr/GC//;
		my $AT = $target_seq =~ tr/AT//;
		my $percentage = 'NA';
		if(($GC + $AT) != 0){
			$percentage = $GC / ($GC + $AT);
		}
		$fh_out->print(join("\t", $chr, $start, $end, $GC, $AT, $percentage) . "\n");
	}
}

