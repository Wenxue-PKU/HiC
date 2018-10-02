#!/usr/bin/perl
# # 2017/02/16 binあたりの制限酵素部位の数を調べる

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV != 8 or $ARGV[0] eq '--help'){
	die "Usage : $0 -f [fasta file] -r [resolution] -o [output file] -t [restriction sequence e.x. GATC for MboI]\n";
}
my %opt;
getopts("f:o:r:t:", \%opt);

my $FILE_fasta = $opt{f};
my $RESOLUTION = $opt{r};
my $FILE_out = $opt{o};
my $recognitionSeq = $opt{t};
$recognitionSeq = uc($recognitionSeq);

my $chr = '';
my $seq = '';
my $fh_fasta = IO::File->new($FILE_fasta) or die "cannot open $FILE_fasta: $!";
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";

while($_ = $fh_fasta->getline()){
	s/\r?\n//;
	if(m/^>(\w+)$/){
		if($seq ne ""){
			&getFragNum($chr, $seq);
		}
		$chr = $1;
		$seq = "";
		next;
	}
	$seq .= $_;
}
&getFragNum($chr, $seq);
$fh_fasta->close();
$fh_out->close();

sub getFragNum{
	my ($chr, $seq) = @_;
	$seq = uc($seq);
	for(my $i=0; $i < length($seq); $i+= $RESOLUTION){
		my $start = $i;
		my $end = $i + $RESOLUTION - 1;
		my $len = $end - $start + 1;
		my $target_seq = substr($seq, $start, $RESOLUTION);

		my $NUM_restriction = 0;
		my $pos = 0;
		my $previous = 0;
		while($pos != -1){
			$pos = index $target_seq, $recognitionSeq, $previous;
			if($pos != -1){
				$NUM_restriction++;
			}
			$previous = $pos + 1;
		}

		$fh_out->print(join("\t", $chr, $start, $end, $NUM_restriction) . "\n");
	}
}


