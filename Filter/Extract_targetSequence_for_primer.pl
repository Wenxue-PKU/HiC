#!/usr/bin/perl
# 2015/04/02 check用のprimerをデザインするための配列を出力する。制限酵素部位の前後だけを出力。真ん中には目印にnnnnnを入力。

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);

if(@ARGV != 4 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [fasta file] -t [recognition sequence]\n";
}

# 制限酵素部位を中心にどれだけ配列をとってくるか
my $TRIM_LEFT = 100;
my $TRIM_RIGHT = 20;


my %opt;
getopts("i:t:", \%opt);
my $FastaFile = $opt{i};
my $recognitionSeq = $opt{t};
$recognitionSeq = uc($recognitionSeq);

my $fh = IO::File->new($FastaFile) or die "cannot open $FastaFile: $!";
my $id = '';
my $seq = '';
while($_ = $fh->getline()){
	s/\r?\n//;
	if(m/^>(\S+)/){
		if($seq ne ''){
			$seq = uc($seq);
			&parseSeq($id, $seq);
		}
		$id = $1;
		$seq = '';
	}else{
		$seq .= $_;
	}
}
$seq = uc($seq);
&parseSeq($id, $seq);
$fh->close();


sub parseSeq{
	my ($chr, $seq) = @_;

	my $pos = 0;
	my $previous = 0;
	my @locationList = (0);
	while($pos != -1){
		$pos = index $seq, $recognitionSeq, $previous;
		if($pos != -1){
			push @locationList, $pos;
		}
		$previous = $pos + 1;
	}
	my $totalLength = length $seq;
	push @locationList, $totalLength;

	for(my $i = 1; $i < @locationList - 1; $i++){
		my $length_before = $locationList[$i] - $locationList[$i-1];
		my $length_after = $locationList[$i+1] - $locationList[$i];


		if($length_before < 100){
			next;
		}
		if($length_after < 100){
			next;
		}

		print ">$chr $locationList[$i] left:$length_before right:$length_after\n";
		my $leftSeq = substr($seq, ($locationList[$i]-$TRIM_LEFT), $TRIM_LEFT);
		my $rightSeq = substr($seq, $locationList[$i], $TRIM_RIGHT);
		my $CombindSeq = $leftSeq . 'nnnnn' . $rightSeq;
		print "$CombindSeq\n";

	}
}