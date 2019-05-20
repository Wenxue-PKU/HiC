#!/usr/bin/perl
# 2016/04/06 originalのマトリックスのフォーマットをHiTCように変換する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if((@ARGV != 4 and @ARGV != 5) or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [original matrix] -o [output matrix] -p \n";
}

my %opt;
getopts("i:o:p", \%opt);
my $FILE_in = $opt{i};
my $FILE_out = $opt{o};
my $FLAG_reverse = $opt{p};

my %chrCnv = (
	'I' => 'chr1',
	'II' => 'chr2',
	'III' => 'chr3'
);

my $fh_in;
if($FILE_in =~ /\.gz/){
	$fh_in = IO::File->new("gzip -dc $FILE_in |") or die "cannot open $FILE_in: $!";
}else{
	$fh_in = IO::File->new($FILE_in) or die "cannot open $FILE_in: $!";
}
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";
my $TITLE = $fh_in->getline();
$TITLE =~ s/\r?\n//;
my @titles = split /\t/, $TITLE;
my @names;
my $binNUM=0;
foreach my $v(@titles){
	my ($chr, $start, $end);
	
	if($v =~ m/(\w+):(\d+):(\d+)/){
		$chr = $1;
		$start= $2;
		$end= $3;
		$binNUM++;
		if(exists $chrCnv{$chr}){
			$chr =  $chrCnv{$chr};
		}
		my $n = sprintf "bin%05d|sample|%s:%d-%d", $binNUM, $chr, $start, $end;
		push @names, $n;
	}
}
$fh_out->print(join("\t", "", @names) . "\n");
my $dataNum = -1;
while($_ = $fh_in->getline()){
	s/\r?\n//;
	$dataNum++;
	my @data = split /\t/;
	shift @data;

	if(defined $FLAG_reverse){
		my @newData;
		foreach my $v (@data){
			if($v eq "NA"){
				push @newData, "NA";
			}else{
				my $v2 = -1 * $v;
				push @newData, $v2;
			}
		}
		@data = @newData;
	}

	$fh_out->print("$names[$dataNum]\t");
	$fh_out->print(join("\t", @data) . "\n");
}
$fh_in->close();
$fh_out->close();
