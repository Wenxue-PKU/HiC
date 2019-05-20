#!/usr/bin/perl
# Domain Caller用にmatrixを変換する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV != 4or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [original matrix] -o [output matrix]\n";
}

my %opt;
getopts("i:o:", \%opt);
my $FILE_in = $opt{i};
my $FILE_out = $opt{o};

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
$fh_in->getline();

my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";
while($_ = $fh_in->getline()){
	s/\r?\n//;
	my @data = split /\t/;
	my $LOC = shift @data;

    if($LOC =~ m/(\w+):(\d+):(\d+)/){
		my $chr = $1;
		my $start= $2;
		my $end= $3+1;
		if(exists $chrCnv{$chr}){
			$chr =  $chrCnv{$chr};
		}
        $fh_out->print(join("\t", $chr, $start, $end) . "\t");
	}

	my @newData;
	foreach my $v (@data){
		if($v eq "NA"){
			push @newData, 0;
		}else{
			push @newData, $v;
		}
	}
	@data = @newData;
	$fh_out->print(join("\t", @data) . "\n");
}
$fh_in->close();
$fh_out->close();
