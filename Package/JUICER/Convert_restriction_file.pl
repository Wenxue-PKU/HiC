#!/usr/bin/perl
# 2020-02-25 Convert restriction site file for JUICER box

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV < 4 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [restriction file (ex. MboI_sites.txt)] -o [output file (MboI_sites.juicer.txt)]\n";
}
my %opt;
getopts("i:o:", \%opt);
my $FILE_in = $opt{i};
my $FILE_out = $opt{o};

my $fh_in = IO::File->new($FILE_in) or die "cannot open $FILE_in: $!";
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";
my $chr_pre = 'dummy';
while($_ = $fh_in->getline()){
	s/\r?\n//;
	if(m/^#/){
		next;
	}
	my ($num, $chr, $position, $before, $after) = split /\t/, $_;
	if($position == 0){
		next;
	}
	if($chr_pre ne $chr){
		if($chr_pre eq 'dummy'){
			$fh_out->print("$chr");
		}else{
			$fh_out->print("\n$chr");
		}
	}else{
		$fh_out->print(" $position");
	}
	$chr_pre = $chr;
}
$fh_in->close();
$fh_out->close();
