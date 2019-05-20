#!/usr/bin/perl
# 2015/10/18 分割して出力したデータベースの中身をカウントする

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;

if(@ARGV != 2 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [read file]\n";
}

my %opt;
getopts("i:", \%opt);
my $FILE_target = $opt{i};

#---------------------------------------
# read file and count
#---------------------------------------
my %data;
my $fh_in = IO::File->new($FILE_target) or die "cannot open $FILE_target: $!";
while($_ = $fh_in->getline()){
	s/\r?\n//;
	$data{$_} += 1;
}
$fh_in->close();



my $FILE_counted = $FILE_target . '_counted';

#---------------------------------------
# 出力する
#---------------------------------------
my $fh_out = IO::File->new($FILE_counted, 'w') or die "cannot write $FILE_counted: $!";
foreach my $key(keys %data){
	$fh_out->print("$key\t$data{$key}\n");
}
$fh_out->close();


rename $FILE_counted, $FILE_target or die "cannot rename file: $!";


