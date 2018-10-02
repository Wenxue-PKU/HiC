#!/usr/bin/perl
# 2018-03-09 convert border information to bed

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV < 4 or $ARGV[0] eq '--help'){
    die "Usage : $0 -i [input border definition file] -o [output bed file] -c [column of border(default:6 border strength, DI:5)]\n";
}

my %opt;
getopts("i:o:c:", \%opt);
my $FILE_in = $opt{i};
my $FILE_out = $opt{o};
my $COLUMN = $opt{c};
unless(defined $COLUMN){
    $COLUMN = 6;
}
$COLUMN = $COLUMN - 1;

my $fh_in = IO::File->new($FILE_in) or die "cannot open $FILE_in: $!";
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";

my $C = '';
my $S = 0;
while($_ = $fh_in->getline()){
    s/\r?\n//;
    my @data = split /\t/;
    my ($chr, $start, $end) = @data[0,1,2];
    my $border = $data[$COLUMN];
    if($border == 1){
        if($S != $start and $C eq $chr){
            $fh_out->print("$C\t$S\t$start\n");
        }
        $S = $end + 1;
        $C = $chr;
    }
}
$fh_in->close();
$fh_out->close();




