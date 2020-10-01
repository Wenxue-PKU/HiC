#!/usr/bin/perl
# 2018-03-09 convert border information to bed

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV < 4 or $ARGV[0] eq '--help'){
    die "Usage : $0 -i [input border definition file] -o [output bed file] -c [column of border(default:6 border strength, DI:5)] -t [1 to get only true tad. 0 (default) for all region between border]\n";
}

my %opt;
getopts("i:o:c:t:", \%opt);
my $FILE_in = $opt{i};
my $FILE_out = $opt{o};
my $COLUMN = $opt{c};
unless(defined $COLUMN){
    $COLUMN = 6;
}
$COLUMN = $COLUMN - 1;
my $FLAG_true = $opt{t};
unless(defined $FLAG_true){
    $FLAG_true = 0;
}



my $fh_in = IO::File->new($FILE_in) or die "cannot open $FILE_in: $!";
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";

my $inTAD = 0;
my $C = '';
my $S = 0;
my $E = 0;
while($_ = $fh_in->getline()){
    s/\r?\n//;
    my @data = split /\t/;
    my ($chr, $start, $end) = @data[0,1,2];
    my $border = $data[$COLUMN];
	my $OKTAD = 1;
	if($FLAG_true == 1){
		$OKTAD = $data[7];
	}
	if($border == 0 and $OKTAD == 1){
		if($inTAD == 0){
			$C = $chr;
			$S = $start;
			$E = $end;
			$inTAD = 1;
		}else{
			$E = $end;
		}
	}
    if($border == 1){
        # if($S != $start and $C eq $chr){
		if($inTAD == 1){
            $fh_out->print("$C\t$S\t$E\n");
        }
		$inTAD = 0;
    }
}
if($inTAD == 1){
	$fh_out->print("$C\t$S\t$E\n");
}
$fh_in->close();
$fh_out->close();




