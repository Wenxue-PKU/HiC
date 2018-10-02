#!/usr/bin/perl
# 2015/06/16 ソートのために１００のファイルに分割する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;



if(@ARGV != 6 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [map file] -x [target organism (pombe or human) -o [file list]\n";
}

my %chromosome_length = (
	'human' => 3095677412,
	'human_EBV' => 3157782322,
	'mouse' => 2725537669,
	'pombe' => 12571820,
);

my %opt;
getopts("i:x:o:", \%opt);
my $FILE_map = $opt{i};
my $TARGET_ORGANISM = $opt{x};
my $FILE_list = $opt{o};
my ($name, $dir, $ext) = &fileparse($FILE_map, '\..*');
my $TMP_PREFIX = $dir . 'tmp_register_' . $name;


unless(exists $chromosome_length{$TARGET_ORGANISM}){
	die "target organism is unknown\n";
}


# ファイルを作成する解像度
my $RESOLUTION = $chromosome_length{$TARGET_ORGANISM} / 100;

# ファイルハンドル
my %FH;
my %NAMES;

#---------------------------------------
# read file and split and save
#---------------------------------------
my $TOTAL_READ_NUMBER = 0;
my $fh_map = IO::File->new($FILE_map) or die "cannot open $FILE_map: $!";
while($_ = $fh_map->getline()){
	s/\r?\n//;
	my ($id, $chr1, $loc1, $direction1, $mapQ1, $resID1, $resLoc1, $uniq1, $chr2, $loc2, $direction2, $mapQ2, $resID2, $resLoc2, $uniq2) = split /\t/;
	if($resID1 eq 'NA' or $resID2 eq 'NA'){
		next;
	}


	# 右と左を比べて小さい方を左にするように統一する
	if($chr1 eq $chr2){
		if($loc1 > $loc2){
			($loc2, $loc1) = ($loc1, $loc2);
			($direction2, $direction1) = ($direction1, $direction2);
			($mapQ2, $mapQ1) = ($mapQ1, $mapQ2);
			($resID2, $resID1) = ($resID1, $resID2);
			($resLoc2, $resLoc1) = ($resLoc1, $resLoc2);
			($uniq2, $uniq1) = ($uniq1, $uniq2);
		}
	}else{
		if($TARGET_ORGANISM eq 'pombe'){
			my $comparison = $chr1 cmp $chr2;
			if($comparison == 1){
				($chr2, $chr1) = ($chr1, $chr2);
				($loc2, $loc1) = ($loc1, $loc2);
				($direction2, $direction1) = ($direction1, $direction2);
				($mapQ2, $mapQ1) = ($mapQ1, $mapQ2);
				($resID2, $resID1) = ($resID1, $resID2);
				($resLoc2, $resLoc1) = ($resLoc1, $resLoc2);
				($uniq2, $uniq1) = ($uniq1, $uniq2);
			}
		}else{
			if(&ComparisonChr($chr1, $chr2) == 1){
				($chr2, $chr1) = ($chr1, $chr2);
				($loc2, $loc1) = ($loc1, $loc2);
				($direction2, $direction1) = ($direction1, $direction2);
				($mapQ2, $mapQ1) = ($mapQ1, $mapQ2);
				($resID2, $resID1) = ($resID1, $resID2);
				($resLoc2, $resLoc1) = ($resLoc1, $resLoc2);
				($uniq2, $uniq1) = ($uniq1, $uniq2);
			}
		}
	}


	my $bin = int($loc1 / $RESOLUTION);

	# もしファイルハンドルが未生成だったら作成（合計で約１００個できる可能性あり）
	unless(exists $FH{$chr1}{$bin}){
		my $file_out = $TMP_PREFIX . '_' . $chr1 . '_' . $bin;
		my $f = IO::File->new($file_out, 'w') or die "cannot write $file_out: $!";
		$FH{$chr1}{$bin} = $f;
		$NAMES{$chr1}{$bin} = $file_out;
	}

	my $LINE = join("\t", $id, $chr1, $loc1, $direction1, $mapQ1, $resID1, $resLoc1, $uniq1, $chr2, $loc2, $direction2, $mapQ2, $resID2, $resLoc2, $uniq2);
	$FH{$chr1}{$bin}->print("$LINE\n");
}
$fh_map->close();


sub ComparisonChr{
	my ($chr1, $chr2) = @_;
	my ($num1, $num2) = ($chr1, $chr2);
	if($chr1 =~ m/chr(\w+)/){
		$num1 = $1;
	}
	if($chr2 =~ m/chr(\w+)/){
		$num2 = $1;
	}
	if($num1 eq 'X'){
		$num1 = 23;
	}
	if($num2 eq 'X'){
		$num2 = 23;
	}
	if($num1 eq 'Y'){
		$num1 = 24;
	}
	if($num2 eq 'Y'){
		$num2 = 24;
	}
	if($num1 eq 'M'){
		$num1 = 25;
	}
	if($num2 eq 'M'){
		$num2 = 25;
	}
	if($num1 eq 'EBV'){
		$num1 = 26;
	}
	if($num2 eq 'EBV'){
		$num2 = 26;
	}

	if($num1 < $num2){
		return -1;
	}else{
		return 1;
	}
}




#---------------------------------------
# ファイルのリストを作成する &&
#---------------------------------------
my $fh_list = IO::File->new($FILE_list, 'w') or die "cannot write $FILE_list: $!";
foreach my $chr(sort {$a cmp $b} keys %FH){
	foreach my $bin(sort {$a <=> $b} keys %{$FH{$chr}}){
		$FH{$chr}{$bin}->close();
		$fh_list->print("$NAMES{$chr}{$bin}\n");
	}
}
$fh_list->close();
