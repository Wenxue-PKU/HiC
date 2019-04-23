#!/usr/bin/perl
# 2017/11/17 Mergeするために、dataForFragDbを分割する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;

use DBI;

if($ARGV[0] eq '--help'){
	die "Usage : $0 -i [dataForFragDB.txt] -l [chromosome length] -o [output sample name]\n";
}


my %opt;
getopts("i:l:o:", \%opt);
my $FILE_dataForFragDB = $opt{i};
my $CHROMOSOME_LENGTH = $opt{l};
my $OUT_NAME = $opt{o};
my $FILE_list = $OUT_NAME . "_list.txt";
my $TMP_PREFIX = "tmpFrag_" . $OUT_NAME;

# ファイルを作成する解像度
my $RESOLUTION = $CHROMOSOME_LENGTH / 100;

# ファイルハンドル
my %FH;
my %NAMES;


my $fh_in = IO::File->new($FILE_dataForFragDB) or die "cannot open $FILE_dataForFragDB: $!";
while($_ = $fh_in->getline()){
	s/\r?\n//;
	my ($chr1, $start1, $end1, $frag1, $chr2, $start2, $end2, $frag2, $count) = split /\t/;

	my $bin = int($start1 / $RESOLUTION);

	# もしファイルハンドルが未生成だったら作成（合計で約１００個できる可能性あり）
	unless(exists $FH{$chr1}{$bin}){
		my $file_out = $TMP_PREFIX . '_' . $chr1 . '_' . $bin;
		my $f = IO::File->new($file_out, 'a') or die "cannot write $file_out: $!";
		$FH{$chr1}{$bin} = $f;
		$NAMES{$chr1}{$bin} = $file_out;
	}
	$FH{$chr1}{$bin}->print("$_\n");
}
$fh_in->close();


#---------------------------------------
# ファイルのリストを作成する &&
#---------------------------------------
my $fh_list = IO::File->new($FILE_list, 'a') or die "cannot write $FILE_list: $!";
foreach my $chr(sort {$a cmp $b} keys %FH){
	foreach my $bin(sort {$a <=> $b} keys %{$FH{$chr}}){
		$FH{$chr}{$bin}->close();
		$fh_list->print("$NAMES{$chr}{$bin}\n");
	}
}
$fh_list->close();

