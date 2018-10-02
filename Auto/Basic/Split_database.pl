#!/usr/bin/perl
# 2015/10/18 filterするためにデータベースのデータのフォーマット変換した後、分割して出力する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;

use DBI;

if((@ARGV != 8 and @ARGV != 10) or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [database] -x [target organism (pombe or human or mouse or pombe_HindIII) -r [restriction site] -o [file list] -m [mapQ threshold (default:10)]\n";
}

my %chromosome_length = (
	'human' => 3095677412,
	'human_EBV' => 3157782322,
	'mouse' => 2725537669,
	'pombe' => 12571820,
);

my %opt;
getopts("i:x:o:r:m:", \%opt);
my $FILE_database = $opt{i};
my $TARGET_ORGANISM = $opt{x};
my $TARGET_RESTRICTION = $opt{r};
my $FILE_list = $opt{o};
my ($name, $dir, $ext) = &fileparse($FILE_database, '\..*');
my $TMP_PREFIX = $dir . 'tmp_database_' . $name;
my $MAPQ_threshold = $opt{m};
unless(defined $MAPQ_threshold){
	$MAPQ_threshold = 10;
}


# ファイルを作成する解像度
my $RESOLUTION = $chromosome_length{$TARGET_ORGANISM} / 100;

#---------------------------------------
# read Hind III file
#---------------------------------------
my $FILE_HINDIII_def;
if($TARGET_RESTRICTION eq 'human'){
	$FILE_HINDIII_def = '/wistar/noma/Data/Human_seq/hg19/MboI_sites.txt';
}elsif($TARGET_RESTRICTION eq 'human_EBV'){
	$FILE_HINDIII_def = '/wistar/noma/Data/Human_seq/hg19_EBV/MboI_sites.txt';
}elsif($TARGET_RESTRICTION eq 'human_HindIII'){
	$FILE_HINDIII_def = '/wistar/noma/Data/Human_seq/hg19/HindIII_sites.txt';
}elsif($TARGET_RESTRICTION eq 'mouse'){
	$FILE_HINDIII_def = '/wistar/noma/Data/Mouse_seq/mm10/MboI_sites.txt';
}elsif($TARGET_RESTRICTION eq 'pombe'){
	$FILE_HINDIII_def = '/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18/MboI_sites.txt';
}elsif($TARGET_RESTRICTION eq 'pombe_HindIII'){
	$FILE_HINDIII_def = '/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18/HindIII_sites.txt';
}else{
	die "please specify correct organism name\n";
}
my %Hinds;
my %Chromosomes;
{
	my $fh_in = IO::File->new($FILE_HINDIII_def) or die "cannot open $FILE_HINDIII_def: $!";
	while($_ = $fh_in->getline()){
		if(m/^#/){
			next;
		}
		s/\r?\n//;
		my ($number, $chr, $pos, $before, $after) = split /\t/;

		# 番号0の断片を作っておく
		if($number == 1){
			my $id0 = $chr . ':0';
			$Hinds{$id0} = "1\t$pos";
		}

		my $id = $chr . ':' . $number;
		my $end = $pos + $after;
		$Hinds{$id} = "$pos\t$end";
		$Chromosomes{$chr} = 1;
	}
	$fh_in->close();
}


# ファイルハンドル
my %FH;
my %NAMES;

my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");
my $sth = $dbh->prepare("select chr1, restNum1, restSide1, position1, direction1, restLoc1, chr2, restNum2, restSide2, position2, direction2, restLoc2 from map
	where uniq1 = 'U' and uniq2 = 'U' and mapQ1 > $MAPQ_threshold and mapQ2 > $MAPQ_threshold and restLoc1 != 'NA' and restLoc2 != 'NA'");
$sth->execute();
while(my $ref = $sth->fetchrow_arrayref()){
	my ($chr1, $restNum1, $restSide1, $position1, $direction1, $restLoc1, $chr2, $restNum2, $restSide2, $position2, $direction2, $restLoc2) = @$ref;

	# もし以下のケースだった場合には、最も近い制限酵素部位の位置が、正しいソースでないので、位置を調整する
	# 向きが-で、align位置 < 制限酵素部位。restrictionnumをマイナス1する。
	# 向きが+で、align位置 > 制限酵素部位。restrictionnumをプラス1する。
	if($direction1 eq '+' and $position1 > $restLoc1){
		$restNum1++;
	}
	if($direction1 eq '-' and $position1 < $restLoc1){
		$restNum1--;
	}
	if($direction2 eq '+' and $position2 > $restLoc2){
		$restNum2++;
	}
	if($direction2 eq '-' and $position2 < $restLoc2){
		$restNum2--;
	}


	if($restSide1 eq 'L'){
		$restNum1--;
	}
	if($restSide2 eq 'L'){
		$restNum2--;
	}

	my $id1 = $chr1 . ':' . $restNum1;
	my $id2 = $chr2 . ':' . $restNum2;

	my ($start1, $end1) = split /\t/, $Hinds{$id1};
	my ($start2, $end2) = split /\t/, $Hinds{$id2};


	my $middle1 = ($start1 + $end1) / 2;
	my $middle2 = ($start2 + $end2) / 2;

	# 20kb以下の距離で、向きが異なる組み合わせについては除去する
	if($chr1 eq $chr2 and abs($middle1 - $middle2) < 20000 and $direction1 ne $direction2){
		next;
	}

	my $bin = int($middle1 / $RESOLUTION);

	# もしファイルハンドルが未生成だったら作成（合計で約１００個できる可能性あり）
	unless(exists $FH{$chr1}{$bin}){
		my $file_out = $TMP_PREFIX . '_' . $chr1 . '_' . $bin;
		my $f = IO::File->new($file_out, 'w') or die "cannot write $file_out: $!";
		$FH{$chr1}{$bin} = $f;
		$NAMES{$chr1}{$bin} = $file_out;
	}
	$FH{$chr1}{$bin}->print("$chr1\t$start1\t$end1\t$restNum1\t$chr2\t$start2\t$end2\t$restNum2\n");
}
$sth->finish();
$dbh->disconnect();


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

