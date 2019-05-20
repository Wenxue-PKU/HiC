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
	die "Usage : $0 -i [database] -l [chromosome length] -o [file list] -m [mapQ threshold (default:10)] -e [enzyme definition file]\n";
}

my %opt;
getopts("i:l:o:m:e:", \%opt);
my $FILE_database = $opt{i};
my $CHROM_LENGTH = $opt{l};
my $FILE_list = $opt{o};
my ($name, $dir, $ext) = &fileparse($FILE_database, '\..*');
my $TMP_PREFIX = $dir . 'tmp_database_' . $name;
my $MAPQ_threshold = $opt{m};
unless(defined $MAPQ_threshold){
	$MAPQ_threshold = 10;
}

my $FILE_ENZYME_def = $opt{e};


# ファイルを作成する解像度
my $RESOLUTION = $CHROM_LENGTH / 100;

#---------------------------------------
# read restriction information
#---------------------------------------
my %Enzymes;
my %Chromosomes;
{
	my $fh_in = IO::File->new($FILE_ENZYME_def) or die "cannot open $FILE_ENZYME_def: $!";
	while($_ = $fh_in->getline()){
		if(m/^#/){
			next;
		}
		s/\r?\n//;
		my ($number, $chr, $pos, $before, $after) = split /\t/;

		# 番号0の断片を作っておく
		if($number == 1){
			my $id0 = $chr . ':0';
			$Enzymes{$id0} = "1\t$pos";
		}

		my $id = $chr . ':' . $number;
		my $end = $pos + $after;
		$Enzymes{$id} = "$pos\t$end";
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

	my ($start1, $end1) = split /\t/, $Enzymes{$id1};
	my ($start2, $end2) = split /\t/, $Enzymes{$id2};


	my $middle1 = ($start1 + $end1) / 2;
	my $middle2 = ($start2 + $end2) / 2;

	# 10kb以下の距離で、向きが異なる組み合わせについては除去する
	if($chr1 eq $chr2 and abs($middle1 - $middle2) < 10000 and $direction1 ne $direction2){
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

