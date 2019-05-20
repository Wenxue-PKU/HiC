#!/usr/bin/perl
# 2019-01-03 intra EBVについてsingle enzyme resolutionのmatrixを作成

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;

use DBI;

if(@ARGV != 6 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [detabase files] -o [output prefix] -e [enzyme site file]\n";
}

my %opt;
getopts("i:o:e:", \%opt);
my $FILE_database = $opt{i};
my $FILE_out_prefix = $opt{o};
my $FILE_enzyme = $opt{e};


# 全データを入れる変数
my %data;
my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");

#---------------------------------------
# collect information
#---------------------------------------
my $sth_data;

# intra-chromosomeのデータのみを取得する
$sth_data = $dbh->prepare("select chr1, start1, end1, fragNum1, chr2, start2, end2, fragNum2, score from fragment where chr1=chr2 and chr1=='EBV';");
$sth_data->execute();
while(my $ref = $sth_data->fetchrow_arrayref()){
	my ($chr1, $start1, $end1, $frag1, $chr2, $start2, $end2, $frag2, $score) = @$ref;

	# 隣同士のfragmentはカウントしない
	if(abs($frag1 - $frag2) < 2){
		next;
	}

	my $middle1 = ($start1 + $end1) / 2;
	my $middle2 = ($start2 + $end2) / 2;
	my $distance = abs($middle1 - $middle2);

	# 20kb以内の距離だった場合には、scoreを２倍にする
	# (20kb以内については、同じ向きのデータしか無いから)
	if($distance < 20000){
		$score = $score * 2;
	}

	# count data（既に左側が小さいということは保証されている）
	my $bin1 = $chr1 . ":" . $start1 . ":" . $end1;
	my $bin2 = $chr2 . ":" . $start2 . ":" . $end2;
	$data{$bin1}{$bin2} += $score;
}
$sth_data->finish();
$dbh->disconnect();



#---------------------------------------
# output
#---------------------------------------
my $FILE_out = $FILE_out_prefix . "EBV.matrix";
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";

# register bins
my @bins;
my $fh_in = IO::File->new($FILE_enzyme) or die "cannot open $FILE_enzyme: $!";
$fh_in->getline();
while($_ = $fh_in->getline()){
	s/\r?\n//;
	my ($num, $chr, $position, $len_before, $len_after) = split /\t/;
	if($chr ne "EBV"){
		next;
	}
	my $start = $position - $len_before;
	if($start == 0){
		$start = 1;
	}
	my $end = $position;
	my $key = $chr . ":" . $start . ":" . $end;
	push @bins, $key;
	$fh_out->print("\t$key");
}
$fh_in->close();
$fh_out->print("\n");

for(my $i = 0; $i < @bins; $i++){
	my @values;
	for(my $j = 0; $j < @bins; $j++){
		if($i < $j){
			my $value = exists $data{$bins[$i]}{$bins[$j]} ? $data{$bins[$i]}{$bins[$j]} : 0;
			push @values, $value;
		}else{
			my $value = exists $data{$bins[$j]}{$bins[$i]} ? $data{$bins[$j]}{$bins[$i]} : 0;
			push @values, $value;
		}
	}
	$fh_out->print("$bins[$i]\t");
	$fh_out->print(join("\t", @values) . "\n");
}
$fh_out->close();



