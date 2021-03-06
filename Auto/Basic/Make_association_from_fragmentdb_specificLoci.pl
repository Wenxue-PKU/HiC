#!/usr/bin/perl
# 2015/07/02 filtered dbからの読み取りに変更
# 2015/03/28 dbファイルからデータを読み取って相互作用ファイルを作成する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;

use DBI;

if((@ARGV != 16 and @ARGV != 14 and @ARGV != 12) or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [detabase files] -o [output prefix] -r [resolution] -d [distance normalize file] -b [black list of fragment] -c [chromsoome] -s [start position] -e [end position]\n";
}


my %opt;
getopts("i:o:r:d:b:c:s:e:", \%opt);
my $FILE_database = $opt{i};
my $FILE_out = $opt{o};
my $Resolution = $opt{r};
my $FILE_distance = $opt{d};
my $FILE_black = $opt{b};
my $chromosome = $opt{c};
my $start = $opt{s};
my $end = $opt{e};

# 全データを入れる変数
my %data;


#---------------------------------------
# Distance curveを読み込む
#---------------------------------------
my %AVE;
if(defined $FILE_distance){
	my $fh_dis = IO::File->new($FILE_distance) or die "cannot open $FILE_distance: $!";
	while($_ = $fh_dis->getline()){
		s/\r?\n//;
		my ($chr1, $chr2, $d, $score, $probability) = split /\t/;
		$AVE{"$chr1\t$chr2\t$d"} = $probability;
	}
	$fh_dis->close();
}

#---------------------------------------
# fragmentのblack listを読み込む
#---------------------------------------
my %Black;
if(defined $FILE_black){
	my $fh_in = IO::File->new($FILE_black ) or die "cannot open $FILE_black: $!";
	while($_ = $fh_in->getline()){
		s/\r?\n//;
		my ($chr, $fragID) = split /\t/;
		$Black{"$chr\t$fragID"} = 1;
	}
	$fh_in->close();
}


my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");



#---------------------------------------
# collect information
#---------------------------------------
my $sth_data;

# intra-chromosomeのデータのみを取得する
$sth_data = $dbh->prepare("select chr1, start1, end1, fragNum1, chr2, start2, end2, fragNum2, score from fragment
	where chr1=chr2 and chr1='$chromosome' and start1 >= $start and start2 >= $start and end1 <= $end and end2 <= $end;");
$sth_data->execute();
while(my $ref = $sth_data->fetchrow_arrayref()){
	my ($chr1, $start1, $end1, $frag1, $chr2, $start2, $end2, $frag2, $score) = @$ref;

	# 隣同士のfragmentはカウントしない
	if(abs($frag1 - $frag2) < 2){
		next;
	}

	# fragmentがblack listに含まれていたら計算しない
	if(exists $Black{"$chr1\t$frag1"}){
		next;
	}
	if(exists $Black{"$chr2\t$frag2"}){
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


	my $distanceAve1 = 1;
	my $distanceAve2 = 1;
	my $distanceAve3 = 1;
	my $distanceAve4 = 1;

	if(defined $FILE_distance){
		# calculate distance average score
		my ($logDistance1, $logDistance2, $logDistance3, $logDistance4) = (-1, -1, -1, -1);

		my $distance1 = abs($start1 - $start2);
		my $distance2 = abs($start1 - $end2);
		my $distance3 = abs($end1 - $start2);
		my $distance4 = abs($end1 - $end2);

		$distance1 = int($distance1 / 100) * 100 + 50;
		$distance2 = int($distance2 / 100) * 100 + 50;
		$distance3 = int($distance3 / 100) * 100 + 50;
		$distance4 = int($distance4 / 100) * 100 + 50;

		if($distance1 < 50000){
			$logDistance1 = $distance1;
		}else{
			$logDistance1 = exp(sprintf("%.3f", log($distance1)));
			$logDistance1 = sprintf("%d", $logDistance1);
		}
		if($distance2 < 50000){
			$logDistance2 = $distance2;
		}else{
			$logDistance2 = exp(sprintf("%.3f", log($distance2)));
			$logDistance2 = sprintf("%d", $logDistance2);
		}
		if($distance3 < 50000){
			$logDistance3 = $distance3;
		}else{
			$logDistance3 = exp(sprintf("%.3f", log($distance3)));
			$logDistance3 = sprintf("%d", $logDistance3);
		}
		if($distance4 < 50000){
			$logDistance4 = $distance4;
		}else{
			$logDistance4 = exp(sprintf("%.3f", log($distance4)));
			$logDistance4 = sprintf("%d", $logDistance4);
		}

		unless(exists $AVE{"$chr1\t$chr2\t$logDistance1"}){
			die "$chr1\t$chr2\t$logDistance1\n";
		}
		unless(exists $AVE{"$chr1\t$chr2\t$logDistance2"}){
			die "$chr1\t$chr2\t$logDistance2\n";
		}
		unless(exists $AVE{"$chr1\t$chr2\t$logDistance3"}){
			die "$chr1\t$chr2\t$logDistance3\n";
		}
		unless(exists $AVE{"$chr1\t$chr2\t$logDistance4"}){
			die "$chr1\t$chr2\t$logDistance4\n";
		}

		$distanceAve1 = $AVE{"$chr1\t$chr2\t$logDistance1"};
		$distanceAve2 = $AVE{"$chr1\t$chr2\t$logDistance2"};
		$distanceAve3 = $AVE{"$chr1\t$chr2\t$logDistance3"};
		$distanceAve4 = $AVE{"$chr1\t$chr2\t$logDistance4"};
	}

	my $bin1a = int($start1/$Resolution) * $Resolution;
	my $bin1b = int($end1/$Resolution) * $Resolution;
	my $bin2a = int($start2/$Resolution) * $Resolution;
	my $bin2b = int($end2/$Resolution) * $Resolution;

	my $id1a = $chr1 . ":" . $bin1a;
	my $id1b = $chr1 . ":" . $bin1b;
	my $id2a = $chr2 . ":" . $bin2a;
	my $id2b = $chr2 . ":" . $bin2b;


	# 4つの組み合わせに均等に分配する
	$score = $score / 4;


	# count data（既に左側が小さいということは保証されている）
	$data{$id1a}{$id2a} += $score / $distanceAve1;
	$data{$id1a}{$id2b} += $score / $distanceAve2;
	$data{$id1b}{$id2a} += $score / $distanceAve3;
	$data{$id1b}{$id2b} += $score / $distanceAve4;

}
$sth_data->finish();
$dbh->disconnect();



#---------------------------------------
# output
#---------------------------------------
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";

# register bins
my @bins;
for(my $i = int($start/$Resolution) * $Resolution; $i <= int($end/$Resolution) * $Resolution; $i += $Resolution){
	push @bins, "$chromosome:$i";
	$fh_out->printf("\t$chromosome:$i:%d", $i + $Resolution - 1);
}
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
	my ($c, $m) = split /:/, $bins[$i];
	$fh_out->printf("%s:%d\t", $bins[$i], $m + $Resolution - 1);
	$fh_out->print(join("\t", @values) . "\n");
}
$fh_out->close();



