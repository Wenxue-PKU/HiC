#!/usr/bin/perl
# 2017/10/04 intra-chromosome以外のbinレベルでのreadの合計値を計算する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;

use DBI;

if((@ARGV != 12 and @ARGV != 10) or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [detabase files xxx_fragment.db] -o [output prefix] -r [resolution] -b [black list of fragment] -m [chromosome1] -n [chromosome2]\n";
}

my %opt;
getopts("i:o:r:b:m:n:", \%opt);
my $FILE_database = $opt{i};
my $FILE_out_prefix = $opt{o};
my $Resolution = $opt{r};
my $FILE_black = $opt{b};
my $chromosome1 = $opt{m};
my $chromosome2 = $opt{n};


if($chromosome1 eq $chromosome2){
	die "Chromosome should be different between two chromosome\n";
}
if(&ComparisonChr($chromosome1, $chromosome2) == 1){
	($chromosome1, $chromosome2) = ($chromosome2, $chromosome1);
}
my @chromosomes = ($chromosome1, $chromosome2);

# 全データを入れる変数
my %data;


#---------------------------------------
# fragmentのblack listを読み込む
#---------------------------------------
my %Black;
if(defined $FILE_black){
	my $fh_in = IO::File->new($FILE_black ) or die "cannot open $FILE_black: $!";
	$fh_in->getline();
	while($_ = $fh_in->getline()){
		s/\r?\n//;
		my ($chr, $fragID) = split /\t/;
		$Black{"$chr\t$fragID"} = 1;
	}
	$fh_in->close();
}


my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");


#---------------------------------------
# check max length
#---------------------------------------
my %MAX_chr;
my $sth_maxCheck1 = $dbh->prepare("select max(end1) from fragment where chr1=?");
my $sth_maxCheck2 = $dbh->prepare("select max(end2) from fragment where chr2=?");
foreach my $chr(@chromosomes){
	$sth_maxCheck1->execute($chr);
	my ($m1) = $sth_maxCheck1->fetchrow_array();
	$sth_maxCheck2->execute($chr);
	my ($m2) = $sth_maxCheck2->fetchrow_array();
	$MAX_chr{$chr} = $m1 < $m2 ? $m2 : $m1;
}
$sth_maxCheck1->finish();
$sth_maxCheck2->finish();

#---------------------------------------
# collect information
#---------------------------------------
# inter-chromosomeのデータのみを取得する
my $sth_data = $dbh->prepare("select chr1, start1, end1, fragNum1, chr2, start2, end2, fragNum2, score from fragment where chr1='$chromosome1' and chr2='$chromosome2';");
$sth_data->execute();
while(my $ref = $sth_data->fetchrow_arrayref()){
	my ($chr1, $start1, $end1, $frag1, $chr2, $start2, $end2, $frag2, $score) = @$ref;

	# fragmentがblack listに含まれていたら計算しない
	if(exists $Black{"$chr1\t$frag1"}){
		next;
	}
	if(exists $Black{"$chr2\t$frag2"}){
		next;
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

	# count data
	$data{$id1a}{$id2a} += $score;
	$data{$id1a}{$id2b} += $score;
	$data{$id1b}{$id2a} += $score;
	$data{$id1b}{$id2b} += $score;
}
$sth_data->finish();
$dbh->disconnect();



#---------------------------------------
# output
#---------------------------------------
foreach my $chr(@chromosomes){
	my $FILE_out = $FILE_out_prefix . $chromosome1 . "_" . $chromosome2 . ".matrix";
	my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";

	# register bins
	my @bins1;
	for(my $i = 0; $i < $MAX_chr{$chromosome1}; $i += $Resolution){
		push @bins1, "$chromosome1:$i";
	}

	my @bins2;
	for(my $i = 0; $i < $MAX_chr{$chromosome2}; $i += $Resolution){
		push @bins2, "$chromosome2:$i";
		$fh_out->printf("\t$chromosome2:$i:%d", $i + $Resolution - 1);
	}
	$fh_out->print("\n");



	for(my $i = 0; $i < @bins1; $i++){
		my @values;
		for(my $j = 0; $j < @bins2; $j++){
			my $value = exists $data{$bins1[$i]}{$bins2[$j]} ? $data{$bins1[$i]}{$bins2[$j]} : 0;
			push @values, $value;
		}
		my ($c, $m) = split /:/, $bins1[$i];
		$fh_out->printf("%s:%d\t", $bins1[$i], $m + $Resolution - 1);
		$fh_out->print(join("\t", @values) . "\n");
	}
	$fh_out->close();


}

sub ComparisonChr{
	my ($chr1, $chr2) = @_;
	if($chr1 eq 'I'){
		$chr1 = 1;
	}
	if($chr2 eq 'I'){
		$chr2 = 1;
	}
	if($chr1 eq 'II'){
		$chr1 = 2;
	}
	if($chr2 eq 'II'){
		$chr2 = 2;
	}
	if($chr1 eq 'III'){
		$chr1 = 3;
	}
	if($chr1 eq 'III'){
		$chr1 = 3;
	}

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

	if($num1 < $num2){
		return -1;
	}else{
		return 1;
	}
}
