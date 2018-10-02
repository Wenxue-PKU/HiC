#!/usr/bin/perl
# 2017/09/27 各Directionのreadの数をカウントしてself-ligationとun-digestがどれくらいまで影響を及ぼすかを調べる

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

use DBI;

if(@ARGV != 4 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [data.db] -o [output file]\n";
}

my %opt;
getopts("i:o:", \%opt);
my $FILE_database = $opt{i};
my $FILE_out = $opt{o};


my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");

#---------------------------------------
# 合計値を計算する
#---------------------------------------
my %DATA;
my $sth = $dbh->prepare("select position1, direction1, position2, direction2 from map
	where uniq1 = 'U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and restLoc1 != 'NA'
	and restLoc2 != 'NA' and chr1=chr2 and abs(position1 - position2) < 1500000");
$sth->execute();
while(my $ref = $sth->fetchrow_arrayref()){
	my ($pos1, $direction1, $pos2, $direction2) = @$ref;
	my $combination = $direction1 . $direction2;
	my $distance = int(($pos2 - $pos1) / 100) * 100;
	$DATA{$distance}{$combination}++;
}
$sth->finish();


#---------------------------------------
# 結果を出力する
#---------------------------------------
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";
$fh_out->print(join("\t", "distance", "++", "+-", "-+", "--") . "\n");
for(my $d = 100; $d < 1000000; $d+=100){
	my @scores;
	foreach my $comb(qw(++ +- -+ --)){
		my $s = exists $DATA{$d}{$comb} ? $DATA{$d}{$comb} : 0;
		push @scores, $s;
	}
	$fh_out->print(join("\t", $d, @scores) . "\n");
}
$fh_out->close();






