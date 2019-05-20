#!/usr/bin/perl
# 2017/05/15 Localの相互作用スコアを定義する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

use DBI;

if(@ARGV != 10 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [database] -o [output file] -c [chromsoome] -s [start position] -e [end position]\n";
}

my %opt;
getopts("i:o:c:s:e:", \%opt);
my $FILE_database = $opt{i};
my $FILE_output = $opt{o};
my $chromosome = $opt{c};
my $start = $opt{s} - 10000;
my $end = $opt{e} + 10000;
my $Resolution = 10;

my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");
my $sth_data = $dbh->prepare("select position1, position2 from map
		where chr1=chr2 and chr1='$chromosome' and position1 between $start and $end and position2 between $start and $end and
		uniq1='U' and uniq2='U' and mapQ1 > 30 and mapQ2 > 30 and direction1 = direction2 and abs(position1 - position2) < 10000");
$dbh->do('BEGIN');
$sth_data->execute();

my %Data;
while(my $ref = $sth_data->fetchrow_arrayref()){
	my ($pos1, $pos2) = @$ref;
	my $bin1 = int($pos1 / $Resolution)*$Resolution;
	my $bin2 = int($pos2 / $Resolution)*$Resolution;

	if($bin1 == $bin2){
		next;
	}

	$Data{$bin1}++;
	$Data{$bin2}++;
}
$dbh->do('COMMIT');
$sth_data->finish();

$dbh->disconnect();


my $fh_out = IO::File->new($FILE_output, 'w') or die "cannot write $FILE_output: $!";
foreach my $pos(sort {$a <=> $b} keys %Data){
	my $start = $pos;
	my $end = $pos + $Resolution - 1;
	my $value = $Data{$pos};
	$fh_out->print("$chromosome\t$start\t$end\t$value\n");
}
$fh_out->close();

