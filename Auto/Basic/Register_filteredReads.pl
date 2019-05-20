#!/usr/bin/perl
# 2016/12/14 スピードアップするために作成

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;


use DBI;

if(@ARGV != 4  or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [dataForFragDb.txt] -o [output databaseName]\n";
}

my %opt;
getopts("i:o:", \%opt);
my $FILE_import = $opt{i};
my $FILE_database = $opt{o};

#---------------------------------------
# create database
#---------------------------------------
my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");

# 高速化の設定
$dbh->do("PRAGMA journal_mode = WAL");
$dbh->do("PRAGMA synchronous = NORMAL");

my $ret = $dbh->do("create table fragment(
	chr1 text,
	start1 integer,
	end1 integer,
	fragNum1 integer,
	chr2 text,
	start2 integer,
	end2 integer,
	fragNum2 integer,
	score integer,
	primary key(chr1, fragNum1, chr2, fragNum2));
");



#---------------------------------------
# register to database
#---------------------------------------
my $q = join(",", ("?")x9);
my $sql = "insert into fragment values($q)";
my $sth = $dbh->prepare($sql);

my $REGISTERED_NUMBER = 0;
my $PREVIOUS_REPORT = 0;

$dbh->do('BEGIN');
my $fh_data = IO::File->new($FILE_import) or die "cannot open $FILE_import : $!";
while($_ = $fh_data->getline()){
	s/\r?\n//;
	my @data = split /\t/;

	# 1万件ごとに一度commitし、100万件ごとに登録したデータの数をレポートする
	if($REGISTERED_NUMBER % 10000 == 0 and $REGISTERED_NUMBER != $PREVIOUS_REPORT){
		$dbh->do('COMMIT');
		$dbh->do('BEGIN');
		if($REGISTERED_NUMBER % 1000000 == 0){
			printf "%d were registered to database [%s]\n", $REGISTERED_NUMBER, &time();
			$PREVIOUS_REPORT = $REGISTERED_NUMBER;
		}
	}
	$sth->execute(@data);
	$REGISTERED_NUMBER++;
}
$fh_data->close();
$dbh->do('COMMIT');
printf "Total %d of data were successfully registered\n", $REGISTERED_NUMBER;

# indexの作成
$sth->finish();
$dbh->do("create index positionIndex on fragment(chr1, start1, chr2, start2)");
$dbh->do("create index chromosomeIndex on fragment(chr1, chr2)");
$dbh->disconnect();


sub time{
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime(time);
	my $time = sprintf ("%d/%0d %02d:%02d",$mon+1, $mday, $hour, $min);
	return $time;
}
