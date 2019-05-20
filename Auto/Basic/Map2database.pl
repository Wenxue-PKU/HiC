#!/usr/bin/perl
# 2015/04/26 ソートして重複を除いてから登録する
# 2015/03/28 map fileをデータベースに登録する(ただし重複データは除く)

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;


use DBI;

if(@ARGV != 4  or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [sorted map file] -o [output databaseName]\n";
}

my %opt;
getopts("i:o:", \%opt);
my $FILE_sortedmap = $opt{i};
my $FILE_database = $opt{o};
my ($name, $dir, $ext) = &fileparse($FILE_sortedmap, '\..*');
my $TMP_PREFIX = $dir . 'tmp_register_' . $name;

#---------------------------------------
# db format
#---------------------------------------
# 1. id (read name)
# 2. chr 1
# 3. position 1
# 4. direction 1
# 5. map quality 1
# 6. restriction # 1
# 7. restriction side 1
# 8. restriction location 1
# 9. chr 2
# 10. position 2
# 11. direction 2
# 12. map quality 2
# 13. restriction # 1
# 14. restriction side 1
# 15. restriction location 1

#---------------------------------------
# create database
#---------------------------------------
my $dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");
my $ret = $dbh->do("create table map(
	id text,
	chr1 text,
	position1 integer,
	direction1 text,
	mapQ1 integer,
	restNum1 integer,
	restSide1 text,
	restLoc1 integer,
	uniq1 text,
	chr2 text,
	position2 integer,
	direction2 text,
	mapQ2 integer,
	restNum2 integer,
	restSide2 text,
	restLoc2 integer,
	uniq2 text,
	primary key(chr1, position1, direction1, chr2, position2, direction2));
");





#---------------------------------------
# register to database
#---------------------------------------
$dbh = DBI->connect("dbi:SQLite:dbname=$FILE_database");
my $q = join(",", ("?")x17);
my $sql = "insert into map values($q)";
my $sth = $dbh->prepare($sql);

my $FINISHED_FILE = 0;
my $REGISTERED_NUMBER = 0;
my $PREVIOUS_REPORT = 0;
my $previous = '';
my %duplication;

$dbh->do('BEGIN');
my $fh_map = IO::File->new($FILE_sortedmap) or die "cannot open $FILE_sortedmap : $!";
while($_ = $fh_map->getline()){
	s/\r?\n//;
	my ($id, $chr1, $loc1, $direction1, $mapQ1, $resID1, $resLoc1, $uniq1, $chr2, $loc2, $direction2, $mapQ2, $resID2, $resLoc2, $uniq2) = split /\t/;

	# 1万件ごとに一度commitし、100万件ごとに登録したデータの数をレポートする
	if($REGISTERED_NUMBER % 10000 == 0 and $REGISTERED_NUMBER != $PREVIOUS_REPORT){
		$dbh->do('COMMIT');
		$dbh->do('BEGIN');
		if($REGISTERED_NUMBER % 1000000 == 0){
			printf "%d were registered to database [%s]\n", $REGISTERED_NUMBER, &time();
			$PREVIOUS_REPORT = $REGISTERED_NUMBER;
		}
	}

	if($resID1 eq 'NA'){
		next;
	}
	if($resID2 eq 'NA'){
		next;
	}

	my $resNum1 = substr($resID1, 0, -1);
	my $resNum2 = substr($resID2, 0, -1);
	my $resSide1 = substr($resID1, -1, 1);
	my $resSide2 = substr($resID2, -1, 1);
	my @regiterData = ($id, $chr1, $loc1, $direction1, $mapQ1, $resNum1, $resSide1, $resLoc1, $uniq1, $chr2, $loc2, $direction2, $mapQ2, $resNum2, $resSide2, $resLoc2, $uniq2);

	if("$chr1\t$loc1\t$direction1" ne $previous){
		foreach my $key(keys %duplication){
			$sth->execute(@{$duplication{$key}});
			$REGISTERED_NUMBER++;
		}
		%duplication = ();
	}
	$duplication{"$chr2\t$loc2\t$direction2"} = \@regiterData;
	$previous = "$chr1\t$loc1\t$direction1";
}
foreach my $key(keys %duplication){
	$sth->execute(@{$duplication{$key}});
	$REGISTERED_NUMBER++;
}
$fh_map->close();
$dbh->do('COMMIT');
printf "Total %d of data were successfully registered\n", $REGISTERED_NUMBER;

# indexの作成
$sth->finish();
$dbh->do("create index chromosomeIndex on map(chr1, chr2)");
$dbh->disconnect();


sub time{
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime(time);
	my $time = sprintf ("%d/%0d %02d:%02d",$mon+1, $mday, $hour, $min);
	return $time;
}
