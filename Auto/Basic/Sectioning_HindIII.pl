#!/usr/bin/perl
# 2014/03/05 changed for human
# 2012/11/15 Section to nearest HindIII recgnition site

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV != 2 or $ARGV[0] eq '--help'){
	die "Usage : $0  -i [Hind III sites file]\n";
}


my %opt;
getopts("i:", \%opt);
my $FILE_HINDIII = $opt{i};


### 20000bp bin size
my $unit = 20000;



### read Hind III sites
my %Hinds;
my %MAX_chr;
{
	my $fh_in = IO::File->new($FILE_HINDIII) or die "cannot open $FILE_HINDIII: $!";
	while($_ = $fh_in->getline()){
		if(m/^#/){
			next;
		}
		s/\r?\n//;
		my ($number, $chr, $pos, $before, $after) = split /\t/;
		my $cate = int($pos / $unit) * $unit;
		my $id = $chr . ':' . $number;
		push @{$Hinds{"$chr\t$cate"}}, $id;
		unless(exists $MAX_chr{$chr}){
			$MAX_chr{$chr} = 0;
		}
		if($MAX_chr{$chr} < $cate){
			$MAX_chr{$chr} = $cate;
		}
	}
	$fh_in->close();
}



### output HindIII number lists
### check the nearest Hind III for empty location
foreach my $chr(keys %MAX_chr){
	for(my $i = 0; $i <= $MAX_chr{$chr} + 5 * $unit; $i += $unit){
		my @lists;
		if(exists $Hinds{"$chr\t$i"}){
			push @lists, @{$Hinds{"$chr\t$i"}};
		}
		my $before = $i - $unit;
		if(exists $Hinds{"$chr\t$before"}){
			push @lists, @{$Hinds{"$chr\t$before"}};
		}
		my $after = $i + $unit;
		if(exists $Hinds{"$chr\t$after"}){
			push @lists, @{$Hinds{"$chr\t$after"}};
		}

		if(@lists > 0){
			print "$chr\t$i\t" . join(",", @lists) . "\n";
		}else{
			print "$chr\t$i\tNA\n";
		}
	}
}
