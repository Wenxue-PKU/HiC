#!/usr/bin/perl
# 2015/06/18 repeatのフラグを新しく作る
# 2015/04/18 bowtie2用にいくつか修正
# 2014/11/04 modify for HiC of S.pombe
# 2014/03/05 changed for human
# 2013/09/11 update for using latest s.pombe sequence
# 2012/12/05 Bug fix. Location was differently reported
# 2012/11/06 Alignment by BWA with several information for HiD or ChIA-PET data
# 2012/10/27 Alignment by BWA without allowing repeat location for HiD (ChIA-PET) data

####
#### cation !!!!!!  (aligned location was shift to HindIII location. Not original place)
####



use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;


#---------------------------------------
# FIX parameter
#---------------------------------------
my $unit = 20000;	# section file bin size (20000bp)

if(@ARGV != 8  or $ARGV[0] eq '--help'){
	die "Usage : $0 -a [sam file 1] -b [sam file 2] -o [output map file] -x [organism]\n";
}


my %opt;
getopts("a:b:o:x:", \%opt);
my $FILE_sam1 = $opt{a};
my $FILE_sam2 = $opt{b};
my $FILE_map = $opt{o};
my $ORGANISM = $opt{x};

my $FILE_HINDIII_index;
my $FILE_HINDIII_def;
if($ORGANISM eq 'human'){
	$FILE_HINDIII_index = '/wistar/noma/Data/Human_seq/hg19/Sectioning_MboI.txt';
	$FILE_HINDIII_def = '/wistar/noma/Data/Human_seq/hg19/MboI_sites.txt';
}elsif($ORGANISM eq 'human_EBV'){
	$FILE_HINDIII_index = '/wistar/noma/Data/Human_seq/hg19_EBV/Sectioning_MboI.txt';
	$FILE_HINDIII_def = '/wistar/noma/Data/Human_seq/hg19_EBV/MboI_sites.txt';
}elsif($ORGANISM eq 'human_HindIII'){
	$FILE_HINDIII_index = '/wistar/noma/Data/Human_seq/hg19/Sectioning_HindIII.txt';
	$FILE_HINDIII_def = '/wistar/noma/Data/Human_seq/hg19/HindIII_sites.txt';
}elsif($ORGANISM eq 'mouse'){
	$FILE_HINDIII_index = '/wistar/noma/Data/Mouse_seq/mm10/Sectioning_MboI.txt';
	$FILE_HINDIII_def = '/wistar/noma/Data/Mouse_seq/mm10/MboI_sites.txt';
}elsif($ORGANISM eq 'pombe'){
	$FILE_HINDIII_def = '/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18/MboI_sites.txt';
	$FILE_HINDIII_index = '/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18/Sectioning_MboI.txt';
}elsif($ORGANISM eq 'pombe_HindIII'){
	$FILE_HINDIII_def = '/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18/HindIII_sites.txt';
	$FILE_HINDIII_index = '/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v1.18/Sectioning_HindIII.txt';
}else{
	die "please specify correct organism name\n";
}


#---------------------------------------
# SAM format
#---------------------------------------
#@SQ	SN:chr1	LN:5579133
#@SQ	SN:chr2	LN:4539804
#@SQ	SN:chr3	LN:2452883
#2_ChIA-PET_linkerAA	0	chr2	1896903	37	25M	*	0	0	ATTAGCCAAGCGGTTTCATGATATC	@DDDBABDFHHEEIIIAHI@9<BHB	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:25
#4_ChIA-PET_linkerAA	16	chr3	2789	0	25M	*	0	0	ATATCCCTCCTCTCTCCACCAAAAG	B??::@FDAFBGHFDFDADA?BB?:	XT:A:R	NM:i:0	X0:i:5	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:25
#5_ChIA-PET_linkerAA	16	chr3	1233976	37	25M	*	0	0	TCCTTTCGAGTATCAACCAGCTTCA	?FCA3AC3IC>BBIFFDF<AA++A1	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:25
#6_ChIA-PET_linkerAA	16	chr3	16411	0	25M	*	0	0	AATATTTTTTTTACTAGGATTTGTG	HBF@IIGIIIIGEIFDFHFDDDDA:	XT:A:R	NM:i:0	X0:i:3	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:25
#7_ChIA-PET_linkerAA	0	chr2	1601563	0	25M	*	0	0	TTTCGATGACTTCCCTAATAAATTA	@FFFDFFFGDHGG>FHIEEHIIII@	XT:A:R	NM:i:0	X0:i:2	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:25	XA:Z:chr2,+1646901,25M,0;


#---------------------------------------
# map format
#---------------------------------------
# 1. id
# 2. chr 1
# 3. position 1
# 4. direction 1
# 5. map quality 1
# 6. hindiii id 1
# 7. hindiii location 1
# 8. chr 2
# 9. position 2
# 10. direction 2
# 11. map quality 2
# 12. hindiii id 2
# 13. hindiii location 2

#---------------------------------------
# read Hind III file
#---------------------------------------
### read Hind III sites
my %Hinds;
{
	my $fh_in = IO::File->new($FILE_HINDIII_def) or die "cannot open $FILE_HINDIII_def: $!";
	while($_ = $fh_in->getline()){
		if(m/^#/){
			next;
		}
		s/\r?\n//;
		my ($number, $chr, $pos, $before, $after) = split /\t/;
		my $id = $chr . ':' . $number;
		$Hinds{$id} = $pos;
	}
	$fh_in->close();
}
### read Hind III index file
my %Hind_index;
{
	my $fh_in = IO::File->new($FILE_HINDIII_index) or die "cannot open $FILE_HINDIII_index: $!";
	while($_ = $fh_in->getline()){
		if(m/^#/){
			next;
		}
		s/\r?\n//;
		my ($chr, $loc, $ids) = split /\t/;
		if($ids ne 'NA'){
			my @lists = split /,/, $ids;
			$Hind_index{$chr}{$loc} = \@lists;
		}
	}
	$fh_in->close();
}



#---------------------------------------
# mapping two files
#---------------------------------------
my $TOTAL_read = 0;
my $Aligned_single = 0;
my $Aligned_both = 0;
my $Not_aligned = 0;

my $fh_sam1 = IO::File->new($FILE_sam1) or die "cannot open $FILE_sam1: $!";
my $fh_sam2 = IO::File->new($FILE_sam2) or die "cannot open $FILE_sam2: $!";
my $fh_map = IO::File->new($FILE_map, 'w') or die "cannot write $FILE_map: $!";
while(my $sam1 = $fh_sam1->getline()){
	my $sam2 = $fh_sam2->getline();
#	if($sam1 =~ m/^\@/){
#		next;
#	}
#	$sam1 =~ s/\r?\n//;
#	$sam2 =~ s/\r?\n//;
	my ($id1, $flag1, $chr1, $position1, $mapQ1, $CIAGR1, $mate_a1, $mate_b1, $mate_c1, $seq1, $quality1, @option1) = split /\t/, $sam1;
	my ($id2, $flag2, $chr2, $position2, $mapQ2, $CIAGR2, $mate_a2, $mate_b2, $mate_c2, $seq2, $quality2, @option2) = split /\t/, $sam2;


	$id1 = substr($id1, 0, -2);
	$id2 = substr($id2, 0, -2);


	if($id1 ne $id2){
		die ("$id1 and $id2 is different\n");
	}

	$TOTAL_read++;
	my $flag_aligned = 0;
	my ($direction1, $direction2) = ('NA','NA');

	if($flag1 == 0){
		$direction1 = '+';
		$flag_aligned++;
	}elsif($flag1 == 16){
		$direction1 = '-';
		$position1 += length($seq1) - 1;
		$flag_aligned++;
	}else{
		$chr1 = 'NA';
		$position1 = 'NA';
	}



	if($flag2 == 0){
		$direction2 = '+';
		$flag_aligned++;
	}elsif($flag2 == 16){
		$direction2 = '-';
		$position2 += length($seq2) - 1;
		$flag_aligned++;
	}else{
		$chr2 = 'NA';
		$position2 = 'NA';
	}

	my ($hinID1, $hinLoc1) = ('NA', 'NA');
	my ($hinID2, $hinLoc2) = ('NA', 'NA');
	if($chr1 ne 'NA'){
		($hinID1, $hinLoc1) = &FindHindIIIinfo($chr1, $position1, $direction1);
	}
	if($chr2 ne 'NA'){
		($hinID2, $hinLoc2) = &FindHindIIIinfo($chr2, $position2, $direction2);
	}


	my $Uniq1 = "U";
	foreach my $v(@option1){
		my $f = index $v, "XS:i";
		if($f != -1){
			$Uniq1 = "R";
			last;
		}
	}
	my $Uniq2 = "U";
	foreach my $v(@option2){
		my $f = index $v, "XS:i";
		if($f != -1){
			$Uniq2 = "R";
			last;
		}
	}


	if($flag_aligned == 0){
		$Not_aligned++;
	}elsif($flag_aligned == 1){
		$Aligned_single++;
	}elsif($flag_aligned == 2){
		$Aligned_both++;
	}

	# idとしてread番号を割り当てる
	$fh_map->print("$TOTAL_read\t$chr1\t$position1\t$direction1\t$mapQ1\t$hinID1\t$hinLoc1\t$Uniq1\t$chr2\t$position2\t$direction2\t$mapQ2\t$hinID2\t$hinLoc2\t$Uniq2\n");
}
$fh_sam1->close();
$fh_sam2->close();
$fh_map->close();

print "OutputFile : $FILE_map\n";
printf "Total read:\t%d\n", $TOTAL_read;
printf "Both aligned:\t%d\n", $Aligned_both;
printf "Single aligned:\t%d\n", $Aligned_single;
printf "Not aligned:\t%d\n", $Not_aligned;



sub FindHindIIIinfo{
	my ($chr, $pos, $direction) = @_;
	my $cate = int($pos / $unit) * $unit;
	if(exists $Hind_index{$chr}{$cate}){
		my @candidates = @{$Hind_index{$chr}{$cate}};

		# 左側の候補と右側の候補を２つ探す
		my $MIN_left = 99999999999;
		my $MIN_right = 99999999999;
		my $MIN_left_id = '';
		my $MIN_right_id = '';
		foreach my $i(@candidates){
			if($pos < $Hinds{$i}){
				my $distance = $Hinds{$i} - $pos;
				if($distance < $MIN_left){
					$MIN_left = $distance;
					$MIN_left_id = $i;
				}
			}else{
				my $distance = $pos - $Hinds{$i};
				if($distance < $MIN_right){
					$MIN_right = $distance;
					$MIN_right_id = $i;
				}
			}
		}

		# もし向きが+なら右の候補を採用し、-なら左の候補を採用する
		my $MIN_id = '';
		if($direction eq '+'){
			$MIN_id = $MIN_right_id;
		}else{
			$MIN_id = $MIN_left_id;
		}

		if($MIN_id eq ''){
			return ('NA', 'NA');
		}

		my ($chrTmp, $hinID) = split /:/, $MIN_id;
		if($direction eq '+'){
			$hinID .= 'L';
		}else{
			$hinID .= 'R';
		}

		return ($hinID, $Hinds{$MIN_id});
	}else{
#		warn ("Nearest restriction cites were not found for $chr : $pos \n");
		return ('NA', 'NA');
	}}


