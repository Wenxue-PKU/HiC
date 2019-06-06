#!/bin/bash
# Modeling by MOGEN software

MOGEN_root=/wistar/noma/Program/MOGEN/

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-i, --in [input file name]
		input data name

	--matrix [input matrix]
		input matrix

	-d, --dir [output directory]
		output directory (Not OUTPUT_FOLDER parameter of MOGEN. This folder output all results)
	
	--NBR_OF_CHR [3]
		number of chromosomes

	--VERBOSE [true/false]
		information during optimization printed out
	
	--W3_inter [0.5]
		NEG_MAX_DIST

	--W3_intra [5.0]
		NEG_MAX_DIST

	--W4_inter [0.5]
		NEG_MIN_DIST
	
	--W4_intra [5.0]
		NEG_MIN_DIST
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hi:d:
LONG=help,in:,matrix:,dir:,NBR_OF_CHR:,VERBOSE:,W3_inter:,W3_intra:,W4_inter:,W4_intra:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [ $? -ne 0 ]; then
	exit 2
fi
eval set -- "$PARSED"

while true; do
	case "$1" in
		-h|--help)
			get_usage
			exit 1
			;;
		-i|--in)
			NAME_in="$2"
			shift 2
			;;
		--matrix)
			FILE_matrix="$2"
			shift 2
			;;
		-d|--dir)
			DIR_OUT="$2"
			shift 2
			;;
		--NBR_OF_CHR)
			NBR_OF_CHR="$2"
			shift 2
			;;
		--VERBOSE)
			VERBOSE="$2"
			shift 2
			;;
		--W3_inter)
			W3_inter="$2"
			shift 2
			;;
		--W3_intra)
			W3_intra="$2"
			shift 2
			;;
		--W4_inter)
			W4_inter="$2"
			shift 2
			;;
		--W4_intra)
			W4_intra="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "Programming error"
			exit 3
			;;
	esac
done

DIR_LIB=$(dirname $0)
TIME_STAMP=$(date +"%Y-%m-%d")

[ ! -n "${NAME_in}" ] && echo "Please specify input file name" && exit 1
[ ! -n "${FILE_matrix}" ] && echo "Please specify input matrix" && exit 1
[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1
NBR_OF_CHR=${NBR_OF_CHR:-3}
VERBOSE=${VERBOSE:-true}
W3_inter=${W3_inter:-0.5}
W3_intra=${W3_intra:-5.0}
W4_inter=${W4_inter:-0.5}
W4_intra=${W4_intra:-5.0}

[ ! -e ${DIR_OUT}/output ] && mkdir ${DIR_OUT}/output
[ ! -e ${DIR_OUT}/input ] && mkdir ${DIR_OUT}/input

#==============================================================
# convert matrix
#==============================================================
[ ! -e ${DIR_OUT}/input/${NAME_in}.txt ] && Rscript --vanilla --slave ${DIR_LIB}/Convert_matrix_for_MOGEN.R -i $FILE_matrix -o ${DIR_OUT}/input/${NAME_in}.txt --break_point ${DIR_OUT}/input/${NAME_in}_break_point.txt




#==============================================================
# make setting file
#==============================================================
cat<<EOF > ${DIR_OUT}/input/pos_max_dist_weight.txt
#the first coefficient is for all (chr1,chr2) where chr1 != chr2
1.0
#coefficients for intrachromosomal contacts of chromosomes
EOF

cat<<EOF > ${DIR_OUT}/input/pos_min_dist_weight.txt
#the first coefficient is for all (chr1,chr2) where chr1 != chr2
2.0
#coefficients for intrachromosomal contacts of chromosomes
EOF

### W3
cat<<EOF > ${DIR_OUT}/input/neg_max_dist_weight.txt
#the first coefficient is for all (chr1,chr2) where chr1 != chr2
$W3_inter
#coefficients for intrachromosomal contacts of chromosomes
EOF

### W4
cat<<EOF > ${DIR_OUT}/input/neg_min_dist_weight.txt
#the first coefficient is for all (chr1,chr2) where chr1 != chr2
$W4_inter
#coefficients for intrachromosomal contacts of chromosomes
EOF

for i in $(seq 1 $NBR_OF_CHR)
do
	echo "$i 1.0" >> ${DIR_OUT}/input/pos_max_dist_weight.txt
	echo "$i 4.0" >> ${DIR_OUT}/input/pos_min_dist_weight.txt
	echo "$i $W3_intra" >> ${DIR_OUT}/input/neg_max_dist_weight.txt
	echo "$i $W4_intra" >> ${DIR_OUT}/input/neg_min_dist_weight.txt
done


cat<<EOF > ${DIR_OUT}/input/parameter.txt
#all distances here are square distance, all 15 parameters are required and no space should be included in the value
#number of structures will be generated
NUM = 1

#number of chromosomes
NBR_OF_CHR = ${NBR_OF_CHR}

#file contains break-points for chromosomes
CHR_UPPER_BOUND_ID_FILE = ${DIR_OUT}/input/${NAME_in}_break_point.txt

#contact with interaction frequency less than this is considered as non-contact, 

#INTRA_IF_THRESHOLD = 0.65
INTRA_IF_THRESHOLD = 80%
#INTER_IF_THRESHOLD = 0.58
INTER_IF_THRESHOLD = 18%

#NOTICE: the following distances are in square
#maximum distance between 2 adjacent points
ADJACENT_DIST = 1.5
#contact distance, points that are in contact should have square distance less than this
#when it is large , the whole structure will be scaled down in optimization and zoom out later
CONTACT_DIST = 6.0
POS_MIN_DIST = 0.2
NEG_MAX_DIST_INTRA = 30
NEG_MAX_DIST_INTER = 150

#increase this parameter to improve contact score, (but will decrease non-contact score)
POS_MAX_DIST_WEIGHT_FILE = ${DIR_OUT}/input/pos_max_dist_weight.txt
#increase this parameter if adjacent points are to close to each other
POS_MIN_DIST_WEIGHT_FILE = ${DIR_OUT}/input/pos_min_dist_weight.txt
#increase this parameter to improve non-contact score, (but will decrease contact score)
NEG_MIN_DIST_WEIGHT_FILE = ${DIR_OUT}/input/neg_min_dist_weight.txt
#increase this parameter to prevent the structure from spanning too much (make the structure smaller)
NEG_MAX_DIST_WEIGHT_FILE = ${DIR_OUT}/input/neg_max_dist_weight.txt


OUTPUT_FOLDER = output
INPUT_FILE = input/${NAME_in}.txt

#set VERBOSE = true for information during optimization printed out
VERBOSE = ${VERBOSE}

#learning rate for the optimization process, increase the learning rate can speed up the optimization process significantly, but sometimes, the optimization may fail
#if the program fails to generate structures, or the distance between 2 consecutive points are too large, try to reduce this learning rate
LEARNING_RATE = 0.001
#during parameter adjustment, increase LEARNING_RATE and decrease MAX_ITERATION, so that "coarse" structures can be quickly generated
MAX_ITERATION = 20000
EOF


#==============================================================
# Run MOGEN
#==============================================================
cd ${DIR_OUT}
java -jar ${MOGEN_root}/bin/3DGenerator.jar ${DIR_OUT}/input/parameter.txt

