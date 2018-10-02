#!/bin/bash
# Prpare BUTLR format matrices

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-d, --directory [data directory]
		directory name of analysis file locate

	-n, --name [sample name]
		sample name

	-o, --log [log directory]
		log file directory

	-r, --resolution [resolution]
		resolution of output (ex. 40kb)

	-a, --genome [hg19|mm10]
		genome assembly version

EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hd:n:o:r:a:
LONG=help,directory:,name:,log:,resolution:,genome:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? -ne 0 ]]; then
    exit 2
fi
eval set -- "$PARSED"

while true; do
    case "$1" in
        -h|--help)
            get_usage
            exit 1
            ;;
        -d|--directory)
            DIR_DATA="$2"
            shift 2
            ;;
        -n|--name)
            NAME="$2"
            shift 2
            ;;
        -o|--log)
            DIR_LOG="$2"
            shift 2
            ;;
        -r|--resolution)
            RESOLUTION="$2"
            shift 2
            ;;
        -a|--genome)
            ASSEMBLY="$2"
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

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${DIR_LOG}" ] && echo "Please specify log directory" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "${ASSEMBLY}" ] && echo "Please specify genome assembly" && exit 1




case $ASSEMBLY in
	hg19) CHROM_SIZE=/wistar/noma/Data/Human_seq/hg19/LENGTH.txt
		  CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY) ;;
	mm10) CHROM_SIZE=/wistar/noma/Data/Mouse_seq/mm10/LENGTH.txt
		  CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY) ;;
	*)    echo "Please specify correct organism"
	      eixt 1 ;;
esac


### Create Temporary directory
[ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION}/tmpBUTLRT ] && mkdir ${DIR_DATA}/${NAME}/${RESOLUTION}/tmpBUTLRT

FILE_MAT_LIST=${DIR_DATA}/${NAME}/${RESOLUTION}/tmpBUTLRT/matrix.list


#==============================================================
# Convert matrices
#==============================================================
for i in $(seq 1 ${#CHRs[@]})
do
	let index=i-1
	CHR=${CHRs[index]}
	if [ -e ${DIR_DATA}/${NAME}/${RESOLUTION}/Raw/${CHR}.rds ]; then
		Rscript --vanilla --slave ${DIR_LIB}/ConvertMatricesForBUTLR.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/Raw/${CHR}.rds --chrom $CHROM_SIZE -o ${DIR_DATA}/${NAME}/${RESOLUTION}/tmpBUTLRT/${CHR}.matrix
		echo -e "$CHR\t${DIR_DATA}/${NAME}/${RESOLUTION}/tmpBUTLRT/${CHR}.matrix" >> $FILE_MAT_LIST
	else
		echo "${CHR}.rds (${NAME}:${RESOLUTION}) is not exists"
	fi
done

# convert to BUTLR
perl ${DIR_LIB}/matrixToButlr.pl -g $CHROM_SIZE -a $ASSEMBLY -m $FILE_MAT_LIST -r ${RESOLUTION} -o ${DIR_DATA}/${NAME}/${NAME}.${RESOLUTION}.btr


rm -r ${DIR_DATA}/${NAME}/${RESOLUTION}/tmpBUTLRT



