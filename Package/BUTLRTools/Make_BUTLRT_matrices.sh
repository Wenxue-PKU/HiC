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

	-r, --resolution [resolution]
		resolution of output (ex. 40kb)
	
	-t, --type [matrices type]
		Raw, ICE, ICE2 etc...

	-x, --ref [ex. hg19]
		organism reference sequence name

	--include [including chromsome list]
		if only specific chromosome should be calculated, specify the list. Separated by ,.

	--exclude [exclude chromsome list]
		list of excluding chromosomes. Separated by ,. Ex. chrM,chrY
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hd:n:r:t:x:
LONG=help,directory:,name:,resolution:,type:,ref:,include:,exclude:
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
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
			;;
		-t|--type)
			TYPE="$2"
			shift 2
			;;
		-x|--ref)
			REF="$2"
			shift 2
			;;
		--include)
			CHR_include="$2"
			shift 2
			;;
		--exclude)
			CHR_exclude="$2"
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
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "${REF}" ] && echo "Please specify ref" && exit 1
[ ! -n "${TYPE}" ] && echo "Please specify matrices type (Raw, ICE etc...)" && exit 1
CHR_include=${CHR_include:-NA}
CHR_exclude=${CHR_exclude:-NA}

#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/../../utils/load_setting.sh -x $REF -r NA

#-----------------------------------------------
# Load chromosome length
#-----------------------------------------------
CHR_TABLE=$(Rscript --vanilla --slave ${DIR_LIB}/../../utils/Chromosome_length.R --in $FILE_CHROME_LENGTH --include $CHR_include --exclude $CHR_exclude)
CHRs=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==1' | tr ',' ' '))
LENGTH=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==2' | tr ',' ' '))
CHRs_list=$(echo ${CHRs[@]} | tr ' ' ',')


### Create Temporary directory
DIR_tmp=$(mktemp -d /tmp/tmp_${NAME}.XXXXX)
FILE_MAT_LIST=${DIR_tmp}/matrix.list
trap "rm -r ${DIR_tmp}" 0

#==============================================================
# Convert matrices
#==============================================================
for i in $(seq 1 ${#CHRs[@]})
do
	let index=i-1
	CHR=${CHRs[index]}
	if [ -e ${DIR_DATA}/${NAME}/${RESOLUTION}/${TYPE}/${CHR}.rds ]; then
		Rscript --vanilla --slave ${DIR_LIB}/ConvertMatricesForBUTLR.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/${TYPE}/${CHR}.rds --chrom $FILE_CHROME_LENGTH -o ${DIR_tmp}/${CHR}.matrix
		echo -e "$CHR\t${DIR_tmp}/${CHR}.matrix" >> $FILE_MAT_LIST
	else
		echo "${CHR}.rds (${NAME}:${RESOLUTION}) is not exists"
	fi
done

# convert to BUTLR
[ ! -e ${DIR_DATA}/${NAME}/BUTLRT ] && mkdir ${DIR_DATA}/${NAME}/BUTLRT
perl ${DIR_LIB}/matrixToButlr.pl -g $FILE_CHROME_LENGTH -a $REF -m $FILE_MAT_LIST -r ${RESOLUTION} -o ${DIR_DATA}/${NAME}/BUTLRT/${NAME}.${RESOLUTION}.btr

exit 0



