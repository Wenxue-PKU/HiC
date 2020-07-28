#!/bin/bash
# Extracting score


get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-n, --name
		job name

	-l, --log
		log file

	-m, --matrix [matrix directory]
		matrix directory

	-i, --in [target combination file template]
		target combination file ex. window10kb_combination_XXX.txt  replace XXX with chromosome name
	
	-o, --out [output file template]
		output file template. replace XXX with chromosome name

	-x, --ref [ex. hg19]
		organism name
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvn:l:m:i:o:x:
LONG=help,version,name:,log:,matrix:,in:,out:,ref:
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
		-v|--version)
			get_version
			exit 1
			;;
		-n|--name)
			NAME="$2"
			shift 2
			;;
		-l|--log)
			FILE_log="$2"
			shift 2
			;;
		-m|--matrix)
			DIR_data="$2"
			shift 2
			;;
		-i|--in)
			TEMPLATE_IN="$2"
			shift 2
			;;
		-o|--out)
			TEMPLATE_OUT="$2"
			shift 2
			;;
		-x|--ref)
			REF="$2"
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
TIME_STAMP=$(date +"%Y-%m-%d_%H.%M.%S")
INPUT_FILES=$@

[ ! -n "${NAME}" ] && echo "Please specify name" && exit 1
[ ! -n "${FILE_log}" ] && echo "Please specify log file" && exit 1
[ ! -n "${DIR_data}" ] && echo "Please specify matrix directory" && exit 1
[ ! -n "${TEMPLATE_IN}" ] && echo "Please specify target file template" && exit 1
[ ! -n "${TEMPLATE_OUT}" ] && echo "Please specify output file template" && exit 1
[ ! -n "${REF}" ] && echo "Please specify ref" && exit 1

#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/../utils/load_setting.sh -x $REF -r NA


#-----------------------------------------------
# Load chromosome length
#-----------------------------------------------
CHR_TABLE=$(Rscript --vanilla --slave ${DIR_LIB}/../utils/Chromosome_length.R --in $FILE_CHROME_LENGTH --exclude "chrM,chrY")
CHRs=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==1' | tr ',' ' '))
LENGTH=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==2' | tr ',' ' '))
CHRs_list=$(echo ${CHRs[@]} | tr ' ' ',')

#==============================================================
# Raw matrixの作成
#==============================================================
for i in $(seq 1 ${#CHRs[@]})
do
	let index=i-1
	CHR=${CHRs[index]}
	FILE_IN=${TEMPLATE_IN/XXX/$CHR}
	FILE_OUT=${TEMPLATE_OUT/XXX/$CHR}
	sbatch --account=nomalab -N 1 -n 2 --job-name=exScore_${NAME}_${CHR} --mem=100G --partition short -o "${FILE_log}_${CHR}.tmp.log" --open-mode append --wrap="echo ${CHR}; Rscript --slave --vanilla ${DIR_LIB}/Extract_targetScore_with_distanceNormalized_Zscore.R -i $FILE_IN -m ${DIR_data}/${CHR}.rds -o $FILE_OUT --format rds"
done

JOB_ID=($(squeue -o "%j %F" -u hidekit | grep -e "exScore_${NAME}_${CHR}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterany:${JOB_ID_string}"
sbatch --account=nomalab -N 1 -n 2 --job-name=exScore_${NAME}_mix --mem=10G --partition short -o "${FILE_log}" --open-mode append --wrap="cat ${FILE_log}_*.tmp.log && rm ${FILE_log}_*.tmp.log"

