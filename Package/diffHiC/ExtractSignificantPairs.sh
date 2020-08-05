#!/bin/bash
# Significantなペアを抽出する(diffHic packageを使う)


get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [sample names (separated by space)]

Fold change will be set2 / set1 (set1 is initial condition. set2 is changed condition)

Description
	-h, --help
		show help

	-v, --version
		show version

	-r, --resolution [resolution]
		resolution

	-g, --group [group numbers (1|2)]
		group setting. separated by , ex (1,1,2,2)

	-d, --data [directory]
		data directory
	
	-o, --out [output directory]
		output directory
	
	-x, --ref [ex. hg19]
		organism name

	--include [including chromsome list]
		if only specific chromosome should be calculated, specify the list. Separated by ,.

	--exclude [exclude chromsome list]
		list of excluding chromosomes. Separated by ,. Ex. chrM,chrY

	--FDR [FDR threshold]
		how to restrict the result. default (FDR < 0.05)
	
	--distance [distance to merge result]
		distance between interaction for merge results (default 0)
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvr:g:d:o:x:
LONG=help,version,resolution:,group:,data:,out:,ref:,include:,exclude:,FDR:,distance:
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
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
			;;
		-g|--group)
			group="$2"
			shift 2
			;;
		-d|--data)
			DIR_DATA="$2"
			shift 2
			;;
		-o|--out)
			DIR_OUT="$2"
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
		--FDR)
			FDR="$2"
			shift 2
			;;
		--distance)
			MERGE_DISTANCE="$2"
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

[ ! -n "${group}" ] && echo "Please specify group" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1
[ ! -n "${REF}" ] && echo "Please specify ref" && exit 1
FDR=${FDR:-0.05}
CHR_include=${CHR_include:-NA}
CHR_exclude=${CHR_exclude:-NA}
MERGE_DISTANCE=${MERGE_DISTANCE:-0}

NAME_LIST="$@"

[ ! -e ${DIR_OUT} ] && mkdir ${DIR_OUT}
cd ${DIR_OUT}
[ ! -e scores ] && mkdir log scores merge
UNIQ_ID=$(echo $DIR_OUT | rev | cut -c 1-12 | rev)

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

#==============================================================
# Significantly different fragment pairsを定義
#==============================================================
for i in $(seq 1 ${#CHRs[@]})
do
	let index=i-1
	CHR=${CHRs[index]}
	FILE_in=$(echo "$NAME_LIST" | xargs -n1 | xargs -I@ sh -c "echo ${DIR_DATA}/\@/${RESOLUTION}/Raw/${CHR}.rds" | xargs | tr ' ' ',')
	sbatch --account=nomalab -n 4 --job-name=si_${UNIQ_ID}_${CHR} --mem=100G --partition short -o "${DIR_OUT}/log/define_significant_pairs_${CHR}.log" --open-mode append --wrap="module load R/4.0.2 && Rscript --vanilla --slave ${DIR_LIB}/Extract_diff_pairs.R -i ${FILE_in} -o ${DIR_OUT}/scores/${CHR}.txt --group ${group}"
done

#==============================================================
# 指定した距離以内のinteractionをmergeする
#==============================================================
for i in $(seq 1 ${#CHRs[@]})
do
	let index=i-1
	CHR=${CHRs[index]}
	JOB_ID=($(squeue -o "%j %F" -u hidekit | grep -e "si_${UNIQ_ID}_${CHR}" | cut -f2 -d' ' | xargs))
	JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
	DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
	sbatch --account=nomalab -N 1 -n 1 --job-name=si2_${UNIQ_ID}_${CHR} $DEPEND --mem=10G --partition short -o "${DIR_OUT}/log/define_significant_pairs_${CHR}.log" --open-mode append --wrap="tail -n +2 ${DIR_OUT}/scores/${CHR}.txt | awk -v OFS='\t' -v T=$FDR '\$11<T{print}' | pgltools formatbedpe | pgltools merge -stdInA -c 11 -o min -d $MERGE_DISTANCE -noH > ${DIR_OUT}/merge/${CHR}.txt"
done


