#!/bin/bash
# Output example Hi-C contact maps (すべてのchrosmosomeの全体図を作成してExcelに出力する)

get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [target sample name(s). separated by space]

Description
	-h, --help
		show help

	-v, --version
		show version
	
	-o, --out [output directory]
		output directory

	-r, --resolution [ex 200kb]
		hic map resolution

	-x, --organism [ex human]
		organism name

	-d, --data [data directory]
		data directory
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvo:r:x:d:
LONG=help,version,out:,resolution:,organism:,data:
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
		-o|--out)
			DIR_OUT="$2"
			shift 2
			;;
		-d|--data)
			DIR_DATA="$2"
			shift 2
			;;
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
			;;
		-x|--organism)
			ORGANISM="$2"
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

[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ $# -lt 1 ] && echo "Please specify target Hi-C sample name(s)" && exit 1

SAMPLES=$@


case $ORGANISM in
	pombe) CHRs=(I II III)
	       LENGTH=(5579133 4539804 2452883) ;;
	human) #CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
	       #LENGTH=(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560 59373566) ;;
		   CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
	       LENGTH=(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560) ;;
	human_EBV) CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY EBV)
			LENGTH=(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560 59373566 171823) ;;
	mouse) #CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
			#LENGTH=(195471971 182113224 160039680 156508116 151834684 149736546 145441459 129401213 124595110 130694993 122082543 120129022 120421639 124902244 104043685 98207768 94987271 90702639 61431566 171031299 91744698) ;;
			CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX)
			LENGTH=(195471971 182113224 160039680 156508116 151834684 149736546 145441459 129401213 124595110 130694993 122082543 120129022 120421639 124902244 104043685 98207768 94987271 90702639 61431566 171031299) ;;
	*) echo "Please specify correct organism"
		eixt 1 ;;
esac

[ ! -e ${DIR_OUT}/img ] && mkdir ${DIR_OUT}/img
[ ! -e ${DIR_OUT}/log ] && mkdir ${DIR_OUT}/log

for i in `seq 1 ${#CHRs[@]}`
do
	let index=i-1
	CHR=${CHRs[index]}
	START=1
	END=${LENGTH[index]}
	for NAME in $SAMPLES
	do
		sbatch -n 4 --job-name=dr_${NAME}_${CHR} $(sq --node) -o "${DIR_OUT}/log/${TIME_STAMP}_map_for_${NAME}_${CHR}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw_matrix.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds --normalize NA --zero NA --na na --chr ${CHR} --start ${START} --end ${END} --unit p --max 0.95 --color red --width 500 -o ${DIR_OUT}/img/${NAME}_${CHR}.png"
	done
done

### Pythonで出力
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "dr" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterany:${JOB_ID_string}"
sbatch -n 4 --job-name=sum $(sq --node) $DEPEND -o "${DIR_OUT}/log/${TIME_STAMP}_makeExcel.log" --open-mode append --wrap="python ${DIR_LIB}/summarize_draw_all_chromosome_result.py -o ${DIR_OUT}/Summary.xlsx --image ${DIR_OUT}/img --name $(echo $SAMPLES | tr ' ' ',') --chromosome $(IFS=,; echo \"${CHRs[*]}\")"

