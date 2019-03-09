#!/bin/bash
# 複数のHiCmapの相関を計算する

get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [sample names separated by space]

Description
	-h, --help
		show help

	-v, --version
		show version

	-o, --out [output directory]
		output directory

	-r, --resolution [ex 200kb]
		hic map resolution

	-c, --chromosome [ex chr1]
		comparison chromsoome name

	--distance [ex 10000000]
		calculation max distance default 10Mb (=10000000)

	-d, --data [data directory]
		data directory
	
	-s, --summary
		after initial run, do this option to make summary report
EOF

}

get_version(){
	echo "${0} version 1.0"
}
SHORT=hvo:r:c:d:s
LONG=help,version,out:,resolution:,distance:,chromosome:,data:,summary
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
		--distance)
			MAX_distance="$2"
			shift 2
			;;
		-c|--chromosome)
			CHR="$2"
			shift 2
			;;
		-s|--summary)
			FLAG_run=TRUE
			shift 1
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
[ ! -n "${CHR}" ] && echo "Please specify chromosome name" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ $# -lt 1 ] && echo "Please specify target Hi-C sample name(s)" && exit 1
MAX_distance=${MAX_distance:-10000000}
FLAG_run=${FLAG_run:-FALSE}

SAMPLES=($@)
SAMPLE_NUM=$#


[ ! -e ${DIR_OUT}/score ] && mkdir ${DIR_OUT}/score
[ ! -e ${DIR_OUT}/log ] && mkdir ${DIR_OUT}/log


if [ "${FLAG_run}" = "FALSE" ]; then
	let maxk=${SAMPLE_NUM}-2
	for k in $(seq 0 $maxk)
	do
		let sj=${k}+1
		let maxj=${SAMPLE_NUM}-1
		for j in $(seq $sj $maxj)
		do
			NAME1=${SAMPLES[$k]}
			NAME2=${SAMPLES[$j]}
			FILE_map1=${DIR_DATA}/${NAME1}/${RESOLUTION}/Raw/${CHR}.rds
			FILE_map2=${DIR_DATA}/${NAME2}/${RESOLUTION}/Raw/${CHR}.rds
			sbatch -n 4 --job-name=hicrep_${NAME1}_${NAME2}_${CHR} $(sq --node) -o "${DIR_OUT}/log/${TIME_STAMP}_correlation_for_${NAME1}_${NAME2}_${CHR}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Correlation_of_two_map_by_HiCRep.R -a $FILE_map1 -b $FILE_map2 --max_distance $MAX_distance > ${DIR_OUT}/score/${NAME1}_${NAME2}_${CHR}.txt"
		done
	done
fi

#==============================================================
# まとめる
#==============================================================
if [ "${FLAG_run}" = "TRUE" ]; then
	FILE_summary=${DIR_OUT}/Correlation_of_${CHR}.txt
	echo "${CHR} ${SAMPLES[@]}" | tr ' ' '\t' > $FILE_summary
	let max=${SAMPLE_NUM}-1
	for k in $(seq 0 $max)
	do
		NAME1=${SAMPLES[$k]}
		echo -n "${NAME1}" >> $FILE_summary
		for j in $(seq 0 $max)
		do
			NAME2=${SAMPLES[$j]}
			if [ $k -eq $j ]; then
				score=1
			elif [ $k -lt $j ]; then
				score=$(cat ${DIR_OUT}/score/${NAME1}_${NAME2}_${CHR}.txt)
			else
				score=$(cat ${DIR_OUT}/score/${NAME2}_${NAME1}_${CHR}.txt)
			fi
			echo -en "\t$score" >> $FILE_summary
		done
		echo >> $FILE_summary
	done
fi
