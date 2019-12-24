#!/bin/bash
# Extract significant associations

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	--max [max distance (ex. 2Mb)]
		default is 2Mb
	
	-d, --data [directory]
		data directory
	
	-o, --out [output directory]
		output directory
	
	-x, --organism [human|mouse]
		organism name
	
	--control [read number]
		minimum required read value for background (default is 1)
	
	--FDR [FDR value]
		FDR cut off. (default is 0.01)

	--local [local fc threshold]
		local fc threshold. (default is 4)

	--fc [fold-change threhold to background]
		fold-change threshold to same distance background. (default is 4)
	
	--background [threshold for background]
		background average of surrounded -200kb, +200kb area (default 4)

	--remove_tmp [TRUE|FALSE]
		remove temporary directory(TRUE) or not(FALSE). (Default is TRUE)
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvd:o:x:
LONG=help,version,max:,data:,out:,organism:,control:,FDR:,local:,fc:,background:,remove_tmp:
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
		--max)
			MAX_distance="$2"
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
		-x|--organism)
			ORGANISM="$2"
			shift 2
			;;
		--control)
			T_CONTROL="$2"
			shift 2
			;;
		--FDR)
			FDR="$2"
			shift 2
			;;
		--local)
			T_LOCAL="$2"
			shift 2
			;;
		--fc)
			T_fc="$2"
			shift 2
			;;
		--background)
			T_back="$2"
			shift 2
			;;
		--remove_tmp)
			FLAG_remove_tmp="$2"
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

[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
MAX_distance=${MAX_distance:-"2Mb"}
MAX_distance=${MAX_distance/Mb/000000}
MAX_distance=${MAX_distance/kb/000}
T_CONTROL=${T_CONTROL:-1}
FDR=${FDR:-0.01}
T_LOCAL=${T_CONTROL:-4}
T_fc=${T_fc:-4}
T_back=${T_back:-4}
FLAG_remove_tmp=${FLAG_remove_tmp:-TRUE}

case $ORGANISM in
	human) CHRs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX" ;;
	mouse) CHRs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX" ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac

[ ! -e ${DIR_OUT} ] && mkdir ${DIR_OUT}
[ ! -e ${DIR_OUT}/log ] && mkdir ${DIR_OUT}/log
[ ! -e ${DIR_OUT}/scores ] && mkdir ${DIR_OUT}/scores
UNIQ_ID=$(echo $DIR_OUT | rev | cut -c 1-20 | rev)

#==============================================================
# Significant fragment pairsを定義
#==============================================================
for CHR in $CHRs
do
	sbatch --account=nomalab -n 1 -N 1 --job-name=si_${UNIQ_ID}_${CHR} --mem=40G $(sq --node --partition short) -o "${DIR_OUT}/log/define_significant_pairs_${CHR}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Define_all_significant_pairs.R -i ${DIR_DATA}/${CHR}.rds -o ${DIR_OUT}/scores/${CHR}_sig.txt --all ${DIR_OUT}/scores/${CHR}_all.txt --max ${MAX_distance} --control $T_CONTROL --FDR $FDR --local $T_LOCAL --fc $T_fc --background $T_back"
done

### 結果をまとめてlocal fold-changeでソートする
JOB_ID=($(squeue -o "%j %F" -u hidekit | grep -e "si_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch --account=nomalab -n 1 --job-name=si2_${UNIQ_ID} $DEPEND --mem=40G $(sq --node --partition short) -o "${DIR_OUT}/merge.log" --open-mode append <<-EOF
#!/bin/sh
cd ${DIR_OUT}/scores
cat chr1_sig.txt | head -n1 > ${DIR_OUT}/significant.txt
ls chr*_sig.txt | xargs -n1 | xargs -n1 -I@ sh -c "cat \@ | tail -n+2" | sort -k13,13r >> ${DIR_OUT}/significant.txt
[ "$FLAG_remove_tmp" = "TRUE" ] && rm -rf ${DIR_OUT}/scores && rm -rf ${DIR_OUT}/log
EOF


