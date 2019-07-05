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
	
	-n, --name
		sample name

	-r, --resolution [resolution (ex. 10kb)]
		resolution (default 10kb)

	--max [max distance (ex. 2Mb)]
		default is 2Mb
	
	-d, --data [directory]
		data directory
	
	-o, --out [output file]
		output file
	
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

	--remove_tmp [TRUE|FALSE]
		remove temporary directory(TRUE) or not(FALSE). (Default is TRUE)
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvn:r:d:o:x:
LONG=help,version,name:,resolution:,max:,data:,out:,organism:,control:,FDR:,local:,fc:,remove_tmp:
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
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
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
			FILE_OUT="$2"
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

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output file" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
MAX_distance=${MAX_distance:-"2Mb"}
MAX_distance=${MAX_distance/Mb/000000}
MAX_distance=${MAX_distance/kb/000}
T_CONTROL=${T_CONTROL:-1}
FDR=${FDR:-0.01}
T_LOCAL=${T_CONTROL:-4}
T_fc=${T_fc:-4}
RESOLUTION=${RESOLUTION:-10kb}
FLAG_remove_tmp=${FLAG_remove_tmp:-TRUE}

case $ORGANISM in
	human) CHRs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX" ;;
	mouse) CHRs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX" ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac

DIR_tmp=${FILE_OUT}_tmpDir
[ ! -e ${DIR_tmp} ] && mkdir ${DIR_tmp} && mkdir ${DIR_tmp}/log ${DIR_tmp}/scores
UNIQ_ID=$(echo $FILE_OUT | rev | cut -c 1-20 | rev)
FILE_log=${FILE_OUT/.txt/.log}

#==============================================================
# Significant fragment pairsを定義
#==============================================================
for CHR in $CHRs
do
	FILE_in=${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds
	sbatch -n 4 -N 2 --job-name=si_${UNIQ_ID}_${CHR} $(sq --node) -o "${DIR_tmp}/log/define_significant_pairs_${CHR}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Define_all_significant_pairs.R -i ${FILE_in} -o ${DIR_tmp}/scores/${CHR}.txt --max ${MAX_distance} --control $T_CONTROL --FDR $FDR --local $T_LOCAL --fc $T_fc"
done

### 結果をまとめてFDRでソートする
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "si_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=si2_${UNIQ_ID} $DEPEND $(sq --node) -o "${FILE_log}" --open-mode append <<-EOF
#!/bin/sh
cd ${DIR_tmp}/scores
cat chr1.txt | head -n1 > $FILE_OUT
ls chr*.txt | xargs -n1 | xargs -n1 -I@ sh -c "cat \@ | tail -n+2" | sort -k12,12g >> $FILE_OUT
ls ${DIR_tmp}/log/define_significant_pairs_*.log | xargs -n1 | xargs -n1 -I@ sh -c "cat @ | cut -d':' -f2 | tr -d ' ' | xargs" | awk -v OFS='\t' -v FDR=$FDR -v con=$T_CONTROL -v fc=$T_fc 'BEGIN{N1=0; N2=0; N3=0; N4=0; N5=0; N6=0}{N1+=\$1; N2+=\$2; N3+=\$3; N4+=\$4; N5+=\$5; N6+=\$6; N7+=\$7;}END{print "Total target: "N1; print "FDR < "FDR": "N2; print "Fold-change > "fc": "N3; print "Control > "con": "N4; print "Surround area check: "N5; print "After 1st filtering: "N6; print "After scan surrounded area: "N7;}'
[ "$FLAG_remove_tmp" = "TRUE" ] && rm -rf ${DIR_tmp}
EOF


