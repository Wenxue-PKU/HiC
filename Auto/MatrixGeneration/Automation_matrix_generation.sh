#!/bin/bash
# 2017/02/17 全自動 Hi-C matrix 作成プログラム

#==============================================================
# コマンド引数の受取り
#==============================================================
SHORT=d:n:x:r:o:t:e:c:
LONG=directory:,name:,organism:,resolution:,log:,intra:,normalization:,raw:
PARSED=`getopt --options $SHORT --longoptions $LONG --name "$0" -- "$@"`
if [[ $? -ne 0 ]]; then
    exit 2
fi
eval set -- "$PARSED"

while true; do
    case "$1" in
        -d|--directory)
            DIR_DATA="$2"
            shift 2
            ;;
        -n|--name)
            NAME="$2"
            shift 2
            ;;
        -x|--organism)
            ORGANISM="$2"
            shift 2
            ;;
        -r|--resolution)
            RESOLUTION_string="$2"
            RESOLUTION=${RESOLUTION_string/kb/000}
            shift 2
            ;;
        -o|--log)
            FILE_LOG="$2"
            shift 2
            ;;
        -t|--intra)
        	# if only intra chromosome TRUE, otherwise FALSE (default TRUE)
            FLAG_INTRA="$2"
            shift 2
            ;;
        -e|--normalization)
            FLAG_NORM="$2"
            shift 2
            ;;
        -c|--raw)
        	# whether making raw matrices or not
            FLAG_RAW="$2"
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

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${FILE_LOG}" ] && echo "Please specify log file" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
FLAG_INTRA=${FLAG_INTRA:-TRUE}
FLAG_NORM=${FLAG_NORM:-TRUE}
FLAG_RAW=${FLAG_RAW:-TRUE}

TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)

case $ORGANISM in
	pombe) DIR_contig=/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v2.30
		   CHRs=(I II III)
	       LENGTH=(5579133 4539804 2452883) ;;
	human) DIR_contig=/wistar/noma/Data/Human_seq/hg19
		   #CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
	       #LENGTH=(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560 59373566) ;;
		   CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
	       LENGTH=(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560) ;;
	human_EBV) DIR_contig=/wistar/noma/Data/Human_seq/hg19_EBV
		   CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY EBV)
	       LENGTH=(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560 59373566 171823) ;;
	mouse) DIR_contig=/wistar/noma/Data/Mouse_seq/mm10
		   CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
	       LENGTH=(195471971 182113224 160039680 156508116 151834684 149736546 145441459 129401213 124595110 130694993 122082543 120129022 120421639 124902244 104043685 98207768 94987271 90702639 61431566 171031299 91744698) ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac


#==============================================================
# matrix用のフォルダの作成
#==============================================================
if [ ! -e ${DIR_DATA}/${NAME} ]; then
	mkdir ${DIR_DATA}/${NAME}
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string} ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/Raw ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}/Raw
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE
fi
if [ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/InterBin ]; then
	mkdir ${DIR_DATA}/${NAME}/${RESOLUTION_string}/InterBin
fi



#==============================================================
# Raw matrixの作成
#==============================================================
if [ $FLAG_RAW = "TRUE" ]; then
	if [ $FLAG_INTRA = "TRUE" ]; then
		PRO_RAW_matrix=${DIR_LIB}/Make_association_from_fragmentdb_onlyIntraChr.pl
	else
		PRO_RAW_matrix=${DIR_LIB}/Make_association_from_fragmentdb_allChromosome.pl
	fi
	JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "badFragment_${NAME}" | cut -f2 -d' ' | xargs))
	JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
	DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
	sbatch -N 1 -n 5 --exclusive=user --job-name=RAW_${NAME}_${RESOLUTION_string} $DEPEND -o ${FILE_LOG} --open-mode append --wrap="cd ${DIR_DATA}; perl $PRO_RAW_matrix -i ${NAME}_fragment.db -o ${NAME}/${RESOLUTION_string}/Raw/  -r ${RESOLUTION} -b ${NAME}_bad_fragment.txt"


	JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "RAW_${NAME}_${RESOLUTION_string}" | cut -f2 -d' ' | xargs))
	JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
	DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
	if [ $FLAG_INTRA = "TRUE" ]; then
		for i in `seq 1 ${#CHRs[@]}`
		do
			let index=i-1
			CHR=${CHRs[index]}
			sbatch -n 1 --job-name=covRaw_${NAME}_${RESOLUTION_string}_${CHR} $DEPEND -o ${FILE_LOG} --open-mode append --wrap="cd ${DIR_DATA}/${NAME}/${RESOLUTION_string}/Raw; Rscript --slave --vanilla ${DIR_LIB}/../../Conv/Convert_matrix_to_object.R -i ${CHR}.matrix"
		done
	else
		sbatch -n 1 --job-name=covRaw_${NAME}_${RESOLUTION_string} $DEPEND -o ${FILE_LOG} --open-mode append --wrap="cd ${DIR_DATA}/${NAME}/${RESOLUTION_string}/Raw; Rscript --slave --vanilla ${DIR_LIB}/../../Conv/Convert_matrix_to_object.R -i ALL.matrix"
	fi
fi


#==============================================================
# binごとのinter-chromosomeのデータを計算
#==============================================================
if [ $FLAG_INTRA = "TRUE" ]; then
	JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "badfragment_${NAME}" | cut -f2 -d' ' | xargs))
	JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
	DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
	[ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/InterBin/${CHRs[0]}.txt ] && sbatch -n 1 --job-name=InterBin_${NAME}_${RESOLUTION_string} $DEPEND -o ${FILE_LOG} --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Make_association_from_fragmentdb_interChromosome_perBin.pl -i ${NAME}_fragment.db -o ${NAME}/${RESOLUTION_string}/InterBin/  -r ${RESOLUTION} -b ${NAME}_bad_fragment.txt"
fi



#==============================================================
# すべてのraw readの合計値を計算する(intra+inter)
#==============================================================
if [ ! -e ${DIR_DATA}/${NAME}/TotalRead.txt ]; then
	ALREDY_WORKING=$(squeue -o "%j %F" -u htanizawa | grep -e "TotalRaw_${NAME}"  | cut -f2 -d' ' | xargs)

	if [ -n "$ALREDY_WORKING" ]; then
		JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "InterBin_${NAME}_${RESOLUTION_string}" -e "RAW_${NAME}_${RESOLUTION_string}" | cut -f2 -d' ' | xargs))
		JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
		DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
		sbatch -n 1 --job-name=TotalRaw_${NAME} $DEPEND -o  ${FILE_LOG} --export=DIR_DATA=${DIR_DATA},NAME=${NAME},DIR_LIB=${DIR_LIB},RESOLUTION=${RESOLUTION_string},ORGANISM=${ORGANISM} --open-mode append ${DIR_LIB}/Count_total_Raw_read.sh
	fi
fi



#==============================================================
# ICE normalization
#==============================================================
if [ $FLAG_NORM = "TRUE" ]; then
	JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "covRaw_${NAME}_${RESOLUTION_string}" -e "InterBin_${NAME}_${RESOLUTION_string}" | cut -f2 -d' ' | xargs))
	JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
	DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
	if [ $FLAG_INTRA = "TRUE" ]; then
		for i in `seq 1 ${#CHRs[@]}`
		do
			let index=i-1
			CHR=${CHRs[index]}
			[ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION_string}/ICE/${CHR}.rds ] && sbatch -N 1 -n 5 --exclusive=user --job-name=ICE_${NAME}_${RESOLUTION_string}_${CHR} $DEPEND -o ${FILE_LOG} --open-mode append --wrap="cd ${DIR_DATA}/${NAME}/${RESOLUTION_string}; Rscript --vanilla --slave ${DIR_LIB}/Bias_normalization.R -i Raw/${CHR}.matrix -o ICE/${CHR}.matrix --inter InterBin/${CHR}.txt --times 30 && Rscript --slave --vanilla ${DIR_LIB}/../../Conv/Convert_matrix_to_object.R -i ICE/${CHR}.matrix"
		done
	else
		sbatch -N 1 -n 5 --exclusive --job-name=ICE_${NAME}_${RESOLUTION_string} $DEPEND -o ${FILE_LOG} --open-mode append --wrap="cd ${DIR_DATA}/${NAME}/${RESOLUTION_string}; Rscript --vanilla --slave ${DIR_LIB}/Bias_normalization.R -i Raw/ALL.matrix -o ICE/ALL.matrix --times 30 && Rscript --slave --vanilla ${DIR_LIB}/../../Conv/Convert_matrix_to_object.R -i ICE/ALL.matrix"
	fi
fi

