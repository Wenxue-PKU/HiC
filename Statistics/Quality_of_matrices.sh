#!/bin/bash
# Quality control program

#==============================================================
# コマンド引数の受取り
#==============================================================
SHORT=d:o:x:
LONG=directory:,output:,organism:
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
        -o|--output)
            FILE_OUT="$2"
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

[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1


TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)

case $ORGANISM in
	pombe) DIR_contig=/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v2.30
		   CHRs=(I II III)
	       LENGTH=(5579133 4539804 2452883) ;;
	human) DIR_contig=/wistar/noma/Data/Human_seq/hg19
		   CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
	       LENGTH=(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560 59373566) ;;
	human_EBV) DIR_contig=/wistar/noma/Data/Human_seq/hg19_EBV
		   CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY EBV)
	       LENGTH=(249250621 243199373 198022430 191154276 180915260 171115067 159138663 146364022 141213431 135534747 135006516 133851895 115169878 107349540 102531392 90354753 81195210 78077248 59128983 63025520 48129895 51304566 155270560 59373566 171823) ;;
	mouse) DIR_contig=/wistar/noma/Data/Mouse_seq/mm10
		   CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
	       LENGTH=(195471971 182113224 160039680 156508116 151834684 149736546 145441459 129401213 124595110 130694993 122082543 120129022 120421639 124902244 104043685 98207768 94987271 90702639 61431566 171031299 91744698) ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac

FILE_tmp_CHR=${FILE_OUT}.tmpCHR
FILE_tmp_NA=${FILE_OUT}.tmpNA
FILE_tmp_SD=${FILE_OUT}.tmpSD
FILE_tmp_Bias=${FILE_OUT}.tmpBias

for i in `seq 1 ${#CHRs[@]}`
do
	let index=i-1
	CHR=${CHRs[index]}
	echo "${CHR}" >> $FILE_tmp_CHR
	Rscript --vanilla --slave ${DIR_LIB}/ICE_success_check.R -i ${DIR_DATA}/${CHR}.rds --method p >> ${FILE_tmp_NA}
	Rscript --vanilla --slave ${DIR_LIB}/Map/Map_property.R -i ${DIR_DATA}/${CHR}.rds --lineMAX TRUE --info SDAve >> ${FILE_tmp_Bias}
done

echo -e "Chromosome\tNA%\tBias\" >  ${FILE_OUT}
paste $FILE_tmp_CHR ${FILE_tmp_NA} ${FILE_tmp_Bias} >> ${FILE_OUT}
rm $FILE_tmp_CHR ${FILE_tmp_NA} ${FILE_tmp_Bias}
