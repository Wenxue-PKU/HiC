#!/bin/sh
# Raw matrix のすべてのreadの数の合計値を出す

### 以下の変数を指定する必要あり
# DIR_DATA, NAME, DIR_LIB, RESOLUTION, ORGANISM


FILE_TMP=${DIR_DATA}/tmp_${NAME}_totalReadCount.txt

[ ! -n "${RESOLUTION}" ] && echo "Please specify RESOLUTION" && exit 1

case $ORGANISM in
	pombe) CHRs=(I II III) ;;
	human) CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY) ;;
	human_EBV) CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY EBV) ;;
	mouse) CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY) ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac


for i in `seq 1 ${#CHRs[@]}`
do
	let index=i-1
	CHR=${CHRs[index]}

	# intra-chromosome (upper triangularの数の合計。対角線も入れる)
	if [ -e ${DIR_DATA}/${NAME}/${RESOLUTION}/Raw/${CHR}.rds ]; then
		Rscript --vanilla --slave ${DIR_LIB}/../../Statistics/Map/Map_property.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/Raw/${CHR}.rds --info SUM >> $FILE_TMP
	else
		echo "${NAME}/${RESOLUTION}/Raw/${CHR}.rds is not exist. Skip it."
	fi

	# inter-chromsoome (同じ組み合わせを２回数えることになるので、半分にする)
	if [ -e ${DIR_DATA}/${NAME}/${RESOLUTION}/InterBin/${CHR}.txt ]; then
		cat ${DIR_DATA}/${NAME}/${RESOLUTION}/InterBin/${CHR}.txt | awk 'BEGIN{S=0}{S+=$2}END{print S/2}' >> $FILE_TMP
	else
		echo "${NAME}/${RESOLUTION}/InterBin/${CHR}.txt is not exist. Skip it."
	fi

done

# すべてを足す
cat $FILE_TMP | awk 'BEGIN{S=0}{S+=$1}END{print S}' > ${DIR_DATA}/${NAME}/TotalRead.txt

rm $FILE_TMP

