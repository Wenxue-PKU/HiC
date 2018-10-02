#!/bin/bash
# Sequenceの終わったファイルをまとめる

TIME_STAMP=$(date +"%Y-%m-%d")
PROJECT=/wistar/bioinfo-nfs/hideki_projects/2017-08-08_HiC_human_senescent
DIR_DATA=${PROJECT}/data/Sequence
DIR_LOG=${PROJECT}/log


# 名前をまとめたファイル  1: userのつけた元々の名前 2: sampleにつける名前 3: fastqの名前
FILE_sampleMap=${DIR_DATA}/sampleMap.txt


#### 下のフラグをONにするとjobを実行。OFFだと見るだけ
FLAG_RUN=OFF


for DATA in $(cat $FILE_sampleMap | xargs -n 3 | sed -e 's/ /:/g')
do
	SAMPLE_NAME=$(echo $DATA | cut -f1 -d':')
	NEW_NAME=$(echo $DATA | cut -f2 -d':')
	FASTQ_NAME=$(echo $DATA | cut -f3 -d':')
	#echo "1:$SAMPLE_NAME, 2:$NEW_NAME, 3:$FASTQ_NAME"

	FLAG_PAIRED_END=FALSE
	STOCK1=()
	STOCK2=()
    for r in $(seq 1 4)
    do
		FILE_ORI1=${FASTQ_NAME}${r}_R1_001.fastq.gz
		FILE_ORI2=${FASTQ_NAME}${r}_R2_001.fastq.gz
		[ -e ${DIR_DATA}/$FILE_ORI1 ] && STOCK1=("${STOCK1[@]}" $FILE_ORI1)
		[ -e ${DIR_DATA}/$FILE_ORI2 ] && STOCK2=("${STOCK2[@]}" $FILE_ORI2) && FLAG_PAIRED_END=TRUE
    done

    [ ${#STOCK1[@]} -eq 0 ] && echo "no original sequence found. Stop script" && exit

    if [ $FLAG_PAIRED_END == "TRUE" ]; then
		FILE_NEW1=${NEW_NAME}_1.fastq
		FILE_NEW2=${NEW_NAME}_2.fastq
		echo "${FILE_NEW1} will be created by "
		echo "${STOCK1[@]}" | xargs -n1 | xargs -n1 -I@ sh -c "echo @"
		[ $FLAG_RUN == "ON" ] && sbatch -n 1 --job-name=zcat_${NEW_NAME}_1 -o "${DIR_LOG}/${TIME_STAMP}_uncompress_${NEW_NAME}_1.log" --open-mode append <<-EOF
		#!/bin/sh
		cd ${DIR_DATA}
		echo "${STOCK1[@]}" | xargs -n1 | xargs -n1 -I@ sh -c "zcat @ >> $FILE_NEW1"
		EOF

		echo "-------------------------------------------------"
		echo "${FILE_NEW2} will be created by "
		echo "${STOCK2[@]}" | xargs -n1 | xargs -n1 -I@ sh -c "echo @"
		echo "-------------------------------------------------"
		[ $FLAG_RUN == "ON" ] && sbatch -n 1 --job-name=zcat_${NEW_NAME}_2 -o "${DIR_LOG}/${TIME_STAMP}_uncompress_${NEW_NAME}_2.log" --open-mode append <<-EOF
		#!/bin/sh
		cd ${DIR_DATA}
		echo "${STOCK2[@]}" | xargs -n1 | xargs -n1 -I@ sh -c "zcat @ >> $FILE_NEW2"
		EOF

    else
		FILE_NEW1=${DIR_DATA}/${NEW_NAME}.fastq
		echo "${FILE_NEW1} will be created by "
		echo "${STOCK1[@]}" | xargs -n1 | xargs -n1 -I@ sh -c "echo @"
		echo "${STOCK1[@]}" | xargs -n1 | xargs -n1 -I@ sh -c "echo @"
		[ $FLAG_RUN == "ON" ] && sbatch -n 1 --job-name=zcat_${NEW_NAME} -o "${DIR_LOG}/${TIME_STAMP}_uncompress_${NEW_NAME}.log" --open-mode append <<-EOF
		#!/bin/sh
		cd ${DIR_DATA}
		echo "${STOCK1[@]}" | xargs -n1 | xargs -n1 -I@ sh -c "zcat @ >> $FILE_NEW1"
		EOF

		echo "-------------------------------------------------"
    fi
done
