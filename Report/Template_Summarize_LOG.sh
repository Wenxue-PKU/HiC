#!/bin/bash

DIR_DATA=/wistar/bioinfo-nfs/hideki_projects/325_20170727_HiC_Rugang_HUman_WuShui/log
FILE_OUT=/wistar/bioinfo-nfs/hideki_projects/325_20170727_HiC_Rugang_HUman_WuShui/out/2017-08-01_read_summary.txt

cd ${DIR_DATA}

tmpfile=$(mktemp tmp.XXXXXX)
tmpfile2=$(mktemp tmp.XXXXXX)

FLAG=0
for FILE in `ls -1 ????-??-??_count_*.log`
do
	TMP=${FILE#????-??-??_count_}
	NAME=${TMP%.log}

	FILE_map=${FILE/count/make_map}

	if [ $FLAG -eq 0 ]
	then
		echo "SAMPLE name" > $FILE_OUT
		tail -n 4 $FILE_map | cut -f1 -d: >> $FILE_OUT
		cut -f1 -d: ${FILE} >> $FILE_OUT
		FLAG=1
	fi

	echo $NAME > $tmpfile
	tail -n 4 $FILE_map | cut -f2 >> $tmpfile
	cut -f2 -d: ${FILE} | sed 's/^ //' >> $tmpfile
	paste $FILE_OUT $tmpfile > $tmpfile2
	mv $tmpfile2 $FILE_OUT
done

rm $tmpfile
