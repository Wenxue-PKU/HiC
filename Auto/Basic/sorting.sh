#!/bin/bash

cd ${DIR_DATA}
TOTAL_LINE=`cat ${NAME}_list.txt | wc -l`
if test $TOTAL_LINE -ge $SLURM_ARRAY_TASK_ID
then
	TARGET=`head -n $SLURM_ARRAY_TASK_ID ${NAME}_list.txt | tail -n 1`
	sort -k 3,3n -k 4,4 $TARGET > ${TARGET}_sorted
	mv ${TARGET}_sorted $TARGET
fi

