#!/bin/bash

cd ${DIR_DATA}
TOTAL_LINE=`cat ${NAME}_list.txt | wc -l`
if test $TOTAL_LINE -ge $SLURM_ARRAY_TASK_ID
then
	TARGET_FILE=`head -n $SLURM_ARRAY_TASK_ID ${NAME}_list.txt | tail -n 1`
	perl ${DIR_LIB}/Count_mixed_count.pl -i ${TARGET_FILE}
fi
