#!/bin/bash
# drawing individual hic map

id=$SLURM_ARRAY_TASK_ID

CHR1=$(sqlite3 ${DB_loc} "select chr1 from loc where id='${id}'")
START1=$(sqlite3 ${DB_loc} "select start1 from loc where id='${id}'")
END1=$(sqlite3 ${DB_loc} "select end1 from loc where id='${id}'")

REGION=$(sqlite3 ${DB_loc} "select name from loc where id='${id}'" 2> /dev/null)
CHR2=$(sqlite3 ${DB_loc} "select chr2 from loc where id='${id}'" 2> /dev/null)
START2=$(sqlite3 ${DB_loc} "select start2 from loc where id='${id}'" 2> /dev/null)
END2=$(sqlite3 ${DB_loc} "select end2 from loc where id='${id}'" 2> /dev/null)
MOVING_AVERAGE=$(sqlite3 ${DB_loc} "select moving_average from loc where id='${id}'" 2> /dev/null)

CHR2=${CHR2:-"${CHR1}"}
START2=${START2:-"${START1}"}
END2=${END2:-"${END1}"}
MOVING_AVERAGE=${MOVING_AVERAGE:-"0"}
REGION=${REGION:-"${id}"}


if [ $RESOLUTION = "NA" ]; then
	RESOLUTION=$(sqlite3 ${DB_loc} "select resolution from loc where id='${id}'")
fi

if [ $FLAG_lineh -eq 0 ]; then
	DRAW_LINE_H=""
else
	DRAW_LINE_H="--lineh_chr ${CHR1} --lineh_pos $(sqlite3 ${DB_loc} "select lineh from loc where id='${id}'")"
fi

if [ $FLAG_linev -eq 0 ]; then
	DRAW_LINE_V=""
else
	DRAW_LINE_V="--linev_chr ${CHR2} --linev_pos $(sqlite3 ${DB_loc} "select linev from loc where id='${id}'")"
fi	

[ ! -e ${DIR_OUT}/img/${REGION}_${NAME}.png ] && Rscript --vanilla --slave ${DIR_LIB}/Draw_matrix.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/${MAP_TYPE}/${CHR1}.rds --normalize NA --zero NA --na na --moving_average ${MOVING_AVERAGE} --chr ${CHR1} --start ${START1} --end ${END1} --chr2 ${CHR2} --start2 ${START2} --end2 ${END2} --unit p --max 0.95 --color $COL --width 500 -o ${DIR_OUT}/img/${REGION}_${NAME}.png --circle $FILE_circles $DRAW_LINE_H $DRAW_LINE_V

