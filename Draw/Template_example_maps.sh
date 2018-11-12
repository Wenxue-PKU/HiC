#!/bin/bash
# Template of Example Hi-C plot

TIME_STAMP=$(date +"%Y-%m-%d")
PROJECT=/wistar/bioinfo-nfs/hideki_projects/
DIR_LIB=${PROJECT}/lib/HiC
DIR_LOG=${PROJECT}/log
DIR_DATA=${PROJECT}/data
DIR_OUT=${PROJECT}/out
DIR_SRC=${PROJECT}/src
DIR_SH=${PROJECT}/sh
DIR_DB=${PROJECT}/db
DIR_DOC=${PROJECT}/doc
DIR_REC=${PROJECT}/received_data
DIR_REP=${PROJECT}/report

FILE_DB=${DIR_DB}/Data.db
CHROM_SIZE=/wistar/noma/Data/Human_seq/hg19/LENGTH.txt


#==============================================================
# 描画する領域をデータベースに登録
#==============================================================
DB_loc=${DIR_OUT}/location.db
file2database.R -i ${DIR_OUT}/location.txt --id TRUE --db ${DB_loc} --table loc


#==============================================================
# Drawing HiC map 
#==============================================================
RESOLUTION=40kb
COL=red
[ ! -e ${DIR_OUT}/img ] && mkdir ${DIR_OUT}/img
for id in $(sqlite3 ${DB_loc} "select id from loc")
do
	CHR=$(sqlite3 ${DB_loc} "select chr from loc where id='${id}'")
	START=$(sqlite3 ${DB_loc} "select start from loc where id='${id}'")
	END=$(sqlite3 ${DB_loc} "select end from loc where id='${id}'")

	for NAME in samples
	do
		sbatch -n 4 --job-name=${id}_${NAME}_hic $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_example_graph_making.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw/Draw_matrix.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds --normalize NA --zero NA --na na --chr ${CHR} --start ${START} --end ${END} --unit p --max 0.95 --color $COL --width 800 -o ${DIR_OUT}/img/${id}_${NAME}_hic.png"

		sbatch -n 4 --job-name=${id}_${NAME}_tad $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_example_graph_making.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw/Draw_borderStrength.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds --chr ${CHR} --start ${START} --end ${END} --width 800 --height 50 --out ${DIR_OUT}/img/${id}_${NAME}_tad.png"

		sbatch -n 4 --job-name=${id}_${NAME}_comp $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_example_graph_making.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw/Draw_PCAscore_from_location.R -i ${DIR_DATA}/${NAME}/Compartment_40kb.txt --chr ${CHR} --start ${START} --end ${END} --fill TRUE --out ${DIR_OUT}/img/${id}_${NAME}_comp.png --ymin \" -60\" --ymax 40 --width 800 --height 150 --line \" -20\""
	done

	#==============================================================
	# axis
	#==============================================================
	sbatch -n 1 --job-name=gh_axis_${id} -o "${DIR_LOG}/${TIME_STAMP}_example_graph_making.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw/Draw_axis.R --chr ${CHR} --start ${START} --end ${END} --width 800 --out ${DIR_OUT}/img/${id}_axis.png"
done


#==============================================================
# Report
#==============================================================
REPORT_TITLE=${TIME_STAMP}_
FILE_md=${DIR_OUT}/${REPORT_TITLE}.md
FILE_html=${DIR_REP}/${REPORT_TITLE}.html

cat <<EOF > $FILE_md
# title of report
<div style="text-align:right">${TIME_STAMP}</div>

## Method

## Result
<style type="text/css">
<!--
.hk_cell{
display : inline-block;
padding: 10px;
vertical-align: top;
font-weight: bold; 
text-align: center;
}
//-->
</style>
EOF

for id in $(seq 1 10)
do
	CHR=$(sqlite3 ${DB_loc} "select chr from loc where id='${id}'")
	START=$(sqlite3 ${DB_loc} "select start from loc where id='${id}'")
	END=$(sqlite3 ${DB_loc} "select end from loc where id='${id}'")

	for NAME in samples
	do
		cat <<-EOF >> ${FILE_md}
		<div class="hk_cell">${NAME} ${CHR}:${START}-${END}</br>
		<img width=400 src="img/${id}_${NAME}_hic.png"/></br>
		<img width=400 src="img/${id}_${NAME}_comp.png"/></br>
		<img width=400 src="img/${id}_${NAME}_tad.png"/></br>
		<img width=400 src="img/${id}_axis.png"/>
		</div>

		EOF
	done
done

cd ${DIR_OUT} && pandoc $FILE_md -s --self-contained -t html5 -c ~/.pandoc/github.css  -o $FILE_html


