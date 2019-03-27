#!/bin/bash
# Output example Hi-C contact maps

get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [target sample name(s). separated by space]

Format of location file should have following columns, 
1. resolution
2. chr
3. start
4. end
optional field is 
5. moving_average
6. name (use for region name)
7. linev and lineh (line location for map)
File should be deliminated by tab. Order is flexible and okay to have extrac columns but should have header as described in above.

Description
	-h, --help
		show help

	-i, --in [location file]
		file defined location of drawing. File format is as described in above
	
	-o, --out [output directory]
		output directory

	-d, --data [data directory]
		data directory

	-c, --color [color list]
		specify color for each samples. separated by space but surrounded with double quatation
	
	--circles [circle location file]
		file for circle drawing
EOF

}

SHORT=hvi:o:d:c:
LONG=help,version,in:,out:,data:,color:,circles:
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
		-i|--in)
			FILE_location="$2"
			shift 2
			;;
		-o|--out)
			DIR_OUT="$2"
			shift 2
			;;
		-d|--data)
			DIR_DATA="$2"
			shift 2
			;;
		-c|--color)
			COLORS=($2)
			shift 2
			;;
		--circles)
			FILE_circles="$2"
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

[ ! -n "${FILE_location}" ] && echo "Please specify location file" && exit 1
[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${COLORS}" ] && echo "Please specify colors" && exit 1
[ $# -lt 1 ] && echo "Please specify target Hi-C sample name(s)" && exit 1
[ ! -e "${FILE_location}" ] && echo "Location file is not exists" && exit 1
FILE_circles=${FILE_circles:-"NULL"}

SAMPLES=$@


### Check optional field
FLAG_moving_average=$(cat ${FILE_location} | head -n1 | grep -c moving_average)
FLAG_name=$(cat ${FILE_location} | head -n1 | grep -c name)
FLAG_lineh=$(cat ${FILE_location} | head -n1 | grep -c lineh)
FLAG_linev=$(cat ${FILE_location} | head -n1 | grep -c linev)

#==============================================================
# 描画する領域をデータベースに登録
#==============================================================
DB_loc=${FILE_location/.txt/.db}
[ ! -e ${DB_loc} ] && file2database.R -i ${FILE_location} --id TRUE --db ${DB_loc} --table loc


#==============================================================
# Drawing HiC map 
#==============================================================
[ ! -e ${DIR_OUT}/img ] && mkdir ${DIR_OUT}/img
[ ! -e ${DIR_OUT}/log ] && mkdir ${DIR_OUT}/log
for id in $(sqlite3 ${DB_loc} "select id from loc")
do
	RESOLUTION=$(sqlite3 ${DB_loc} "select resolution from loc where id='${id}'")
	CHR=$(sqlite3 ${DB_loc} "select chr from loc where id='${id}'")
	START=$(sqlite3 ${DB_loc} "select start from loc where id='${id}'")
	END=$(sqlite3 ${DB_loc} "select end from loc where id='${id}'")

	if [ $FLAG_moving_average -eq 0 ]; then
		MOVING_AVERAGE=0
	else
		MOVING_AVERAGE=$(sqlite3 ${DB_loc} "select moving_average from loc where id='${id}'")
	fi

	if [ $FLAG_lineh -eq 0 ]; then
		DRAW_LINE_H=""
	else
		DRAW_LINE_H="--lineh_chr ${CHR} --lineh_pos $(sqlite3 ${DB_loc} "select lineh from loc where id='${id}'")"
	fi

	if [ $FLAG_linev -eq 0 ]; then
		DRAW_LINE_V=""
	else
		DRAW_LINE_V="--linev_chr ${CHR} --linev_pos $(sqlite3 ${DB_loc} "select linev from loc where id='${id}'")"
	fi	


	i=0
	for NAME in $SAMPLES
	do
		COL=${COLORS[$i]}
		let i=${i}+1
		sbatch -n 4 --job-name=${id}_${NAME}_hic $(sq --node) -o "${DIR_OUT}/log/${TIME_STAMP}_map_for_${id}_${NAME}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw_matrix.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds --normalize NA --zero NA --na na --moving_average ${MOVING_AVERAGE} --chr ${CHR} --start ${START} --end ${END} --unit p --max 0.95 --color $COL --width 800 -o ${DIR_OUT}/img/${id}_${NAME}_hic1.png --circle $FILE_circles $DRAW_LINE_H $DRAW_LINE_V && convert -rotate -45 -crop 1134x300-167+100 -resize 800x ${DIR_OUT}/img/${id}_${NAME}_hic1.png ${DIR_OUT}/img/${id}_${NAME}_hic2.png"

		# sbatch -n 4 --job-name=${id}_${NAME}_tad $(sq --node) -o "${DIR_OUT}/log/${TIME_STAMP}_tad_for_${id}_${NAME}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw_borderStrength.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds --chr ${CHR} --start ${START} --end ${END} --width 800 --height 50 --out ${DIR_OUT}/img/${id}_${NAME}_tad.png"

		# sbatch -n 4 --job-name=${id}_${NAME}_comp $(sq --node) -o "${DIR_OUT}/log/${TIME_STAMP}_comp_for_${id}_${NAME}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw_PCAscore_from_location.R -i ${DIR_DATA}/${NAME}/Compartment_40kb.txt --chr ${CHR} --start ${START} --end ${END} --fill TRUE --out ${DIR_OUT}/img/${id}_${NAME}_comp.png --ymin \" -60\" --ymax 40 --width 800 --height 150 --line \" -20\""
	done

	#==============================================================
	# axis
	#==============================================================
	sbatch -n 1 --job-name=gh_${id}_axis $(sq --node) -o "${DIR_OUT}/log/${TIME_STAMP}_axis_for_${id}_${CHR}_${START}_${END}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Draw_axis.R --chr ${CHR} --start ${START} --end ${END} --width 800 --out ${DIR_OUT}/img/${id}_axis.png"
done
