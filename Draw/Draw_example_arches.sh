#!/bin/bash
# Output example arches

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Format of location file should have following columns, 
1. chr
2. start
3. end

Description
	-h, --help
		show help

	-n, --name [name]
		name of sample
	
	-i, --in [arch drawing information]
		arch drawing information. At least chr1, start1, end1, chr2, start2, end2 are required
	
	-o, --out [output directory]
		output directory

	-l, --location [location file]
		file defined location of drawing. File format is as described in above
	
	-m, --min [minimum distance to draw]
		minimum distance for drawing. default 50000 (50kb)

	-c, --color [color]
		color for drawing. default grey30
EOF

}

SHORT=hn:i:o:l:m:c:
LONG=help,name:,in:,out:,location:,min:,color:
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
		-n|--name)
			NAME="$2"
			shift 2
			;;
		-i|--in)
			FILE_arch="$2"
			shift 2
			;;
		-o|--out)
			DIR_OUT="$2"
			shift 2
			;;
		-l|--location)
			FILE_location="$2"
			shift 2
			;;
		-m|--min)
			MIN="$2"
			shift 2
			;;
		-c|--color)
			COLOR="$2"
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

[ ! -n "${NAME}" ] && echo "Please specify name" && exit 1
[ ! -n "${FILE_location}" ] && echo "Please specify location file" && exit 1
[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1
[ ! -e "${FILE_location}" ] && echo "Location file is not exists" && exit 1
MIN=${MIN:-50000}
COLOR=${COLOR:-grey30}

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
	CHR=$(sqlite3 ${DB_loc} "select chr from loc where id='${id}'")
	START=$(sqlite3 ${DB_loc} "select start from loc where id='${id}'")
	END=$(sqlite3 ${DB_loc} "select end from loc where id='${id}'")

	xvfb-run Rscript --vanilla --slave ${DIR_LIB}/Draw_arch.R -i ${FILE_arch} --start $START --end $END --min ${MIN} --color ${COLOR} -o ${DIR_OUT}/img/${id}_${NAME}_arch.png

	# sbatch -n 4 --job-name=${id}_${NAME} $(sq --node) -o "${DIR_OUT}/log/${TIME_STAMP}_arch_graph_for_${id}_${NAME}.log" --open-mode append --wrap="Rscript2 --vanilla --slave ${DIR_LIB}/Draw_arch.R -i ${FILE_arch} --start $START --end $END --min ${MIN} --color ${COLOR} -o ${DIR_OUT}/img/${id}_${NAME}_arch.png"
done
