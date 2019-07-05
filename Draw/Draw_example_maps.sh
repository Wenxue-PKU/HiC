#!/bin/bash
# Output example Hi-C contact maps

get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [target sample name(s). separated by space]

Format of location file should have following columns, 
1. chr1
2. start1
3. end1
optional field are
4. chr2
5. start2
6. end2
7. resolution (if different from global resolution)
8. moving_average
9. name (region name)
10. linev (line location for map)
11. lineh (line location for map)
File should be deliminated by tab. Order is flexible and okay to have extrac columns but should have header as described in above.

Description
	-h, --help
		show help

	-i, --in [location file]
		file defined location of drawing. File format is as described in above

	-r, --resolution [resolution]
		resolution of Hi-C map. (ex. 10kb). If each map need different resolution, write in location file.

	-t, --type [map type ex (ICE2)]
		data map type Raw/ICE/ICE2. default (ICE)

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

SHORT=hvi:r:t:o:d:c:
LONG=help,version,in:,resolution:,type:,out:,data:,color:,circles:,
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
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
			;;
		-t|--type)
			MAP_TYPE="$2"
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
RESOLUTION=${RESOLUTION:-"NA"}
MAP_TYPE=${MAP_TYPE:-"ICE"}

SAMPLES=$@


### Check optional field
FLAG_lineh=$(cat ${FILE_location} | head -n1 | grep -c lineh)
FLAG_linev=$(cat ${FILE_location} | head -n1 | grep -c linev)

#==============================================================
# 描画する領域をデータベースに登録
#==============================================================
DB_loc=${FILE_location/.txt/.db}
[ ! -e ${DB_loc} ] && file2database.R -i ${FILE_location} --id TRUE --db ${DB_loc} --table loc
NUM_data=$(cat $FILE_location | wc -l)

#==============================================================
# Drawing HiC map 
#==============================================================
[ ! -e ${DIR_OUT}/img ] && mkdir ${DIR_OUT}/img

i=0
for NAME in $SAMPLES
do
	COL=${COLORS[$i]}
	let i=${i}+1
	sbatch -N 2 -n 4 --array=1-${NUM_data} --job-name=dr_${NAME} $(sq --node) -o "${DIR_OUT}/drawing_map_for_${NAME}.log" --export=NAME="${NAME}",DIR_DATA="${DIR_DATA}",DIR_OUT="${DIR_OUT}",DB_loc="${DB_loc}",FLAG_lineh="${FLAG_lineh}",FLAG_linev="${FLAG_linev}",RESOLUTION="${RESOLUTION}",COL="${COL}",MAP_TYPE="${MAP_TYPE}",DIR_LIB="${DIR_LIB}",FILE_circles="${FILE_circles}" --open-mode append ${DIR_LIB}/lib/Drawing_individual_map.sh
done

