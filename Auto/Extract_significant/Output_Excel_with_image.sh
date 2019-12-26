#!/bin/bash
# Make Excel file with image for significant output


get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-d, --data [HiC map directory]
		data directory of HiC map

	-i, --in [significant file]
		Output file from significant calculation

	-o, --out [output Excel file]
		output Excel file name
	
	-t, --title [title of excel tab]
		excel tab name 
	
	--max [maximum number of output]
		maximum number of output (default: 5000) 0 for no limit
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvd:i:o:t:
LONG=help,version,data:,in:,out:,title:,max:
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
		-d|--data)
			DIR_DATA="$2"
			shift 2
			;;
		-i|--in)
			FILE_IN="$2"
			shift 2
			;;
		-o|--out)
			FILE_OUT="$2"
			shift 2
			;;
		-t|--title)
			EXCEL_TITLE="$2"
			shift 2
			;;
		--max)
			OUT_MAX="$2"
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
INPUT_FILES=$@

[ ! -n "${DIR_DATA}" ] && echo "Please specify hic map directory" && exit 1
[ ! -n "${FILE_IN}" ] && echo "Please specify input significant file" && exit 1
[ ! -e "${FILE_IN}" ] && echo "There is no such file" && exit 1
EXCEL_TITLE=${EXCEL_TITLE:-"Significant"}
OUT_MAX=${OUT_MAX:-5000}

if hash module 2>/dev/null; then
	module load racs-eb 2> /dev/null
	module load pandas/0.21.0-intel-2017a-Python-2.7.13 2> /dev/null
fi

DIR_tmp=$(mktemp -d /tmp/tmp_sig.XXXXX)
FILE_sig=${DIR_tmp}/target.txt

### limit output
if [ $OUT_MAX -ne 0 ]; then
	let n=$OUT_MAX+1
	head -n $n ${FILE_IN} > ${FILE_sig}
else
	cp ${FILE_IN} ${FILE_sig}
fi


### location
echo "name chr1 start1 end1 chr2 start2 end2" | tr ' ' '\t' > ${DIR_tmp}/location.txt
awk -v OFS='\t' 'NR>1{s1=$2-200000; e1=$2+200000; s2=$5-200000; e2=$5+200000; name=NR-1; print name,$1,s1,e1,$4,s2,e2}' ${FILE_sig} >> ${DIR_tmp}/location.txt


### Draw graph
Rscript --vanilla --slave ${DIR_LIB}/../../Draw/Draw_matrix_batch.R --in ${DIR_DATA} --out ${DIR_tmp} --location ${DIR_tmp}/location.txt --color red --width 100 


### Output Excel
python ${DIR_LIB}/Summarize_significant_pairs.py -i ${FILE_sig} -o ${FILE_OUT} --title ${EXCEL_TITLE} --image ${DIR_tmp}

rm -r ${DIR_tmp}

