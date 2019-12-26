#!/bin/bash
# Draw example images with circles

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

	-o, --out [output directory]
		output directory
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvd:i:o:
LONG=help,version,data:,in:,out:
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
			DIR_OUT="$2"
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
[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1

DIR_tmp=$(mktemp -d /tmp/tmp_sig.XXXXX)

### circle
awk -v OFS='\t' 'NR>1{printf "%s:%d:%d\t%s:%d:%d\n",$1,$2,$3,$4,$5,$6}' $FILE_IN > ${DIR_tmp}/circle.txt

### Upstream 4Mb, Downstream 4Mbを描画する(上位10位の場所)
CHECK_WIDTH=4000000
awk -v OFS='\t' -v W=$CHECK_WIDTH '
	NR==1{print "name",$1,$2,$3}
	(NR>1&&NR<12){m=($2+$5)/2; if(m-W < 1){m=W+1}; print NR-1,$1,m-W,m+W}
	' $FILE_IN > ${DIR_tmp}/location.txt

### Draw graph
Rscript --vanilla --slave ${DIR_LIB}/../../Draw/Draw_matrix_batch.R --in ${DIR_DATA} --out ${DIR_OUT} --location ${DIR_tmp}/location.txt --color red --width 1000 --circle ${DIR_tmp}/circle.txt


rm -rf $DIR_tmp