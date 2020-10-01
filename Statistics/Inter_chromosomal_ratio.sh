#!/bin/bash
# inter-chromosomeal percentage

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	r, --resolution
		bin size for analysis (ex. 10kb)

	-i, --in [map file]
		map file

	-o, --out [output file]
		output file

EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvr:i:o:
LONG=help,version,resolution:,in:,out:
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
		-r|--resolution)
			RESOLUTION_string="$2"
			shift 2
			;;
		-i|--in)
			FILE_map="$2"
			shift 2
			;;
		-o|--out)
			FILE_out="$2"
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
TIME_STAMP=$(date +"%Y-%m-%d_%H.%M.%S")
INPUT_FILES=$@

[ ! -n "${FILE_map}" ] && echo "Please specify input map file" && exit 1
[ ! -n "${FILE_out}" ] && echo "Please specify output file" && exit 1
[ ! -n "${RESOLUTION_string}" ] && echo "Please specify resolution" && exit 1
[ ! -e ${DIR_out} ] && mkdir ${DIR_out}

RESOLUTION=${RESOLUTION_string/kb/000}


#==============================================================
# count score
#==============================================================
echo "chr start target_chr count" | tr ' ' '\t' > $FILE_out

zcat $FILE_map | tail -n +2 | awk -v OFS='\t' -v BIN=$RESOLUTION '{
	bin1=int($3/BIN)*BIN;
	bin2=int($10/BIN)*BIN;
	INTER[$2"\t"bin1"\t"$9]+=1;
	INTER[$9"\t"bin2"\t"$2]+=1;
}END{
	for(x in INTER){
		print x,INTER[x]
	}
}' >> ${FILE_out}








