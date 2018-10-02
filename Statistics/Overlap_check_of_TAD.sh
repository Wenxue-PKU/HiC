#!/bin/bash
# compare TAD domains
# Tags : TAD

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-a, --file1 [domain file1]
		domain file name

	-b, --file2 [domain file2]
		domain file name

	-s, --sufix
		sufix of input file name (default .txt)

	-o, --out [output file]
		output file name

	-d, --distance [distance]
		distance to recognize as same (default 1000)
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hva:b:d:o:s:
LONG=help,version,file1:,file2:,distance:,out:,sufix:
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
		-a|--file1)
			FILE_BED1="$2"
			shift 2
			;;
		-b|--file2)
			FILE_BED2="$2"
			shift 2
			;;
		-s|--sufix)
			SUFIX="$2"
			shift 2
			;;
		-o|--out)
			FILE_OUT="$2"
			shift 2
			;;
		-d|--distance)
			DISTANCE="$2"
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

[ ! -n "${FILE_BED1}" ] && echo "Please specify file1" && exit 1
[ ! -n "${FILE_BED2}" ] && echo "Please specify file2" && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output file name" && exit 1
DISTANCE=${DISTANCE:-1000}
SUFIX=${SUFIX:-.txt}

NAME1=$(basename $FILE_BED1 $SUFIX)
NAME2=$(basename $FILE_BED2 $SUFIX)

#==============================================================
# bed fileを作成
#==============================================================
cat ${FILE_BED1} | awk -v OFS='\t' '{print $1,$2,$2; print $1,$3,$3}' | sort -k1,1 -k2,2n | uniq > ${FILE_OUT}_tmp1.txt
cat ${FILE_BED2} | awk -v OFS='\t' '{print $1,$2,$2; print $1,$3,$3}' | sort -k1,1 -k2,2n | uniq > ${FILE_OUT}_tmp2.txt


#==============================================================
# common border number
#==============================================================
NUM1=$(cat ${FILE_OUT}_tmp1.txt | wc -l)
NUM2=$(cat ${FILE_OUT}_tmp2.txt | wc -l)
COM1=$(bedtools window -w $DISTANCE -u -a ${FILE_OUT}_tmp1.txt -b ${FILE_OUT}_tmp2.txt | wc -l)
COM2=$(bedtools window -w $DISTANCE -u -b ${FILE_OUT}_tmp1.txt -a ${FILE_OUT}_tmp2.txt | wc -l)
P1=$(echo "scale=2; $COM1/$NUM1*100" | bc)
P2=$(echo "scale=2; $COM2/$NUM2*100" | bc)

printf "%s\t%s\n" "Combinations" "$NAME1 vs $NAME2" > ${FILE_OUT}
printf "%s\t%d\n" "Border in $NAME1" $NUM1 >> ${FILE_OUT}
printf "%s\t%d\n" "Common border in $NAME1" $COM1 >> ${FILE_OUT}
printf "%s\t%.2f\n" "% of common in $NAME1" ${P1} >> ${FILE_OUT}
printf "%s\t%d\n" "Border in $NAME2" $NUM2 >> ${FILE_OUT}
printf "%s\t%d\n" "Common border in $NAME2" $COM2 >> ${FILE_OUT}
printf "%s\t%.2f\n" "% of common in $NAME2" ${P2} >> ${FILE_OUT}

rm ${FILE_OUT}_tmp1.txt ${FILE_OUT}_tmp2.txt
