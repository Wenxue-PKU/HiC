#!/bin/bash
# Creating TAD

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-d, --directory [data directory]
		directory name of analysis file locate

	-n, --name [sample name]
		sample name

	-x, --ref [ex. hg19]
		organism name

	--chr_ALL [TRUE/FALSE]
		use ALL.rds for data

	--include [including chromsome list]
		if only specific chromosome should be calculated, specify the list. Separated by ,.

	--exclude [exclude chromsome list]
		list of excluding chromosomes. Separated by ,. Ex. chrM,chrY

	-r, --resolution [40kb]
		default : 40kb

	-f, --force [FALSE(default)|TRUE]
		force to perform analysis even matrices is in-complete or not (default)

	-e, --method [boderStrength|DI]
		default : borderStrength


EOF

}


SHORT=hd:n:x:r:f:e:
LONG=help,directory:,name:,ref:,chr_ALL:,resolution:,force:,method:,include:,exclude:
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
		-d|--directory)
			DIR_DATA="$2"
			shift 2
			;;
		-n|--name)
			NAME="$2"
			shift 2
			;;
		-x|--ref)
			REF="$2"
			shift 2
			;;
		--chr_ALL)
			FLAG_chr_ALL="$2"
			shift 2
			;;
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
			;;
		-f|--force)
			FLAG_FORCE="$2"
			shift 2
			;;
		-e|--method)
			METHOD="$2"
			shift 2
			;;
		--include)
			CHR_include="$2"
			shift 2
			;;
		--exclude)
			CHR_exclude="$2"
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

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${REF}" ] && echo "Please specify ref" && exit 1
RESOLUTION=${RESOLUTION:-40kb}
FLAG_FORCE=${FLAG_FORCE:-FALSE}
METHOD=${METHOD:-borderStrength}
CHR_include=${CHR_include:-NA}
CHR_exclude=${CHR_exclude:-NA}
FLAG_chr_ALL=${FLAG_chr_ALL:-FALSE}

TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)

#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/../../utils/load_setting.sh -x $REF -r NA

#-----------------------------------------------
# Load chromosome length
#-----------------------------------------------
CHR_TABLE=$(Rscript --vanilla --slave ${DIR_LIB}/../../utils/Chromosome_length.R --in $FILE_CHROME_LENGTH --include $CHR_include --exclude $CHR_exclude)
CHRs=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==1' | tr ',' ' '))
LENGTH=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==2' | tr ',' ' '))
CHRs_list=$(echo ${CHRs[@]} | tr ' ' ',')

if [ "${FLAG_chr_ALL}" == "FALSE" ];then
	MISSING=()
	for i in $(seq 1 ${#CHRs[@]})
	do
		let index=i-1
		CHR=${CHRs[index]}
		[ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds ] && MISSING=("${MISSING[@]}" $CHR)
	done
	if [ ${#MISSING[@]} -gt 0 ]; then
		echo "${NAME} (${RESOLUTION}) matrices of ${MISSING[@]} are not exists"
		[ "$FLAG_FORCE" != "TRUE" ] && echo "quit job" && exit 1
	fi
fi

DIR_TAD=${DIR_DATA}/${NAME}/tmp_TAD_${RESOLUTION}
mkdir $DIR_TAD
for i in $(seq 1 ${#CHRs[@]})
do
	let index=i-1
	CHR=${CHRs[index]}
	if [ "${FLAG_chr_ALL}" == "TRUE" ]; then
		FILE_input=${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/ALL.rds
	else
		FILE_input=${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds
	fi
	if [ "$METHOD" = "borderStrength" ]; then
		Rscript --vanilla --slave ${DIR_LIB}/../../Draw/Draw_borderStrength.R -i ${FILE_input} --location ${DIR_TAD}/${CHR}.txt --chr ${CHR}
	else
		[ ! -e ${DIR_DATA}/${NAME}/${RESOLUTION}/HiTC ] && mkdir ${DIR_DATA}/${NAME}/${RESOLUTION}/HiTC
		Rscript --vanilla --slave ${DIR_LIB}/../../Draw/Draw_Dixon_directionaryIndex.R -i ${FILE_input} -t ${DIR_DATA}/${NAME}/${RESOLUTION}/HiTC/${CHR}.rds --location ${DIR_TAD}/${CHR}.txt --chr ${CHR} -r ${RESOLUTION}
	fi
done

if [ "$METHOD" = "borderStrength" ]; then
	cat ${DIR_TAD}/*.txt | sort -k1,1 -k2,2n > ${DIR_DATA}/${NAME}/TAD_${RESOLUTION}_score.txt
	perl ${DIR_LIB}/Convert_border2bed.pl -i ${DIR_DATA}/${NAME}/TAD_${RESOLUTION}_score.txt -o ${DIR_DATA}/${NAME}/TAD_${RESOLUTION}.txt -c 6 -t 1
else
	cat ${DIR_TAD}/*.txt | sort -k1,1 -k2,2n > ${DIR_DATA}/${NAME}/DI_${RESOLUTION}_score.txt
	perl ${DIR_LIB}/Convert_border2bed.pl -i ${DIR_DATA}/${NAME}/DI_${RESOLUTION}_score.txt -o ${DIR_DATA}/${NAME}/DI_${RESOLUTION}.txt -c 5
fi
rm -r ${DIR_TAD}






