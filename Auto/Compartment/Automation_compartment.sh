#!/bin/bash
# Creating Compartment score

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

	-a, --ref [hg19|mm10]
		assemble name

	--include [including chromsome list]
		if only specific chromosome should be calculated, specify the list. Separated by ,.

	--exclude [exclude chromsome list]
		list of excluding chromosomes. Separated by ,. Ex. chrM,chrY

	-r, --resolution [200kb|500kb]
		default : 200kb

	-f, --force [FALSE(default)|TRUE]
		force to perform analysis even matrices is in-complete or not (default)

EOF

}


SHORT=hd:n:a:r:f:
LONG=help,directory:,name:,ref:,include:,exclude:,resolution:,force:
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
        -a|--ref)
            ASSEMBLE="$2"
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
        -r|--resolution)
            RESOLUTION="$2"
            shift 2
            ;;
        -f|--force)
            FLAG_FORCE="$2"
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
[ ! -n "${ASSEMBLE}" ] && echo "Please specify assemble version" && exit 1
RESOLUTION=${RESOLUTION:-200kb}
FLAG_FORCE=${FLAG_FORCE:-FALSE}
CHR_include=${CHR_include:-NA}
CHR_exclude=${CHR_exclude:-NA}

TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)

#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/../../utils/load_setting.sh -x $ASSEMBLE -r NA


#-----------------------------------------------
# Load chromosome length
#-----------------------------------------------
CHR_TABLE=$(Rscript --vanilla --slave ${DIR_LIB}/../../utils/Chromosome_length.R --in $FILE_CHROME_LENGTH --include $CHR_include --exclude $CHR_exclude)
CHRs=($(echo $CHR_TABLE | xargs -n1 | awk 'NR==1' | tr ',' ' '))


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


DIR_COMP=${DIR_DATA}/${NAME}/tmp_Compartment_${RESOLUTION}
[ ! -e ${DIR_COMP} ] && mkdir $DIR_COMP
for i in $(seq 1 ${#CHRs[@]})
do
	let index=i-1
	CHR=${CHRs[index]}
	[ ! -e ${DIR_COMP}/${CHR}.txt ] && Rscript --vanilla --slave ${DIR_LIB}/PCA_analysis.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds --geneDensity ${DIR_contig}/GENE_density_${RESOLUTION}.txt --location ${DIR_COMP}/${CHR}.txt
done

cat ${DIR_COMP}/*.txt | sort -k1,1 -k2,2n > ${DIR_DATA}/${NAME}/Compartment_${RESOLUTION}.txt

if [ "${RESOLUTION}" = "40kb" ]; then
    Rscript --vanilla --slave ${DIR_LIB}/Modify_compartment.R -i ${DIR_DATA}/${NAME}/Compartment_${RESOLUTION}.txt -o ${DIR_COMP}/Corrected_Compartment.txt
    mv ${DIR_COMP}/Corrected_Compartment.txt ${DIR_DATA}/${NAME}/Compartment_${RESOLUTION}.txt
fi
rm -r ${DIR_COMP}




