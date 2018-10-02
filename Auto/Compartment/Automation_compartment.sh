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

	-a, --assemble [hg19|mm10]
		assemble name

	-r, --resolution [200kb|500kb]
		default : 200kb

	-f, --force [FALSE(default)|TRUE]
		force to perform analysis even matrices is in-complete or not (default)

EOF

}


SHORT=hd:n:a:r:f:
LONG=help,directory:,name:,assemble:,resolution:,force:
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
        -a|--assemble)
            ASSEMBLE="$2"
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

TIME_STAMP=$(date +"%Y-%m-%d")
DIR_LIB=$(dirname $0)

case $ASSEMBLE in
	hg19) DIR_contig=/wistar/noma/Data/Human_seq/hg19
#		  CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY) ;;
          CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX) ;;
	mm10) DIR_contig=/wistar/noma/Data/Mouse_seq/mm10
#		  CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY) ;;
		  CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX) ;;
	*)    echo "Please specify correct organism"
	      eixt 1 ;;
esac



MISSING=()
for i in `seq 1 ${#CHRs[@]}`
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
[ -e ${DIR_COMP} ] && rm -r ${DIR_COMP}
mkdir $DIR_COMP
for i in `seq 1 ${#CHRs[@]}`
do
	let index=i-1
	CHR=${CHRs[index]}
	Rscript --vanilla --slave ${DIR_LIB}/PCA_analysis.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds --geneDensity ${DIR_contig}/GENE_density_${RESOLUTION}.txt --location ${DIR_COMP}/${CHR}.txt
done

cat ${DIR_COMP}/*.txt | sort -k1,1 -k2,2n > ${DIR_DATA}/${NAME}/Compartment_${RESOLUTION}.txt

if [ "${RESOLUTION}" = "40kb" ]; then
    Rscript --vanilla --slave ${DIR_LIB}/Modify_compartment.R -i ${DIR_DATA}/${NAME}/Compartment_${RESOLUTION}.txt -o ${DIR_COMP}/Corrected_Compartment.txt
    mv ${DIR_COMP}/Corrected_Compartment.txt ${DIR_DATA}/${NAME}/Compartment_${RESOLUTION}.txt
fi
rm -r ${DIR_COMP}




