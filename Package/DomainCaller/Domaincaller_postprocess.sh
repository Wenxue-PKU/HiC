#!/bin/bash
# Post processing of Domainc calling

get_usage(){
    cat <<EOF

Usage : $0 [OPTION]

Description
    -h, --help
        show help

    -v, --version
        show version

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

get_version(){
    echo "${0} version 1.0"
}

SHORT=hvd:n:a:r:
LONG=help,version,directory:,name:,assemble:,resolution:
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
            RESOLUTION_string="$2"
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

DIR_LIB=$(dirname $0)
DIR_LIB_perl=/wistar/noma/Program/DomainCaller/perl_scripts
TIME_STAMP=$(date +"%Y-%m-%d")

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${ASSEMBLE}" ] && echo "Please specify assemble version" && exit 1
RESOLUTION_string=${RESOLUTION_string:-40kb}
FLAG_FORCE=${FLAG_FORCE:-FALSE}

RESOLUTION=${RESOLUTION_string/kb/000}

case $ASSEMBLE in
	hg19) DIR_contig=/wistar/noma/Data/Human_seq/hg19
#		  CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY) ;;
		#   CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX) ;;
          CHRs=(chr1 chr2 chr3 chr4) ;;
	mm10) DIR_contig=/wistar/noma/Data/Mouse_seq/mm10
#		  CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY) ;;
		  CHRs=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX) ;;
	*)    echo "Please specify correct organism"
	      eixt 1 ;;
esac

CHROM_SIZE=${DIR_contig}/LENGTH.txt
DIR_DOMAIN=${DIR_DATA}/${NAME}/${RESOLUTION_string}/Domain


### Matlabの結果を確認
MISSING=()
for i in `seq 1 ${#CHRs[@]}`
do
	let index=i-1
	CHR=${CHRs[index]}
	[ ! -e ${DIR_DOMAIN}/${CHR}.hmm.txt ] && MISSING=("${MISSING[@]}" $CHR)
done
if [ ${#MISSING[@]} -gt 0 ]; then
	echo "Matlab result of ${NAME} ${MISSING[@]} are not exists"
	[ "$FLAG_FORCE" != "TRUE" ] && echo "quit job" && exit 1
fi


for i in `seq 1 ${#CHRs[@]}`
do
	let index=i-1
	CHR=${CHRs[index]}
    perl ${DIR_LIB_perl}/file_ends_cleaner.pl ${DIR_DOMAIN}/${CHR}.hmm.txt ${DIR_DOMAIN}/${CHR}.DI.txt | ${DIR_LIB_perl}/converter_7col.pl > ${DIR_DOMAIN}/${CHR}.DI.7col.txt
    perl ${DIR_LIB_perl}/hmm_probablity_correcter.pl ${DIR_DOMAIN}/${CHR}.DI.7col.txt 2 0.99 ${RESOLUTION} | perl ${DIR_LIB_perl}/hmm-state_caller.pl $CHROM_SIZE ${CHR} | perl ${DIR_LIB_perl}/hmm-state_domains.pl > ${DIR_DOMAIN}/${CHR}.domain.txt
done

cat ${DIR_DOMAIN}/*.domain.txt | sort -k1,1 -k2,2n > ${DIR_DATA}/${NAME}/Domain_${RESOLUTION_string}.txt


