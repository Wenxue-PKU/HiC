#!/bin/bash
# gene densityを調べる

get_usage(){
	cat <<EOF

Usage : $0 [OPTION] -r [resolution] -x [organism]

Description
	-h, --help
		show help

	-v, --version
		show version

	-r, --resolution [xxx kb]
		resolution. Add kb at the end

	-x, --organism [organism]
		pombe, human, mouse

EOF
}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvr:x:
LONG=help,version,resolution:,organism:
PARSED=`getopt --options SHORT --longoptions $LONG --name "$0" -- "$@"`
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
            RESOLUTION="$2"
            shift 2
            ;;
        -x|--organism)
            ORGANISM="$2"
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

RESOLUTION_NUM=${RESOLUTION/kb/000}

case $ORGANISM in
	pombe) DIR_contig=/wistar/noma/Data/S.Pombe_seq/pombase_ASM294v2.30 ;;
	human) DIR_contig=/wistar/noma/Data/Human_seq/hg19 ;;
	mouse) DIR_contig=/wistar/noma/Data/Mouse_seq/mm10 ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac

FILE_genes=${DIR_contig}/contig/Gene_annotation_from_UCSC.bed
FILE_SECTION=${DIR_contig}/${RESOLUTION}_section.bed

DIR_PWD=$(dirname $0)
Rscript --vanilla --slave ${DIR_PWD}/Make_section.R --resolution $RESOLUTION_NUM --organism $ORGANISM --out $FILE_SECTION


# bed toolで数を数える
FILE_geneDensity=${DIR_contig}/GENE_density_${RESOLUTION}.txt
bedtools intersect -c -wa -a $FILE_SECTION -b $FILE_genes | awk '{print $1":"$2":"$3"\t"$4}' > $FILE_geneDensity
rm $FILE_SECTION



