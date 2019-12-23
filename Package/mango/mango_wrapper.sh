#!/bin/bash
# Mango analysis

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-o, --out [output directory]
		output directory name

	-n, --name [sample name]
		sample name
		
	-x, --ref [ex. hg19]
		organism name
	
	--fastq1 [fastq]
		fasta file name
	
	--fastq2 [fastq]
		fasta file name
	
	--chrominclude [include chromosome.]
		chromosome list separated by commma,

	--maxDistance [maximim distance between pairs]
		maximum distance between pairs (default: 1000000)
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvd:o:n:x:
LONG=help,version,out:,name:,ref:,fastq1:,fastq2:,chrominclude:,maxDistance:
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
		-o|--out)
			DIR_OUT="$2"
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
		--fastq1)
			FILE_fastq1="$2"
			shift 2
			;;
		--fastq2)
			FILE_fastq2="$2"
			shift 2
			;;
		--chrominclude)
			CHROM_INCLUDE="$2"
			shift 2
			;;
		--maxDistance)
			MAX_DISTANCE="$2"
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

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1
[ ! -n "${FILE_fastq1}" ] && echo "Please specify input fastq1 file" && exit 1
[ ! -n "${FILE_fastq2}" ] && echo "Please specify input fastq2 file" && exit 1
[ ! -n "${REF}" ] && echo "Please specify ref" && exit 1
CHROM_INCLUDE=${CHROM_INCLUDE:-NULL}
MAX_DISTANCE=${MAX_DISTANCE:-1000000}

#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/../../utils/load_setting.sh -x $REF -r NA


if hash module 2>/dev/null; then
	module load miniconda
	conda activate mango-deps
fi

DIR_tmp=$(mktemp -d /tmp/tmp_${NAME}.XXXXX)

cd ${DIR_tmp}
ln -s ${FILE_fastq1} ${NAME}_1.same.fastq
ln -s ${FILE_fastq2} ${NAME}_2.same.fastq


### Mango
Rscript --vanilla --slave ${HOME}/Software/mango/mango.R --fastq1 ${FILE_fastq1} --fastq2 ${FILE_fastq2} --prefix ${NAME} --chrominclude ${CHROM_INCLUDE} --bedtoolsgenome $FILE_CHROME_LENGTH --stages 2:5 --outdir ${DIR_tmp} --bowtieref $BOWTIE_INDEX --keepempty TRUE --reportallpairs TRUE --maxinteractingdist $MAX_DISTANCE


rm ${NAME}_?.same.fastq
rm *.sam *.bedpe *.tagAlign
mv ${DIR_tmp}/* ${DIR_OUT}

rm -r ${DIR_tmp}

exit 0
