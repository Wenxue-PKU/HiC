#!/bin/bash
# download fastq file from SRA

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-n, --name [sample name]
		sample name

	-s, --sra [SRA ID]
		SRA ID (ex. SRR2601012)

	-d, --directory [directory]
		download directory

	-p, --pair [TRUE|FALSE]
		paired-end file (TRUE) or not (FALSE). Default (FALSE)
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvn:s:d:p:
LONG=help,version,name:,sra:,directory:,pair:
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
		-n|--name)
			NAME="$2"
			shift 2
			;;
		-s|--sra)
			SRA="$2"
			shift 2
			;;
		-d|--directory)
			DIR_DATA="$2"
			shift 2
			;;
		-p|--pair)
			FLAG_PAIR="$2"
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

[ ! -n "${NAME}" ] && echo "Please specify NAME" && exit 1
[ ! -n "${SRA}" ] && echo "Please specify SRA id" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify directory" && exit 1
[ ! -e ${DIR_DATA} ] && echo "Data directory is not exists" && exit 1
FLAG_PAIR=${FLAG_PAIR:-FALSE}

URL=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${SRA:0:6}/${SRA}/${SRA}.sra
echo "Download file from $URL"

cd ${DIR_DATA}; 
while [ 1 ]; do
    wget -O ${NAME}.sra ${URL} --no-verbose --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 0 --continue
    if [ $? = 0 ]; then break; fi; # check return value, break if successful (0)
    sleep 1s;
done;


if [ "$FLAG_PAIR" = "TRUE" ]; then
	/applications/sratoolkit/current/fastq-dump --split-files ${NAME}.sra
else
	/applications/sratoolkit/current/fastq-dump ${NAME}.sra
fi

