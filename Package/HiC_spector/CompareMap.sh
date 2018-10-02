#!/bin/bash
# compare two matrix

LOOT=/wistar
[ -e /cscb ] && LOOT=/cscb

DIR_contig=${LOOT}/noma/Data/S.Pombe_seq/pombase_ASM294v1.18
FILE_restriction=${DIR_contig}/MboI_sites_for_juicer.txt
FILE_chromsize=${DIR_contig}/chrom_sizes_for_juicer.txt
PRO_PRE=${LOOT}/noma/Program/juicer_tools_0.7.0.jar



get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	-t, --tmp_dir [directory to keep temp files]
		directory to place tmp files

EOF
}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvt:
LONG=help,version,tmp_dir:
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
        -t|--tmp_dir)
            DIR_TMP="$2"
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



#==============================================================
# output



#==============================================================
# convert to .hic file
#==============================================================
java -Xmx2g -jar $PRO_PRE pre -f $FILE_restriction -r 5000,10000,20000 $FILE_preMAP1 $FILE_HiC1 $FILE_chromsize






