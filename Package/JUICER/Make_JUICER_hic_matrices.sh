#!/bin/bash
# Prepare Juicer's hic matrices

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-i, --in [map.gz file]
		NAME.map.gz file

	-o, --out [output hic file]
		sample name

	-r, --resolution [resolution list]
		resolution list separated by ,. For example 5000,10000,20000,100000 or 2f,5f,10f for fragment resolution

	-f, --restriction [HindIII|MboI|MboI-HinfI]
		name for restriction if prepare fragment resolution data (currently not implemented)

	-q, --mapQ [mapQ threshold 30 (default)]
		MapQ threshold. Defulat is 30.

	-x, --ref [ex. hg19]
		organism reference sequence name
	
	--old_map
		if using old hic pipeline
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hi:o:r:f:q:x:
LONG=help,in:,out:,resolution:,restriction:,mapQ:,ref:,old_map
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
		-i|--in)
			FILE_MAP="$2"
			shift 2
			;;
		-o|--out)
			FILE_OUT="$2"
			shift 2
			;;
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
			;;
		-f|--restriction)
			RESTRICTION="$2"
			shift 2
			;;
		-q|--mapQ)
			THRESHOLD_MAPQ="$2"
			shift 2
			;;
		-x|--ref)
			REF="$2"
			shift 2
			;;
		--old_map)
			FLAG_oldmap=TRUE
			exit 1
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

[ ! -n "${FILE_MAP}" ] && echo "Please specify map file" && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output hic file name" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "$REF{}" ] && echo "Please specify ref" && exit 1
THRESHOLD_MAPQ=${CHR_include:-30}
RESTRICTION=${RESTRICTION:-NA}
FLAG_oldmap=${FLAG_oldmap:-FALSE}

# wget http://hicfiles.s3.amazonaws.com/internal/juicebox_tools/8.5.16/juicebox_tools.jar
PROGRAM_JUICER=${HOME}/Software/juicebox_tools.jar

[ ! -e $PROGRAM_JUICER ] && echo "juicer program not found" && exit 1


#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/../../utils/load_setting.sh -x $REF -r $RESTRICTION

### Create Temporary directory
DIR_tmp=$(mktemp -d /tmp/tmp_hic_juicer.XXXXXX)
FILE_data=${DIR_tmp}/data_for_hic.txt
FILE_hic=${DIR_tmp}/matrix.list
trap "rm -r ${DIR_tmp}" 0


#==============================================================
# Convert matrices
#==============================================================
if [ $FLAG_oldmap = "TRUE" ]; then
	### for old custom hic pipeline (use XXX_sort.map.gz)
	zcat ${FILE_MAP} | awk '$8=="U" && $15=="U"{
		if($4=="+"){str1=0}else{str1=1}
		if($11=="+"){str2=0}else{str2=1}
		$6=gsub("L","",$6); $6=gsub("R","",$6);
		$13=gsub("L","",$13); $13=gsub("R","",$13);
		print $1,str1,$2,$3,$6,str2,$9,$10,$13,$5,$12
	}' | sort -k3,3d -k7,7d > ${FILE_data}
else
	zcat ${FILE_MAP} | awk 'NR>1 && $8=="U" && $15=="U"{
		if($4=="+"){str1=0}else{str1=1}
		if($11=="+"){str2=0}else{str2=1}
		print $1,str1,$2,$3,$6,str2,$9,$10,$13,$5,$12
	}' | sort -k3,3d -k7,7d > ${FILE_data}
fi

java -Xmx8g -jar $PROGRAM_JUICER pre -r ${RESOLUTION} -q $THRESHOLD_MAPQ -t ${DIR_tmp} ${FILE_data} ${FILE_hic} $FILE_CHROME_LENGTH
mv $FILE_hic ${FILE_OUT}

exit 0



