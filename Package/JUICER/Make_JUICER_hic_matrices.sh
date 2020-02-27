#!/bin/bash
# Prepare Juicer's hic matrices
# Hideki Tanizawa (rafysta@gmail.com)

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-i, --in [map.gz file]
		NAME.map.gz file. For muptiple input file, make space between and name and quote. ex. ("a.map.gz b.map.gz")

	-o, --out [output hic file]
		sample name

	-r, --resolution [resolution list]
		resolution list separated by ,. For example 5000,10000,20000,100000 or 2f,5f,10f for fragment resolution

	-f, --restriction [HindIII|MboI|MboI-HinfI]
		name for restriction if prepare fragment resolution data

	-q, --mapQ [mapQ threshold 30 (default)]
		MapQ threshold. Defulat is 30.
	
	-c, --chromosome [only calculate specific chromosome]
		Only calculate the specified chromosome's intra matrices

	-x, --ref [ex. hg19]
		organism reference sequence name
	
	--old_map
		if using old hic pipeline
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hi:o:r:f:q:c:x:
LONG=help,in:,out:,resolution:,restriction:,mapQ:,chromosome:,ref:,old_map
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
		-c|--chromosome)
			TARGET_CHR="$2"
			shift 2
			;;
		-x|--ref)
			REF="$2"
			shift 2
			;;
		--old_map)
			FLAG_oldmap="TRUE"
			shift 1
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

[ ! -n "${FILE_MAP}" ] && echo "Please specify map file" && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output hic file name" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "${REF}" ] && echo "Please specify ref" && exit 1
THRESHOLD_MAPQ=${CHR_include:-30}
RESTRICTION=${RESTRICTION:-NA}
FLAG_oldmap=${FLAG_oldmap:-FALSE}
TARGET_CHR=${TARGET_CHR:-NA}

#-----------------------------------------------
# Load setting
#-----------------------------------------------
source ${DIR_LIB}/../../utils/load_setting.sh -x $REF -r $RESTRICTION


#-----------------------------------------------
# All calculation will be done in tmp directory
#-----------------------------------------------
DIR_tmp=$(mktemp -d /tmp/tmp_hic_juicer.XXXXXX)
FILE_data=data_for_hic.txt
FILE_hic=matrix.hic
trap "rm -r ${DIR_tmp}" 0
cd ${DIR_tmp}


### Juicer restriction site file
FLAG_RESOLUTION=$(echo ${RESOLUTION} | grep -c f)
FILE_enzyme_def=${FILE_enzyme_def/.txt/.juicer.txt}
if [ $FLAG_RESOLUTION -eq 1 ]; then
	[ "$RESTRICTION" = "NA" ] && echo "Please specify restriction file if calculate fragment resolution" && exit 1
	[ ! -e $FILE_enzyme_def ] && echo "Restriction file not found. Please prepare it using Prepare_JUICER_restriction_files_for_all_data.sh" && exit 1
	cp $FILE_enzyme_def restriction_sites.txt
fi
FILE_enzyme_def=restriction_sites.txt


### JUICER program
# Download from https://github.com/aidenlab/juicer/wiki/Download
PROGRAM_JUICER=${HOME}/Software/juicebox_tools.jar
[ ! -e $PROGRAM_JUICER ] && echo "juicer program not found" && exit 1
cp $PROGRAM_JUICER juicebox_tools.jar
PROGRAM_JUICER=juicebox_tools.jar


### Chromosome file
if [ "$TARGET_CHR" != "NA" ]; then
	cat $FILE_CHROME_LENGTH | grep "$TARGET_CHR" > chromosome.txt
else
	cp $FILE_CHROME_LENGTH chromosome.txt
fi
FILE_CHROME_LENGTH=chromosome.txt


#==============================================================
# Convert map file to JUICER input format
#==============================================================
if [ "$TARGET_CHR" != "NA" ]; then
	if [ "$FLAG_oldmap" = "TRUE" ]; then
		### for old custom hic pipeline (use XXX_sort.map.gz)
		zcat ${FILE_MAP} | awk -v tc=$TARGET_CHR '$8=="U" && $15=="U" && $2==tc && $9==tc{
			if($4=="+"){str1=0}else{str1=1}
			if($11=="+"){str2=0}else{str2=1}
			$6=gsub("L","",$6); $6=gsub("R","",$6);
			$13=gsub("L","",$13); $13=gsub("R","",$13);
			print $1,str1,$2,$3,$6,str2,$9,$10,$13,$5,$12
		}' | sort -k3,3d -k7,7d > ${FILE_data}
	else
		zcat ${FILE_MAP} | awk -v tc=$TARGET_CHR 'NR>1 && $8=="U" && $15=="U" && $2==tc && $9==tc{
			if($4=="+"){str1=0}else{str1=1}
			if($11=="+"){str2=0}else{str2=1}
			print $1,str1,$2,$3,$6,str2,$9,$10,$13,$5,$12
		}' | sort -k3,3d -k7,7d > ${FILE_data}
	fi
else
	if [ "$FLAG_oldmap" = "TRUE" ]; then
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
fi


#==============================================================
# JUICER calculation
#==============================================================
if hash module 2>/dev/null; then
	module load java
fi

if [ "$TARGET_CHR" != "NA" ]; then
	if [ $FLAG_RESOLUTION -eq 0 ]; then
		java -Xms512m -Xmx2048m -jar $PROGRAM_JUICER pre -r ${RESOLUTION} -d -c $TARGET_CHR -q $THRESHOLD_MAPQ -t ${DIR_tmp} ${FILE_hic} $FILE_CHROME_LENGTH
	else
		java -Xms512m -Xmx2048m -jar $PROGRAM_JUICER pre -r ${RESOLUTION} -d -c $TARGET_CHR -f $FILE_enzyme_def -q $THRESHOLD_MAPQ -t ${DIR_tmp} ${FILE_data} ${FILE_hic} $FILE_CHROME_LENGTH
	fi
else
	if [ $FLAG_RESOLUTION -eq 0 ]; then
		java -Xms512m -Xmx2048m -jar $PROGRAM_JUICER pre -r ${RESOLUTION} -q $THRESHOLD_MAPQ -t ${DIR_tmp} ${FILE_data} ${FILE_hic} $FILE_CHROME_LENGTH
	else
		java -Xms512m -Xmx2048m -jar $PROGRAM_JUICER pre -r ${RESOLUTION} -f $FILE_enzyme_def -q $THRESHOLD_MAPQ -t ${DIR_tmp} ${FILE_data} ${FILE_hic} $FILE_CHROME_LENGTH
	fi
fi

mv $FILE_hic ${FILE_OUT}

exit 0



