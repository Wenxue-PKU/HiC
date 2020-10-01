#!/bin/bash
# intra/inter ratio

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version

	r, --resolution
		bin size for analysis (ex. 10kb)

	-i, --in [map file]
		map file

	-o, --out [output directory]
		output directory

	-c, --chromosome [chromosome length file]
		chromosome file
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvr:i:o:c:
LONG=help,version,resolution:,in:,out:,chromosome:
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
		-r|--resolution)
			RESOLUTION_string="$2"
			shift 2
			;;
		-i|--in)
			FILE_map="$2"
			shift 2
			;;
		-o|--out)
			DIR_out="$2"
			shift 2
			;;
		-c|--chromosome)
			CHROME="$2"
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
TIME_STAMP=$(date +"%Y-%m-%d_%H.%M.%S")
INPUT_FILES=$@

[ ! -n "${FILE_map}" ] && echo "Please specify input map file" && exit 1
[ ! -n "${DIR_out}" ] && echo "Please specify output directory" && exit 1
[ ! -n "${CHROME}" ] && echo "Please specify chromsome length file" && exit 1
[ ! -n "${RESOLUTION_string}" ] && echo "Please specify resolution" && exit 1
[ ! -e ${DIR_out} ] && mkdir ${DIR_out}

RESOLUTION=${RESOLUTION_string/kb/000}


#==============================================================
# count score
#==============================================================
zcat $FILE_map | tail -n +2 | awk -v OFS='\t' -v BIN=$RESOLUTION '$2==$9&&($10-$3)>10000{
	bin1=int($3/BIN)*BIN;
	bin2=int($10/BIN)*BIN;
	INTRA[$2"\t"bin1]+=1;
	INTRA[$9"\t"bin2]+=1;
}$2==$9&&($10-$3)>10000&&($10-$3)<100000{
	bin1=int($3/BIN)*BIN;
	bin2=int($10/BIN)*BIN;
	LOCAL[$2"\t"bin1]+=1;
	LOCAL[$9"\t"bin2]+=1;
}$2!=$9{
	bin1=int($3/BIN)*BIN;
	bin2=int($10/BIN)*BIN;
	INTER[$2"\t"bin1]+=1;
	INTER[$9"\t"bin2]+=1;
}END{
	for(x in INTRA){
		if( !(x in LOCAL)){LOCAL[x]=0}
		if( !(x in INTER)){INTER[x]=0}
		print x,LOCAL[x],INTRA[x],INTER[x]
	}
}' | sort -k1,1 -k2,2n >> ${DIR_out}/tmp_bin_${RESOLUTION_string}.txt


### filter chromosome size
echo "chr start end local intra inter" | tr ' ' '\t' > ${DIR_out}/count_bin_${RESOLUTION_string}.txt
awk -v OFS='\t' -v BIN=$RESOLUTION 'NR==FNR{
	LEN[$1]=$2;
	next;
}{
	end=$2+BIN-1
	if(LEN[$1]<end){
		end=LEN[$1]
	}
	print $1,$2,end,$3,$4,$5
}' $CHROME ${DIR_out}/tmp_bin_${RESOLUTION_string}.txt >> ${DIR_out}/count_bin_${RESOLUTION_string}.txt
rm ${DIR_out}/tmp_bin_${RESOLUTION_string}.txt


#==============================================================
# convert to bigwig
#==============================================================
cd ${DIR_out}
module load easybuild 2> /dev/null && module load Kent_tools

### inter percentage
# awk -v OFS='\t' 'NR>1&&($5+$6)!=0{
# 	score=$6/($5+$6);
# 	print $1,$2,$3,score
# }' count_bin_${RESOLUTION_string}.txt > inter_percent.bedgraph
# bedGraphToBigWig inter_percent.bedgraph $CHROME inter_percent_bin_${RESOLUTION_string}.bw && rm inter_percent.bedgraph

### inter/less than 100kb
# awk -v OFS='\t' 'NR>1&&$4!=0{
# 	score=$6/$4;
# 	print $1,$2,$3,score
# }' count_bin_${RESOLUTION_string}.txt > ratio_inter_local_100kb.bedgraph
# bedGraphToBigWig ratio_inter_local_100kb.bedgraph $CHROME ratio_inter_local_100kb_bin_${RESOLUTION_string}.bw && rm ratio_inter_local_100kb.bedgraph


### less than 100kb/inter
awk -v OFS='\t' 'NR>1&&$6!=0{
	score=$4/$6;
	print $1,$2,$3,score
}' count_bin_${RESOLUTION_string}.txt > ratio_local_100kb_inter_bin_${RESOLUTION_string}.bedgraph
bedGraphToBigWig ratio_local_100kb_inter_bin_${RESOLUTION_string}.bedgraph $CHROME ratio_local_100kb_inter_bin_${RESOLUTION_string}.bw && rm ratio_local_100kb_inter_bin_${RESOLUTION_string}.bedgraph







