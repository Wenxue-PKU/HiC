#!/bin/bash
# EBV 4C generation

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
	
	-d, --data [directory]
		data directory
	
	-r, --resolution [resolution ex. 10kb]
		resolution for analysis (default : 10kb)

	-o, --out [output directory]
		output directory	

EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvn:d:r:o:
LONG=help,version,name:,data:,resolution:,out:
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
		-d|--data)
			DIR_DATA="$2"
			shift 2
			;;
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
			;;
		-o|--out)
			DIR_OUT="$2"
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
[ ! -n "${DIR_DATA}" ] && echo "Please specify data" && exit 1
[ ! -n "${DIR_OUT}" ] && echo "Please specify output directory" && exit 1
RESOLUTION=${RESOLUTION:-10kb}


CHROM_SIZE=/wistar/noma/Data/Human_seq/hg19_EBV/LENGTH.txt

cd ${DIR_OUT}

[ ! -e window_${RESOLUTION}.bed ] && bedtools makewindows -w ${RESOLUTION/kb/000} -g $CHROM_SIZE -i winnum > window_${RESOLUTION}.bed

#==============================================================
# windowごとの合計scoreを計算する
#==============================================================
RESOLUTION=10kb
awk -v OFS='\t' -v R=${RESOLUTION/kb/000} '
	$1=="EBV"{
		m=int(($6+$7)/2/R)*R; 
		count[$5"\t"m]+=$9
	}
	$5=="EBV"{
		m=int(($2+$3)/2/R)*R;
		count[$1"\t"m]+=$9
	}
	END{
		for (i in count){
			print i, count[i];
		}
	}' ${DIR_DATA}/${NAME}_dataForFragDb.txt | LC_COLLATE=C sort -k1,1 -k2,2n | awk -v OFS='\t' -v R=${RESOLUTION/kb/000} '{e=$2+R-1; print $1,$2,e,$3}' > ${NAME}_score.bedgraph



#==============================================================
# bigwigの作成
#==============================================================
awk -v OFS='\t' '
	NR==FNR{
		len[$1]=$2; next
	}
	{
		if($3 < len[$1]){
			print;
		}
	}' $CHROM_SIZE ${NAME}_score.bedgraph > ${NAME}_score_trimed.bedgraph
/applications/bedgraphtobigwig/current/bedGraphToBigWig ${NAME}_score_trimed.bedgraph $CHROM_SIZE ${NAME}_4C.bw


#==============================================================
# define significant peaks
#==============================================================
Rscript --vanilla --slave ${DIR_LIB}/Pvalue.R -i ${NAME}_score_trimed.bedgraph -o ${NAME}_pval.bedgraph

### macs2でpeak calling
grep -v "EBV" ${NAME}_pval.bedgraph > ${NAME}_pval_noEBV.bedgraph
macs2 bdgpeakcall -i ${NAME}_pval_noEBV.bedgraph -c 5 --no-trackline -l 20000 -g 10000 -o ${NAME}_macs2.txt && cat ${NAME}_macs2.txt | awk -v OFS='\t' '{POS=$2+$10; print $1,$2,$3,"peak"NR,$5,POS}' | sort -k1,1 -k2,2n > ${NAME}_peaks.bed
rm ${NAME}_pval_noEBV.bedgraph

cat ${NAME}_peaks.bed | awk -v OFS='\t' '{print $1,$2,$3,$4}' > ${NAME}_peaks_for_bb.txt
/cm/shared/wistar/kentsource/339/bin/bedToBigBed ${NAME}_peaks_for_bb.txt $CHROM_SIZE ${NAME}_peaks.bb
rm ${NAME}_peaks_for_bb.txt

### 10kb windowごとにpeakがあるかどうかを出力
cat window_${RESOLUTION}.bed | grep -v EBV | cut -f1-3 | bedtools intersect -c -a stdin -b ${NAME}_peaks.bed > ${NAME}_peaks.bedgraph
