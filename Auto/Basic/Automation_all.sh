#!/bin/bash
# Automation Hi-C analysis

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

	-x, --organism [human|human_EBV|pombe|mouse]
		organism name

	-r, --restriction [HindIII|MboI|MboI-HinfI]
		name for restriction

	-o, --log [log directory]
		log file directory

	-m, --mapq [mapq threshold (default:10)]
		threshold mapQ to make map

	-q, --fastqc
		TRUE for doing fastqc analysis or FALSE for not doing (default TRUE)

EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvd:n:x:r:o:m:q:
LONG=help,version,directory:,name:,organism:,restriction:,log:,mapq:,fastqc:
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
		-x|--organism)
			ORGANISM="$2"
			shift 2
			;;
		-r|--restriction)
			RESTRICTION="$2"
			shift 2
			;;
		-o|--log)
			DIR_LOG="$2"
			shift 2
			;;
		-m|--mapq)
			MAPQ_THRESHOLD="$2"
			shift 2
			;;
		-q|--fastqc)
			FLAG_fastqc="$2"
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
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
[ ! -n "${RESTRICTION}" ] && echo "Please specify restriction" && exit 1
[ ! -n "${DIR_LOG}" ] && echo "Please specify log directory" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
MAPQ_THRESHOLD=${MAPQ_THRESHOLD:-10}
FLAG_fastqc=${FLAG_fastqc:-TRUE}

case $ORGANISM in
	pombe)	BOWTIE_TARGET=pombe
			BOWTIE2_INDEXES=/wistar/bioinfo-nfs/hideki_projects/Genome/data/pombe/2018
			CHROM_LENGTH=12571820
			FILE_CHROME_LENGTH=/wistar/bioinfo-nfs/hideki_projects/Genome/data/pombe/2018/LENGTH.txt
			case $RESTRICTION in 
				MboI)	FILE_enzyme_index=/wistar/bioinfo-nfs/hideki_projects/Genome/data/pombe/2018/Sectioning_MboI.txt
						FILE_enzyme_def=/wistar/bioinfo-nfs/hideki_projects/Genome/data/pombe/2018/MboI_sites.txt ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	human)	BOWTIE_TARGET=hg19
			BOWTIE2_INDEXES=/wistar/noma/Data/Human_seq/hg19
			CHROM_LENGTH=3095677412
			FILE_CHROME_LENGTH=/wistar/noma/Data/Human_seq/hg19/LENGTH.txt
			case $RESTRICTION in 
				HindIII)	FILE_enzyme_index=/wistar/noma/Data/Human_seq/hg19/Sectioning_HindIII.txt
							FILE_enzyme_def=/wistar/noma/Data/Human_seq/hg19/HindIII_sites.txt ;;
				MboI)	FILE_enzyme_index=/wistar/noma/Data/Human_seq/hg19/Sectioning_MboI.txt
						FILE_enzyme_def=/wistar/noma/Data/Human_seq/hg19/MboI_sites.txt ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	human_EBV)	BOWTIE_TARGET=hg19_EBV
			BOWTIE2_INDEXES=/wistar/noma/Data/Human_seq/hg19_EBV
			CHROM_LENGTH=3157782322
			FILE_CHROME_LENGTH=/wistar/noma/Data/Human_seq/hg19_EBV/LENGTH.txt
			case $RESTRICTION in 
				MboI)	FILE_enzyme_index=/wistar/noma/Data/Human_seq/hg19_EBV/Sectioning_MboI.txt
						FILE_enzyme_def=/wistar/noma/Data/Human_seq/hg19_EBV/MboI_sites.txt;;
				MboI-HinfI)	FILE_enzyme_index=/wistar/noma/Data/Human_seq/hg19_EBV/Sectioning_MboI-HinfI.txt
						FILE_enzyme_def=/wistar/noma/Data/Human_seq/hg19_EBV/MboI-HinfI_sites.txt ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	mouse)	BOWTIE_TARGET=mm10
			BOWTIE2_INDEXES=/wistar/noma/Data/Mouse_seq/mm10
			CHROM_LENGTH=2725537669
			FILE_CHROME_LENGTH=/wistar/noma/Data/Mouse_seq/mm10/LENGTH.txt
			case $RESTRICTION in 
				MboI)	FILE_enzyme_index=/wistar/noma/Data/Mouse_seq/mm10/Sectioning_MboI.txt
						FILE_enzyme_def=/wistar/noma/Data/Mouse_seq/mm10/MboI_sites.txt ;;
				*)	echo "$RESTRICTION is not registered for $ORGANISM"
					exit ;;
			esac
			;;
	*)	echo "Please specify correct organism"
		eixt 1 ;;
esac


# fastqcディレクトリが存在しなかったら作成
[ "$FLAG_fastqc" = "TRUE" ] && [ ! -e "${DIR_DATA}/fastqc" ] && mkdir "${DIR_DATA}/fastqc"


#-----------------------------------------------
# Alignment
#-----------------------------------------------
sbatch -n 12 --job-name=aln_${NAME} $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_Alignment_${NAME}_1.log" --export=NAME="${NAME}_1",DIR_LIB="${DIR_LIB}",DIR_DATA="${DIR_DATA}",BOWTIE_TARGET="${BOWTIE_TARGET}",BOWTIE2_INDEXES="${BOWTIE2_INDEXES}" --open-mode append ${DIR_LIB}/Alignment_with_trimming.sh
sbatch -n 12 --job-name=aln_${NAME} $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_Alignment_${NAME}_2.log" --export=NAME="${NAME}_2",DIR_LIB="${DIR_LIB}",DIR_DATA="${DIR_DATA}",BOWTIE_TARGET="${BOWTIE_TARGET}",BOWTIE2_INDEXES="${BOWTIE2_INDEXES}" --open-mode append ${DIR_LIB}/Alignment_with_trimming.sh


#-----------------------------------------------
# fastqc
#-----------------------------------------------
[ "$FLAG_fastqc" = "TRUE" ] &&  [ ! -e ${DIR_DATA}/fastqc/${NAME}_fastqc ] && sbatch -n 12 --job-name=fastqc_${NAME}_1 $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_fastqc_${NAME}_1.log" --open-mode append --wrap="cd ${DIR_DATA}; /applications/fastqc/current/fastqc -o fastqc/ --nogroup -t 12 ${NAME}_1.fastq" && sbatch -n 12 --job-name=fastqc_${NAME}_2 -o "${DIR_LOG}/${TIME_STAMP}_fastqc_${NAME}_2.log" --open-mode append --wrap="cd ${DIR_DATA}; /applications/fastqc/current/fastqc -o fastqc/ --nogroup -t 12 ${NAME}_2.fastq"



### assign to nearest restriction enzyme site
# output: <NAME>.map
# +の向きの場合、そのまま、-の向きの場合、aglinした部位からreadの長さ分だけ足した値に修正する
# RepeatとUniqueは、XS:iがあるか無いかで判断
# 最も近い制限酵素部位を出力する
# +の向きの場合、制限酵素からの位置に関わらずL,-の向きの場合,Rと出力する
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "aln_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=map_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_make_map_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Assign_nearest_enzymeSites.pl -a ${NAME}_1.sam -b ${NAME}_2.sam -o ${NAME}.map -e ${FILE_enzyme_def} -d ${FILE_enzyme_index}"



#-----------------------------------------------
# Sort map file
#-----------------------------------------------
### split files
# 染色体全体を100に分割し、左のreadの場所に応じてファイルを約100個出力する。ファイルの名前を<NAME>_list.txtとして出力
# 片方でもalignできなかったもの(制限酵素を特定できなかったものも含めて)は除去
# 右と左を比べて左の方が小さくなるように入れ替える
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "map_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=split_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_sort_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Split_MapFile.pl -i ${NAME}.map -l ${CHROM_LENGTH} -o ${NAME}_list.txt"



### sorting each files
# <NAME>_list.txtに書いてあるファイルを読みこんで、sortコマンド(sort -k 3,3n -k 4,4)を実行する(chromosomeは既に同じことが保証されているので、座標だけをソートする)
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "split_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch --array=1-200 --job-name=sort_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_sort_${NAME}.log" --export=NAME="${NAME}",DIR_DATA="${DIR_DATA}" --open-mode append ${DIR_LIB}/sorting.sh


### merge files
#　<NAME>_list.txtに記載していあるファイルを読みこんで結合する
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "sort_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=merge_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_sort_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; cat \`cat ${NAME}_list.txt\` > ${NAME}_sort.map"


### remove temporary files
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "merge_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=remove_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_sort_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; rm \`cat ${NAME}_list.txt\` ; rm ${NAME}_list.txt"


#-----------------------------------------------
# map fileをdatabaseに登録
#-----------------------------------------------
### register to database
# <NAME>_sort.mapを読みこんで、まったく同じエントリーは除去して、データベースに登録する
sbatch -N 1 -n 1 --exclusive=user --job-name=db_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_db_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Map2database.pl -i ${NAME}_sort.map -o ${NAME}.db"


#-----------------------------------------------
# HiC mapの解像度
#-----------------------------------------------
sbatch -n 1 --job-name=resolution_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_HiCmap_resolution_${NAME}.log" --open-mode append --wrap="bash ${DIR_LIB}/../../Statistics/HiCmap_resolution.sh -L $CHROM_LENGTH --count ${DIR_DATA}/${NAME}_count_for_resolution.txt ${DIR_DATA}/${NAME}_sort.map"

#-----------------------------------------------
# read数を調べる
#-----------------------------------------------
### count reads
# 元々のread数、ユニークなread数、quality > 10のread数、inter, intra-chromosome (including inter-arm), inter-armのread数をそれぞれデータベースから計算する
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "db_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=count_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_count_${NAME}.log" --export=SAMPLE="${NAME}",DIR_DATA="${DIR_DATA}" --open-mode append ${DIR_LIB}/Count_reads.sh

#-----------------------------------------------
# <NAME>.dbから<NAME>_fragment.dbを作る
#-----------------------------------------------
### split database
# database から、ユニークで、mapQ > 10 (defaultの閾値)なreadを抽出する
# 制限酵素部位が見つからないデータは削除する
# 制限酵素部位と、alignした部位を比較して、+なのに制限酵素部位の右にある場合は次の制限酵素に移動しなおし、-なのに制限酵素の左にある場合は前の制限酵素に移動しなおす
# 制限酵素部位ではなく、制限酵素断片のIDに直す(場所がLの場合はIDを１つ減らす）
# 10kb以内の距離で、向きが異なるペアについては除去する
# 染色体を１００に分割し、左の断片の中心位置のbinを基準に異なるファイルに出力する。出力したファイルのリストを<NAME>_list.txtとして出力する
sbatch -n 1 --job-name=DBsp_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_create_fragmentdb_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Split_database.pl -i ${NAME}.db -l ${CHROM_LENGTH} $(sq --node) -o ${NAME}_list.txt -m ${MAPQ_THRESHOLD} -e ${FILE_enzyme_def}"


### count duplicated number
# <NAME>_list.txtを順に読みこんで同じエントリーのデータが何個重複しているかを数え、数のデータを最後のカラムに加えて出力する
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "DBsp_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch --array=1-200 --job-name=DBcou_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_create_fragmentdb_${NAME}.log" --export=NAME="${NAME}",DIR_DATA="${DIR_DATA}",DIR_LIB="${DIR_LIB}" --open-mode append ${DIR_LIB}/Editing_filterReadsCount_jobArray.sh


### merge files
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "DBcou_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=DBmer_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_create_fragmentdb_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; cat \`cat ${NAME}_list.txt\` > ${NAME}_dataForFragDb.txt"


### remove temporary files
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "DBmer_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=DBrem_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_create_fragmentdb_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; rm \`cat ${NAME}_list.txt\` ; rm ${NAME}_list.txt"


### register to fragment.db
sbatch -N 1 -n 1 --job-name=DBreg_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_create_fragmentdb_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Register_filteredReads.pl -i ${NAME}_dataForFragDb.txt -o ${NAME}_fragment.db"



#-----------------------------------------------
# Distance curveの計算
#-----------------------------------------------
### distance curve calculation
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "DBreg_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=distance_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_distanceCurve_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; [ ! -e ${NAME}_distance.txt ] &&  perl ${DIR_LIB}/Create_distanceNormalize_data.pl -i ${NAME}_fragment.db -l ${FILE_CHROME_LENGTH} -o ${NAME}_distance.txt;"



#-----------------------------------------------
# DNAの量を見積もる
#-----------------------------------------------
### 見積もったDNAの量を<NAME>_DNA_amount.bedとして出力する
# 10kb以下の距離で、向きが同じでないもの（self ligation + undigest + 本当のHi-C)から、向きが同じもの(本当のHi-C)を引いた値を計算する
# 1kbのbinで数を数えて、結果をbedファイルとして出力する
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "db_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=DNA_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_DNA_amount_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Count_DNA_amount.pl -i ${NAME}.db -o ${NAME}_DNA_amount.bed"



#==============================================================
# Fragmentごとのread数の合計を調べて出力する(fragmentの長さも同時に出力)
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "DBreg_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=fragmentpro_${NAME} ${DEPEND} $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_fragment_property_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Fragment_property.pl -i ${NAME}_fragment.db > ${NAME}_fragment.txt"



#==============================================================
# Fragmentのblack listを作成
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "fragmentpro_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=badFragment_${NAME} ${DEPEND} $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_define_bad_fragment_${NAME}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Define_bad_fragment_threshold.R -i ${DIR_DATA}/${NAME}_fragment.txt --png ${DIR_DATA}/${NAME}_fragment.png -o ${DIR_DATA}/${NAME}_bad_fragment.txt --name ${NAME}"



#==============================================================
# inter-chromosomeの相互作用を調べる
#==============================================================
if [ ! -e ${DIR_DATA}/${NAME} ]; then
	mkdir ${DIR_DATA}/${NAME}
fi
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "badFragment_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=Inter_${NAME} $DEPEND $(sq --node) -o "${DIR_LOG}/${TIME_STAMP}_InterChromosome_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Make_association_from_fragmentdb_interChromosome.pl -i ${NAME}_fragment.db -o ${NAME}/InterChromosome.matrix -b ${NAME}_bad_fragment.txt"

