#!/bin/bash
# 複数のデータを結合させる

get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [input name lists. Space separated]

Description
	-h, --help
		show help

	-v, --version
		show version

	-d, --directory [data directory]
		directory name of analysis file locate

	-x, --organism [human|pombe|mouse]
		organism name

	-n, --name [merged file name]
		marged file's name

	-o, --log [log directory]
		log file directory

EOF

}

get_version(){
	echo "${0} version 1.0"
}


SHORT=hvd:x:o:n:
LONG=help,version,directory:,organism:,log:,name:
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
        -x|--organism)
            ORGANISM="$2"
            shift 2
            ;;
        -n|--name)
            NAME="$2"
            shift 2
            ;;
        -o|--log)
            DIR_LOG="$2"
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

[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
[ ! -n "${DIR_LOG}" ] && echo "Please specify log directory" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${NAME}" ] && echo "Please specify output file name" && exit 1

INPUT_FILES=$@




#-----------------------------------------------
# 複数の<NAME>_fragment.txtを読みこんで、まとめて左側の座標ごとに100のファイルに分類する
#-----------------------------------------------
sbatch -n 1  --job-name=Mg_Sp_${NAME} -o "${DIR_LOG}/${TIME_STAMP}_split_dataForFragDB_${NAME}.log" --open-mode append <<EOF
#!/bin/sh
cd ${DIR_DATA};
for IN in $INPUT_FILES
do
	perl ${DIR_LIB}/Split_dataForFragDb.pl -i \${IN}_dataForFragDb.txt -x ${ORGANISM} -o ${NAME}
done

EOF


#-----------------------------------------------
# <NAME>_list.txtから重複を除く
#-----------------------------------------------
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "Mg_Sp_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=Mg_Dup_${NAME} ${DEPEND} -o "${DIR_LOG}/${TIME_STAMP}_remove_duplicated_list_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; cat ${NAME}_list.txt | sort | uniq > tmplist_${NAME}_list.txt && mv tmplist_${NAME}_list.txt ${NAME}_list.txt"


#-----------------------------------------------
# count duplicated number
#-----------------------------------------------
# <NAME>_list.txtを順に読みこんで同じエントリーのデータが何個重複しているかを数え、数のデータを最後のカラムに加えて出力する
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "Mg_Dup_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch --array=1-200 --job-name=Mg_Cou_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_counting_number_for_merge_${NAME}.log" --export=NAME="${NAME}",DIR_DATA="${DIR_DATA}",DIR_LIB="${DIR_LIB}" --open-mode append ${DIR_LIB}/Count_mixed_count_jobArray.sh



### merge files
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "Mg_Cou_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=Mg_Mix_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_merge_counted_file_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; cat \`cat ${NAME}_list.txt\` > ${NAME}_dataForFragDb.txt"


### remove temporary files
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "Mg_Mix_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=Mg_Del_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_merge_counted_file_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; rm \`cat ${NAME}_list.txt\` ; rm ${NAME}_list.txt"




# DIR_LIBのパスを上のディレクトリに変更
DIR_LIB=${DIR_LIB}/../Basic


### register to fragment.db
sbatch -n 7 --job-name=Mg_Reg_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_register_merged_data_${NAME}.log" --export SQLITE_TMPDIR="/tmp" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Register_filteredReads.pl -i ${NAME}_dataForFragDb.txt -o ${NAME}_fragment.db"




#-----------------------------------------------
# read数を調べる
#-----------------------------------------------
### count reads
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "Mg_Reg_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=count_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_count_${NAME}.log" --export=SAMPLE="${NAME}",DIR_DATA="${DIR_DATA}" --open-mode append ${DIR_LIB}/../merge_aso_db/merge_aso_db/Count_read_from_fragment.sh



#-----------------------------------------------
# Distance curveの計算
#-----------------------------------------------
### distance curve calculation
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "Mg_Reg_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=distance_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_distanceCurve_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; [ ! -e ${NAME}_distance.txt ] &&  perl ${DIR_LIB}/Create_distanceNormalize_data.pl -i ${NAME}_fragment.db -x ${ORGANISM} -o ${NAME}_distance.txt;"



#==============================================================
# Fragmentごとのread数の合計を調べて出力する(fragmentの長さも同時に出力)
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "Mg_Reg_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=fragmentpro_${NAME} ${DEPEND} -o "${DIR_LOG}/${TIME_STAMP}_fragment_property_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Fragment_property.pl -i ${NAME}_fragment.db > ${NAME}_fragment.txt"



#==============================================================
# Fragmentのblack listを作成
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "fragmentpro_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=badFragment_${NAME} ${DEPEND} -o "${DIR_LOG}/${TIME_STAMP}_define_bad_fragment_${NAME}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/Define_bad_fragment_threshold.R -i ${DIR_DATA}/${NAME}_fragment.txt --png ${DIR_DATA}/${NAME}_fragment.png -o ${DIR_DATA}/${NAME}_bad_fragment.txt --name ${NAME}"



#==============================================================
# inter-chromosomeの相互作用を調べる
#==============================================================
if [ ! -e ${DIR_DATA}/${NAME} ]; then
	mkdir ${DIR_DATA}/${NAME}
fi
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "badFragment_${NAME}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=Inter_${NAME} $DEPEND -o "${DIR_LOG}/${TIME_STAMP}_InterChromosome_${NAME}.log" --open-mode append --wrap="cd ${DIR_DATA}; perl ${DIR_LIB}/Make_association_from_fragmentdb_interChromosome.pl -i ${NAME}_fragment.db -o ${NAME}/InterChromosome.matrix -b ${NAME}_bad_fragment.txt"
