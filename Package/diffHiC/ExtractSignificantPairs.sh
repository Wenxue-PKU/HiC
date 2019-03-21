#!/bin/bash
# Significantなペアを抽出する(diffHiC packageを使う)


get_usage(){
	cat <<EOF

Usage : $0 [OPTION] [sample names (separated by space)]

Fold change will be set2 / set1 (set1 is initial condition. set2 is changed condition)

Description
	-h, --help
		show help

	-v, --version
		show version

	-r, --resolution [resolution]
		resolution

	-g, --group [group numbers (1|2)]
		group setting. separated by , ex (1,1,2,2)

	-d, --data [directory]
		data directory
	
	-o, --out [output file]
		output file
	
	-x, --organism [human|mouse]
		organism name

	-t, --title [title]
		title appeared in excel tab

	--draw_range [ex 200000]
		half size of example map width. resolution x 20 seems good. (for 10kb resolution, draw 200kb half size = 200000)
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvr:g:d:o:x:t:
LONG=help,version,resolution:,group:,data:,out:,organism:,title:,draw_range:
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
			RESOLUTION="$2"
			shift 2
			;;			
		-g|--group)
			group="$2"
			shift 2
			;;
		-d|--data)
			DIR_DATA="$2"
			shift 2
			;;
		-o|--out)
			FILE_OUT="$2"
			shift 2
			;;
		-x|--organism)
			ORGANISM="$2"
			shift 2
			;;
		-t|--title)
			EXCEL_tab="$2"
			shift 2
			;;
		--draw_range)
			DRAW_WIDTH="$2"
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

[ ! -n "${group}" ] && echo "Please specify group" && exit 1
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output file" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
[ ! -n "${DRAW_WIDTH}" ] && echo "Please specify draw half width" && exit 1
EXCEL_tab=${EXCEL_tab:-"Sheet1"}

NAME_LIST="$@"

case $ORGANISM in
	human) CHRs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX" ;;
	mouse) CHRs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX" ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac

DIR_tmp=${FILE_OUT}_tmpDir
[ ! -e ${DIR_tmp} ] && mkdir ${DIR_tmp} && mkdir ${DIR_tmp}/log ${DIR_tmp}/scores ${DIR_tmp}/hic ${DIR_tmp}/img
UNIQ_ID=$(echo $FILE_OUT | rev | cut -c 1-12 | rev)
FILE_excel=${FILE_OUT/.txt/.xlsx}
FILE_top100=${FILE_OUT/.txt/_top100.xlsx}
FILE_log=${FILE_OUT/.txt/.log}

# #==============================================================
# # Significant fragment pairsを定義
# #==============================================================
# for CHR in $CHRs
# do
# 	FILE_in=$(echo "$NAME_LIST" | xargs -n1 | xargs -I@ sh -c "echo ${DIR_DATA}/\@/${RESOLUTION}/Raw/${CHR}.rds" | xargs | tr ' ' ',')
# 	sbatch -n 4 --job-name=si_${UNIQ_ID}_${CHR} -o "${DIR_tmp}/log/define_significant_pairs_${CHR}.log" --open-mode append --wrap="Rscript2 --vanilla --slave ${DIR_LIB}/Extract_diff_pairs.R -i ${FILE_in} -o ${DIR_tmp}/scores/${CHR}.txt --group ${group}"
# done

# ### 結果をまとめてFDRでソートする
# JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "si_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
# JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
# DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
# sbatch -n 1 --job-name=si2_${UNIQ_ID} $DEPEND -o "${DIR_tmp}/log/mix_result.log" --open-mode append <<-EOF
# #!/bin/sh
# cd ${DIR_tmp}/scores
# echo "chr1 start1 end1 chr2 start2 end2 distance logFC logCPM Pvalue FDR" | tr ' ' '\t' > ../pvalue.txt
# ls chr*.txt | xargs -n1 | xargs -n1 -I@ sh -c "cat \@ | tail -n+2" | sort -k11,11g | awk -v OFS='\t' '{m1=int((\$2+\$3)/2); m2=int((\$5+\$6)/2); dis=m1-m2; print \$1,\$2,\$3,\$4,\$5,\$6,dis,\$7,\$8,\$10,\$11}' >> ../pvalue.txt
# cat ../pvalue.txt | awk -v OFS='\t' 'NR>1{print \$1,\$2,\$3,0,\$4,\$5,\$6,0}' > ../pairs.txt
# EOF


# #==============================================================
# # Significant pairsのスコアを出力する
# #==============================================================
# JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "si2_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
# JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
# DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
# for NAME in $NAME_LIST
# do
# 	sbatch -n 4 --job-name=gs_${UNIQ_ID}_${NAME} $DEPEND -o "${DIR_tmp}/log/get_hicScore.log" --open-mode append --wrap="cd ${DIR_tmp} && Rscript --vanilla --slave ${DIR_LIB}/../../Filter/Exctact_scores_of_indicated_pairs.R --in pairs.txt --dir ${DIR_DATA}/${NAME}/${RESOLUTION}/Raw --out ${DIR_tmp}/hic/${NAME}.txt"
# done

#==============================================================
# 取得してきたスコアをまとめる
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "gs_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=hic_${UNIQ_ID} $DEPEND -o "/dev/null" --open-mode append <<-EOF
#!/bin/sh
cd ${DIR_tmp}/hic
echo "$NAME_LIST" | tr ' ' '\t' > ../hic.txt
echo "$NAME_LIST" | xargs -n1 | xargs -I@ sh -c "echo \@.txt" | xargs | xargs paste >> ../hic.txt
cd ${DIR_tmp}
paste pvalue.txt hic.txt > $FILE_OUT
python ${DIR_LIB}/Summarize_significant_pairs.py -i ${FILE_OUT} -o ${FILE_excel} -g "${group}" -t "${EXCEL_tab}"
EOF


#==============================================================
# TOP 100についてHi-C mapを描画する (distance > 100kb)
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "si2_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
DB_loc=${DIR_tmp}/location.db
FILE_location=${DIR_tmp}/location.txt
sbatch -n 1 --job-name=drw_${UNIQ_ID} $DEPEND -o "${DIR_tmp}/log/drawGraph.log" --open-mode append <<-EOF
#!/bin/sh
echo "chr1 start1 end1 chr2 start2 end2" | tr ' ' '\t' > "$FILE_location"
cat ${DIR_tmp}/pvalue.txt | awk -v OFS='\t' -v halfW=${DRAW_WIDTH} 'NR>1&&\$7>100000{print \$1,\$2-halfW,\$3+halfW,\$4,\$5-halfW,\$6+halfW}' | head -n 101 >> $FILE_location
file2database.R -i ${FILE_location} --id TRUE --db ${DB_loc} --table loc
for NAME in $NAME_LIST
do
	for id in \$(seq 1 100)
	do
		CHR1=\$(sqlite3 ${DB_loc} "select chr1 from loc where id='\${id}'")
		START1=\$(sqlite3 ${DB_loc} "select start1 from loc where id='\${id}'")
		END1=\$(sqlite3 ${DB_loc} "select end1 from loc where id='\${id}'")
		CHR2=\$(sqlite3 ${DB_loc} "select chr2 from loc where id='\${id}'")
		START2=\$(sqlite3 ${DB_loc} "select start2 from loc where id='\${id}'")
		END2=\$(sqlite3 ${DB_loc} "select end2 from loc where id='\${id}'")
		sbatch -n 4 --job-name=drw_${UNIQ_ID}_\${id} -o "${DIR_tmp}/log/drawGraph_\${id}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/../../Draw/Draw_matrix.R -i ${DIR_DATA}/\${NAME}/${RESOLUTION}/Raw/\${CHR1}.rds --normalize NA --zero NA --na na --chr \${CHR1} --start \${START1} --end \${END1} --chr2 \${CHR2} --start2 \${START2} --end2 \${END2} --unit p --max 0.95 --color red --width 500 -o ${DIR_tmp}/img/rank_\${id}_\${NAME}.png"
	done
done
EOF

#==============================================================
# 取得してきたスコアと画像をまとめる
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "drw_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
FILE_excel_data=${DIR_tmp}/data_for_excel.txt
sbatch -n 1 --job-name=xls_${UNIQ_ID} $DEPEND -o "${FILE_log}" --open-mode append <<-EEE
#!/bin/sh
cat ${DIR_tmp}/pvalue.txt | head -n1 > $FILE_excel_data
cat ${DIR_tmp}/pvalue.txt | awk -v OFS='\t' 'NR>1&&\$7>100000' | head -n 101 >> $FILE_excel_data
cd ${DIR_tmp}/img
python ${DIR_LIB}/Summarize_top100.py -i $FILE_excel_data -o ${FILE_top100} --image ${DIR_tmp}/img --samples $(echo "$NAME_LIST" | tr ' ' ',')
find ${DIR_tmp}/log -type f -empty -delete
# rm -r ${DIR_tmp}
EEE