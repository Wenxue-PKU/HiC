#!/bin/bash
# Extract significant associations

get_usage(){
	cat <<EOF

Usage : $0 [OPTION]

Description
	-h, --help
		show help

	-v, --version
		show version
	
	-n, --name
		sample name

	-r, --resolution [resolution (ex. 10kb)]
		resolution

	--max [max distance (ex. 2Mb)]
		default is 2Mb
	
	-d, --data [directory]
		data directory
	
	-o, --out [output file]
		output file
	
	-x, --organism [human|mouse]
		organism name
EOF

}

get_version(){
	echo "${0} version 1.0"
}

SHORT=hvn:r:d:o:x:
LONG=help,version,name:,resolution:,max:,data:,out:,organism:
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
		-r|--resolution)
			RESOLUTION="$2"
			shift 2
			;;
		--max)
			MAX_distance="$2"
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
[ ! -n "${RESOLUTION}" ] && echo "Please specify resolution" && exit 1
[ ! -n "${DIR_DATA}" ] && echo "Please specify data directory" && exit 1
[ ! -n "${FILE_OUT}" ] && echo "Please specify output file" && exit 1
[ ! -n "${ORGANISM}" ] && echo "Please specify organism" && exit 1
MAX_distance=${MAX_distance:-"2Mb"}
MAX_distance=${MAX_distance/Mb/000000}
MAX_distance=${MAX_distance/kb/000}

case $ORGANISM in
	human) CHRs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX" ;;
	mouse) CHRs="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX" ;;
	*)     echo "Please specify correct organism"
	       eixt 1 ;;
esac

DIR_tmp=${FILE_OUT}_tmpDir
[ ! -e ${DIR_tmp} ] && mkdir ${DIR_tmp} && mkdir ${DIR_tmp}/log ${DIR_tmp}/scores ${DIR_tmp}/img
UNIQ_ID=$(echo $FILE_OUT | rev | cut -c 1-12 | rev)
FILE_excel=${FILE_OUT/.txt/.xlsx}
FILE_log=${FILE_OUT/.txt/.log}

#==============================================================
# Significant fragment pairsを定義
#==============================================================
for CHR in $CHRs
do
	FILE_in=${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/${CHR}.rds
	sbatch -n 4 --job-name=si_${UNIQ_ID}_${CHR} -o "${DIR_tmp}/log/define_significant_pairs_${CHR}.log" --open-mode append --wrap="Rscript2 --vanilla --slave ${DIR_LIB}/Define_all_significant_pairs.R -i ${FILE_in} -o ${DIR_tmp}/scores/${CHR}.txt --max ${MAX_distance}"
done

### 結果をまとめてP-valueでソートする
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "si_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
sbatch -n 1 --job-name=si2_${UNIQ_ID} $DEPEND -o "${FILE_log}" --open-mode append <<-EOF
#!/bin/sh
cd ${DIR_tmp}/scores
cat chr1.txt | head -n1 > $FILE_OUT
ls chr*.txt | xargs -n1 | xargs -n1 -I@ sh -c "cat \@ | tail -n+2" | sort -k11,11g >> $FILE_OUT
cat ${DIR_tmp}/log/define_significant_pairs_*.log | awk -v OFS='\t' '{Nsig+=\$1; Nall++\$2}END{print Nsig,Nall}'
rm chr*.txt
EOF


#==============================================================
# TOP 50についてHi-C mapを描画する
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "si2_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
DB_loc=${DIR_tmp}/location.db
FILE_location=${DIR_tmp}/top1000.txt
sbatch -n 1 --job-name=drw_${UNIQ_ID} $DEPEND -o "${DIR_tmp}/log/drawGraph.log" --open-mode append <<-EOF
#!/bin/sh
echo "chr1 start1 end1 chr2 start2 end2" | tr ' ' '\t' > "$FILE_location"
cat $FILE_OUT | head -n 1001 | awk -v OFS='\t' 'NR>1{print \$1,\$2-200000,\$3+200000,\$4,\$5-200000,\$6+200000}' >> $FILE_location
file2database.R -i ${FILE_location} --id TRUE --db ${DB_loc} --table loc

for id in \$(sqlite3 ${DB_loc} "select id from loc")
do
	CHR1=\$(sqlite3 ${DB_loc} "select chr1 from loc where id='\${id}'")
	START1=\$(sqlite3 ${DB_loc} "select start1 from loc where id='\${id}'")
	END1=\$(sqlite3 ${DB_loc} "select end1 from loc where id='\${id}'")
	CHR2=\$(sqlite3 ${DB_loc} "select chr2 from loc where id='\${id}'")
	START2=\$(sqlite3 ${DB_loc} "select start2 from loc where id='\${id}'")
	END2=\$(sqlite3 ${DB_loc} "select end2 from loc where id='\${id}'")
	sbatch -n 4 --job-name=drw_${UNIQ_ID}_\${id} -o "${DIR_tmp}/log/drawGraph_\${id}.log" --open-mode append --wrap="Rscript --vanilla --slave ${DIR_LIB}/../../Draw/Draw_matrix.R -i ${DIR_DATA}/${NAME}/${RESOLUTION}/ICE/\${CHR1}.rds --normalize NA --zero NA --na na --chr \${CHR1} --start \${START1} --end \${END1} --chr2 \${CHR2} --start2 \${START2} --end2 \${END2} --unit p --max 0.95 --color red --width 500 -o ${DIR_tmp}/img/rank_\${id}.png"
done
EOF

#==============================================================
# 取得してきたスコアをまとめる
#==============================================================
JOB_ID=($(squeue -o "%j %F" -u htanizawa | grep -e "drw_${UNIQ_ID}" | cut -f2 -d' ' | xargs))
JOB_ID_string=$(IFS=:; echo "${JOB_ID[*]}")
DEPEND=""; [ -n "$JOB_ID_string" ] && DEPEND="--dependency=afterok:${JOB_ID_string}"
FILE_excel_data=${DIR_tmp}/data_for_excel.txt
sbatch -n 1 --job-name=xls_${UNIQ_ID} $DEPEND -o "${FILE_log}" --open-mode append <<-EOF
#!/bin/sh
echo -ne "id\t" > $FILE_excel_data
cat $FILE_OUT | head -n1 >> $FILE_excel_data
cat $FILE_OUT | head -n 1001 | awk -v OFS='\t' 'NR>1{rank=NR-1; print rank,\$0}' >> $FILE_excel_data
cd ${DIR_tmp}/img
python ${DIR_LIB}/Summarize_significant_pairs.py -i $FILE_excel_data -o ${FILE_excel} -t "${NAME}" --image ${DIR_tmp}/img
# rm -r ${DIR_tmp}
EOF
