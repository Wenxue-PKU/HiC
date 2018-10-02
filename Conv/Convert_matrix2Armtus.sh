#!/bin/bash
# Armatusようにフォーマットを変換する

Usage="Usage : $0 [input.matrix] [out prefix] "

FILE_IN=$1
PREFIX=$2

if [ "$PREFIX" = "" ]; then
	echo $Usage
	exit
fi 

cat $FILE_IN | tail -n +2 | sed -e 's/:/\t/g' | awk -v OFS='\t' -v PRE=$PREFIX '{$3=$3+1; print > PRE"_"$1".matrix"}'

gzip ${PREFIX}_*.matrix

for CHR in I II III
do
	armatus -i ChIA-PET_Cut14_mix_20kb_Raw.matrix.gz -g .5 -o test -m
done

