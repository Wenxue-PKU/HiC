#!/bin/bash

cd ${DIR_DATA}

echo -n "no PCR duplicate: "
sqlite3 ${SAMPLE}.db "select count(id) from map;"
echo -n "no repeat: "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U';"
echo -n "enough quality(MapQ>30): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 ;"
echo -n "inter-chromosome: "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1!=chr2;"
echo -n ">10kb: "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) > 10000;"

echo -n "<10kb(same direction): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 10000 and direction1 = direction2;"
echo -n "<10kb(+ +): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 10000 and direction1 = '+' and direction2 = '+';"
echo -n "<10kb(- -): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 10000 and direction1 = '-' and direction2 = '-';"

echo -n "<10kb(different direction): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 10000 and direction1 != direction2;"
echo -n "<10kb(+ -): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 10000 and direction1 = '+' and direction2 = '-';"

echo -n "<10kb(- +): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 10000 and direction1 = '-' and direction2 = '+';"











