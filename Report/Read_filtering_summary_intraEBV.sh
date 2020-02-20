#!/bin/bash
# intra EBV associations

cd ${DIR_DATA}

echo -n "no PCR duplicate: "
sqlite3 ${SAMPLE}.db "select count(id) from map where chr1 == 'EBV' and chr2 == 'EBV';"
echo -n "no repeat: "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and chr1 == 'EBV' and chr2 == 'EBV';"
echo -n "enough quality(MapQ>30): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30  and chr1 == 'EBV' and chr2 == 'EBV';"

echo -n ">10kb: "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1 == 'EBV' and chr2 == 'EBV' and abs(position1 - position2) > 10000;"

echo -n "<10kb(same direction): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1 == 'EBV' and chr2 == 'EBV' and abs(position1 - position2) < 10000 and direction1 = direction2;"
echo -n "<10kb(+ +): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1 == 'EBV' and chr2 == 'EBV' and abs(position1 - position2) < 10000 and direction1 = '+' and direction2 = '+';"
echo -n "<10kb(- -): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1 == 'EBV' and chr2 == 'EBV' and abs(position1 - position2) < 10000 and direction1 = '-' and direction2 = '-';"

echo -n "<10kb(different direction): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1 == 'EBV' and chr2 == 'EBV' and abs(position1 - position2) < 10000 and direction1 != direction2;"
echo -n "<10kb(+ -): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1 == 'EBV' and chr2 == 'EBV' and abs(position1 - position2) < 10000 and direction1 = '+' and direction2 = '-';"

echo -n "<10kb(- +): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1 == 'EBV' and chr2 == 'EBV' and abs(position1 - position2) < 10000 and direction1 = '-' and direction2 = '+';"