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
echo -n ">20kb: "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) > 20000;"

echo -n "<20kb(same direction): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 20000 and direction1 = direction2;"
echo -n "<20kb(+ +): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 20000 and direction1 = '+' and direction2 = '+';"
echo -n "<20kb(- -): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 20000 and direction1 = '-' and direction2 = '-';"

echo -n "<20kb(different direction): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 20000 and direction1 != direction2;"
echo -n "<20kb(+ -): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 20000 and direction1 = '+' and direction2 = '-';"

echo -n "<20kb(- +): "
sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and abs(position1 - position2) < 20000 and direction1 = '-' and direction2 = '+';"


# pombeのinter-arm
#echo -n "inter-arm: "
#sqlite3 ${SAMPLE}.db "select count(id) from map where uniq1='U' and uniq2 = 'U' and mapQ1 > 30 and mapQ2 > 30 and chr1=chr2 and ( ( chr1='I' and position1 < 3753687 and position2 > 3789421 ) or ( chr1='II' and position1 < 1602264 and position2 > 1644747 ) or ( chr1='III' and position1 < 1070904 and position2 > 1137003 ));"










