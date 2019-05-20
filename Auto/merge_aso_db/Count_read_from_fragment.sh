#!/bin/bash

cd ${DIR_DATA}

echo -n "intra-chromosome: "
sqlite3 ${SAMPLE}_fragment.db "select sum(score) from fragment where chr1 = chr2;"

echo -n "inter-chromosome: "
sqlite3 ${SAMPLE}_fragment.db "select sum(score) from fragment where chr1 != chr2;"

