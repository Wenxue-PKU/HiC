#!/bin/bash

cd ${DIR_DATA}

[ ! -e ${SAMPLE}_fragment.db ] && echo "${SAMPLE}_fragment.db was not found" && exit 1

echo -n "intra-chromosome: "
sqlite3 ${SAMPLE}_fragment.db "select sum(score) from fragment where chr1 = chr2;"

echo -n "inter-chromosome: "
sqlite3 ${SAMPLE}_fragment.db "select sum(score) from fragment where chr1 != chr2;"

