#!/bin/bash
N=$1
DS_PATH=$2
DIR_SEP=$3
LS_PATH=$4
K=1
N=`expr $N - $K`
for i in $(seq 0 $N);
do
	python find_reactions.py --ds-path ${DS_PATH} --sep-file ${DIR_SEP}sep_cmpds_${i}.pkl --ls-path ${LS_PATH} --ir-path reactions_${i}.txt --output reactions_list_${i}.pkl &
done;

