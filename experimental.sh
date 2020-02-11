#!/bin/bash

FINAL_PATH=$1
METHOD=$2

for i in $(seq -f "%02g" 1 30)
do
	echo Executing ${i}th experiment
	mkdir $FINAL_PATH/experimento${i}
	time ./execute.sh /dev/shm/input /dev/shm/output $METHOD &&
	mv /dev/shm/output/*.csv $FINAL_PATH/experimento${i}
        mv /dev/shm/output/ET24h.tif $FINAL_PATH/experimento${i}
	rm /dev/shm/output/*
	echo Completing ${i}th experiment
done
