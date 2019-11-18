#!/bin/bash

for i in $(seq 1 30)
do
	echo Executing ${i}th experiment
	mkdir $OUTPUTPATH
	echo Output repository created and starting execution
	time ./execute.sh /dev/shm/input /dev/shm/output &&
	mv /dev/shm/output/*.csv /home/itallo/ExperimentoRAM/experimento${i}
        mv /dev/shm/output/*.txt /home/itallo/ExperimentoRAM/experimento${i}
        mv /dev/shm/output/ET24h.tif /home/itallo/ExperimentoRAM/experimento${i}
	rm /dev/shm/output/*
	echo Completing ${i}th experiment
done
cd ..
