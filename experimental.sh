#!/bin/bash

for i in $(seq 1 30)
do
	echo Executing ${i}th experiment
	OUTPUTPATH=DataExecution/execution${i}
	mkdir $OUTPUTPATH
	echo Output repository created and starting execution
	time ./execute.sh $OUTPUTPATH &&
	echo Completing ${i}th experiment
done
cd ..
