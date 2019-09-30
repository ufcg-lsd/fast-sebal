#!/bin/bash

make
for i in $(seq 1 30)
do
	echo Executing ${i}th experiment
	OUTPUTPATH=DataExecution/execution${i}
	mkdir $OUTPUTPATH
	echo Output repository created and starting execution
	bash ./execute.sh $OUTPUTPATH &&
	echo Execution finished
	echo 'Tarzing TIFs'
	tar -cvzf $OUTPUT_PATH/execution${i}.tar.gz $OUTPUT_PATH/*.tif
	rm $OUTPUT_PATH/*.tif
	rm DataExecution/*.tif
	echo Completing ${i}th experiment
done
cd ..
make clean
