#!/bin/bash

INPUT_PATH=$1
OUTPUT_PATH=$2
METHOD=$3

make clean
make
./run $INPUT_PATH/_B2_converted.tif $INPUT_PATH/_B3_converted.tif $INPUT_PATH/_B4_converted.tif $INPUT_PATH/_B5_converted.tif $INPUT_PATH/_B6_converted.tif $INPUT_PATH/_B7_converted.tif $INPUT_PATH/_B10_converted.tif $INPUT_PATH/MTL.txt $INPUT_PATH/tal_converted.tif $INPUT_PATH/station.csv $OUTPUT_PATH $INPUT_PATH/land_cover_final.tif -meth=$METHOD -nan=-3.39999995214436425e+38 > $OUTPUT_PATH/out.csv &
PID=$(pidof ./run)
ps -aux | grep $PID
sh scripts/collect-cpu-usage.sh $PID > $OUTPUT_PATH/cpu.csv &
sh scripts/collect-memory-usage.sh $PID > $OUTPUT_PATH/mem.csv &
sh scripts/collect-disk-usage.sh $PID > $OUTPUT_PATH/disk.csv
