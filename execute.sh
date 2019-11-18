#!/bin/bash

INPUT_PATH=$1
OUTPUT_PATH=$2

make clean
make
./run $INPUT_PATH/_B2_converted.tif $INPUT_PATH/_B3_converted.tif $INPUT_PATH/_B4_converted.tif $INPUT_PATH/_B5_converted.tif $INPUT_PATH/_B6_converted.tif $INPUT_PATH/_B7_converted.tif $INPUT_PATH/_B10_converted.tif $INPUT_PATH/MTL.txt $INPUT_PATH/tal_converted.tif $INPUT_PATH/station.csv $OUTPUT_PATH -dist=0.98330 > $OUTPUT_PATH/out.txt &
PID=$(pidof ./run)
ps -aux | grep $PID
sh scripts/collect-cpu-usage.sh $PID > $OUTPUT_PATH/cpu.csv &
sh scripts/collect-memory-usage.sh $PID > $OUTPUT_PATH/mem.csv &
sh scripts/collect-disk-usage.sh $PID > $OUTPUT_PATH/disk.csv
