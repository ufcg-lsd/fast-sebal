#!/bin/bash

OUTPUT_PATH=$1

make clean
make
./run input/_B2_converted.tif input/_B3_converted.tif input/_B4_converted.tif input/_B5_converted.tif input/_B6_converted.tif input/_B7_converted.tif input/_B10_converted.tif input/MTL.txt input/tal_converted.tif input/station.csv $OUTPUT_PATH -dist=0.98330 > $OUTPUT_PATH/out.txt &
PID=$(pidof ./run)
ps -aux | grep $PID
sh scripts/collect-cpu-usage.sh $PID > $OUTPUT_PATH/cpu.csv &
sh scripts/collect-memory-usage.sh $PID > $OUTPUT_PATH/mem.csv &
sh scripts/collect-disk-usage.sh $PID > $OUTPUT_PATH/disk.csv
