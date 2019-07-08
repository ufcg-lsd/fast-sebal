#!/bin/bash

make clean
make
./run input/_B2_converted.tif input/_B3_converted.tif input/_B4_converted.tif input/_B5_converted.tif input/_B6_converted.tif input/_B7_converted.tif input/_B10_converted.tif input/MTL.txt tal_converted.tif input/station.csv Testes/TesteClean/ -dist=0.98330 > Testes/TesteClean/outC &
PID=$(pidof ./run)
ps -aux | grep $PID
sh scripts/collect-cpu-usage.sh $PID > cpu.csv &
sh scripts/collect-memory-usage.sh $PID > mem.csv &
sh scripts/collect-disk-usage.sh $PID > disk.csv
