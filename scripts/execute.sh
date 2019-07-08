#!/bin/bash

make clean
make
./run input/_B2_converted.tif input/_B3_converted.tif input/_B4_converted.tif input/_B5_converted.tif input/_B6_converted.tif input/_B7_converted.tif input/_B10_converted.tif input/MTL.txt tal_converted.tif input/station.csv Testes/TesteClean/ -dist=0.98330 > Testes/TesteClean/outC &
PID=$$
sh collect-cpu-usage.sh $PID > /home/ubuntu/TDir/cpu.csv &
sh collect-memory-usage.sh $PID > /home/ubuntu/TDir/mem.csv &
sh collect-disk-usage.sh $PID > /home/ubuntu/TDir/disk.csv
