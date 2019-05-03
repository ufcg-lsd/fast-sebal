#!/bin/bash

DIR=$1
MASK_PATH=$2
END="_masked.tif"

for i in $(ls $DIR | grep nc); do
    FILE_NAME_FULL=$DIR$i
    FILE_NAME=$(echo $i | cut -d'.' -f1)
    OUTPUT_FILE=$DIR$FILE_NAME$END
    gdalwarp -q -cutline $MASK_PATH -tr 0.000272 0.000271 -of GTiff $FILE_NAME_FULL $OUTPUT_FILE
    OUTPUT_FILE_CONVERTED=$DIR$FILE_NAME"_converted.tif"
    ./test/conv $OUTPUT_FILE $OUTPUT_FILE_CONVERTED
    #rm $OUTPUT_FILE_TIF
done
