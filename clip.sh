#!/bin/bash

DIR=$1
MASK_PATH=$2
END="_masked.nc"
ENDT="_masked.tif"

for i in $(ls $DIR | grep nc); do
    FILE_NAME_FULL=$DIR$i
    FILE_NAME=$(echo $i | cut -d'.' -f1)
    OUTPUT_FILE_NC=$DIR$FILE_NAME$END
    gdalwarp -cutline $MASK_PATH $FILE_NAME_FULL $OUTPUT_FILE_NC
    OUTPUT_FILE_TIF=$DIR$FILE_NAME$ENDT
    gdal_translate -of GTiff $OUTPUT_FILE_NC $OUTPUT_FILE_TIF
    rm $OUTPUT_FILE_NC
    OUTPUT_FILE_TIF_CONVERTED=$DIR$FILE_NAME"_converted.tif"
    ./test/conv $OUTPUT_FILE_TIF $OUTPUT_FILE_TIF_CONVERTED
    rm $OUTPUT_FILE_TIF
done