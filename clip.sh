#!/bin/bash

DIR=$1
MASK_PATH=$2
END="_masked.nc"
ENDT="_masked.tif"

for i in $(ls $DIR | grep nc); do
    FILE_NAME=$(echo $i | cut -d'.' -f1)
    OUTPUT_FILE_NC="$FILE_NAME$END"
    gdalwarp -cutline $MASK_PATH $i $OUTPUT_FILE_NC
    OUTPUT_FILE_TIF="$FILE_NAME$ENDT"
    gdal_translate -of GTiff $OUTPUT_FILE_NC $OUTPUT_FILE_TIF
    rm $OUTPUT_FILE_NC
done