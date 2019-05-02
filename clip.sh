#!/bin/bash

DIR=$1
MASK_PATH=$2
END="_masked.tif"

for i in $(ls $DIR | grep tif); do
    OUTPUT_FILE="$(echo $i | cut -d'.' -f1)$END"
    gdalwarp -cutline $MASK_PATH $i $OUTPUT_FILE
    rm $i
    mv $OUTPUT_FILE $i;
done