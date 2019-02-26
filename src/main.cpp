#include "utils.h"
#include <tiffio.h>
#include <string>
#include <iostream>

using namespace std;

void setup(TIFF* new_tif, TIFF* base_tif){
    uint32 image_width, image_length;

    TIFFGetField(base_tif, TIFFTAG_IMAGEWIDTH,      &image_width);
    TIFFGetField(base_tif, TIFFTAG_IMAGELENGTH,     &image_length);
    
    TIFFSetField(new_tif, TIFFTAG_IMAGEWIDTH     , image_width); 
    TIFFSetField(new_tif, TIFFTAG_IMAGELENGTH    , image_length);
    TIFFSetField(new_tif, TIFFTAG_BITSPERSAMPLE  , 64);
    TIFFSetField(new_tif, TIFFTAG_SAMPLEFORMAT   , 3);
    TIFFSetField(new_tif, TIFFTAG_COMPRESSION    , 1);
    TIFFSetField(new_tif, TIFFTAG_PHOTOMETRIC    , 1);
    TIFFSetField(new_tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(new_tif, TIFFTAG_ROWSPERSTRIP   , 1);
    TIFFSetField(new_tif, TIFFTAG_RESOLUTIONUNIT , 1);
    TIFFSetField(new_tif, TIFFTAG_XRESOLUTION    , 1);
    TIFFSetField(new_tif, TIFFTAG_YRESOLUTION    , 1);
    TIFFSetField(new_tif, TIFFTAG_PLANARCONFIG   , PLANARCONFIG_CONTIG);
}

int main(int argc, char *argv[]){
    string metadata_path = argv[9];
    MTL mtl = MTL(metadata_path);

    TIFF *read_bands[9], *write_bands[9];
    for(int i = 1; i < 9; i++){
        string path_tiff_base = argv[i];
        read_bands[i] = TIFFOpen(path_tiff_base.c_str(), "rm");

        string path_tiff_write = path_tiff_base.substr(0, path_tiff_base.size() - 4) + "_write.tif";
        write_bands[i] = TIFFOpen(path_tiff_write.c_str(), "w8m");

        if(!read_bands[i]){
            cerr << "Open problem" << endl;
            return 1;
        } 

        //Check write bands (ERROR)

        setup(write_bands[i], read_bands[i]);
    }

    if(!analisy_shadow(read_bands, write_bands, mtl.number_sensor)){
        for(int i = 1; i < 9; i++){
            TIFFClose(read_bands[i]);
            TIFFClose(write_bands[i]);
        }
        cerr << "Shadow problem" << endl;
        exit(5);
    }

    //calculate products (NDVI, EVI, ...)
    //calculate ET, ET24h

    for(int i = 1; i < 9; i++){
        TIFFClose(read_bands[i]);
        TIFFClose(write_bands[i]);
    }

    return 0;
}