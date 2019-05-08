#include "tiffio.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>

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
};

double diff_tifs(string compare_path_tiff1, string compare_path_tiff2, string diff_tiff_path){

    TIFF *tif1 = TIFFOpen(compare_path_tiff1.c_str(), "rm");
    TIFF *tif2 = TIFFOpen(compare_path_tiff2.c_str(), "rm");
    
    TIFF *diff = TIFFOpen(diff_tiff_path.c_str(), "w8m");
    setup(diff, tif1);

    uint32 heigth_band, width_band;
    TIFFGetField(tif1, TIFFTAG_IMAGEWIDTH, &width_band);
    TIFFGetField(tif1, TIFFTAG_IMAGELENGTH, &heigth_band);

    uint32 heigth_band_2, width_band_2;
    TIFFGetField(tif2, TIFFTAG_IMAGEWIDTH, &width_band_2);
    TIFFGetField(tif2, TIFFTAG_IMAGELENGTH, &heigth_band_2);

    if(heigth_band != heigth_band_2 || width_band != width_band_2){
        cerr << "Diference dimension tif!" << endl;
        exit(2);
    }


    double tif1_line[width_band];
    double tif2_line[width_band];
    double diff_line[width_band];

    double max_diff_relative = 0;
    double relative_error_diff;
    for(int line = 0; line < heigth_band; line++){

        if(TIFFReadScanline(tif1, tif1_line, line) < 0){
            cerr << "Read problem" << endl;
            exit(3);
        }

        if(TIFFReadScanline(tif2, tif2_line, line) < 0){
            cerr << "Read problem" << endl;
            exit(3);
        }

        for(int col = 0; col < width_band; col++){

            if (fabs(tif1_line[col]) < 1e-7 && fabs(tif2_line[col]) < 1e-7) {
                diff_line[col] = 0;
                continue;
            }
            
            diff_line[col] = fabs(tif1_line[col] - tif2_line[col]);
            relative_error_diff = fabs(diff_line[col]/tif2_line[col]) * 100;
            if(relative_error_diff > max_diff_relative){
                printf("TIF C %.7lf\n", tif1_line[col]);
                printf("TIF R %.7lf\n", tif2_line[col]);
                cout << diff_line[col] << endl;
                cout << relative_error_diff << endl;
                max_diff_relative = relative_error_diff;
            }
        }

        if (TIFFWriteScanline(diff, diff_line, line) < 0){
            cerr << "Write problem!" << endl;
            exit(4);
        }
    }

    TIFFClose(tif1);
    TIFFClose(tif2);
    TIFFClose(diff);

    return max_diff_relative;
};

int main(int argc, char *argv[]){
    string compare_path_tiff1 = argv[1];
    string compare_path_tiff2 = argv[2];
    string diff_tiff_path = argv[3];

    double max_diff = diff_tifs(compare_path_tiff1, compare_path_tiff2, diff_tiff_path);
    printf("Max diff value: %.10lf\n", max_diff);

    return 0;
}
