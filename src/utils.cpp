#include "utils.h"

bool analisy_shadow(TIFF* read_bands[], TIFF* write_bands[], int number_sensor){
    int mask = set_mask(number_sensor);
    uint32 heigth_band, width_band;
    long long quant_pixels_valid = 0, quant_pixel_zero = 0;

    PixelReader pixel_read_bands[9];
    tdata_t line_bands[9];
    double *line_write_bands[9];

    for(int i = 1; i < 9; i++){
        uint16 sample_band;
        
        TIFFGetField(read_bands[i], TIFFTAG_SAMPLEFORMAT, &sample_band);
        TIFFGetField(read_bands[i], TIFFTAG_IMAGEWIDTH, &width_band);

        unsigned short byte_size_band = TIFFScanlineSize(read_bands[i]) / width_band;

        line_bands[i] = _TIFFmalloc(TIFFScanlineSize(read_bands[i]));
        line_write_bands[i] = new double[width_band];

        pixel_read_bands[i] = PixelReader(sample_band, byte_size_band, line_bands[i]);
    }

    TIFFGetField(read_bands[1], TIFFTAG_IMAGEWIDTH, &width_band);
    TIFFGetField(read_bands[1], TIFFTAG_IMAGELENGTH, &heigth_band);

    for(int line = 0; line < heigth_band; line ++){
        for(int i = 1; i < 9; i++)
            read_line_tiff(read_bands[i], line_bands[i], line);

        for(int col = 0; col < width_band; col ++){
            for(int i = 1; i < 9; i++)
                line_write_bands[i][col] = pixel_read_bands[i].read_pixel(col);

            if(fabs(line_write_bands[8][col] - 0) <= EPS){
                quant_pixel_zero ++;
                for(int i = 1; i < 9; i++)
                    line_write_bands[i][col] = NaN;
            }
            else if(fabs(line_write_bands[8][col] - mask) <= EPS || fabs(line_write_bands[8][col] - 20480) <= EPS)
                quant_pixels_valid ++;
            else{
                for(int i = 1; i < 9; i++)
                    line_write_bands[i][col] = NaN;
            }
        }

        for(int i = 1; i < 9; i++)
            write_line_tiff(write_bands[i], line_write_bands[i], line);
        
    }

    return (((double)quant_pixels_valid)/(heigth_band*width_band - quant_pixel_zero)) >= 0.01;
};

int set_mask(int number_sensor){
    if(number_sensor != 8)
        return 672;
    else
        return 2720;
};

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

void check_open_tiff(TIFF* tif){
    if(!tif){
        cerr << "Open tiff problem" << endl;
        exit(1);
    }
};

void read_line_tiff(TIFF* tif, double tif_line[], int line){
    if(TIFFReadScanline(tif, tif_line, line) < 0){
        cerr << "Read problem" << endl;
        exit(3);
    }
};

void read_line_tiff(TIFF* tif, tdata_t tif_line, int line){
    if(TIFFReadScanline(tif, tif_line, line) < 0){
        cerr << "Read problem" << endl;
        exit(3);
    }
};

double read_position_tiff(TIFF* tif, int col, int line){
    uint32 width_band;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width_band);

    double tif_line[width_band];

    read_line_tiff(tif, tif_line, line);

    return tif_line[col];
};
/*
//Generate small tiff 30x30
void fill_tiff(TIFF** tif, double min, double max){
    
    double tif_line[30];

    for(int line = 0; line < 30; line++){

        for(int col = 0; col < 30; col ++){
           tif_line[col] = getRandomDouble(min, max);
        }

        if (TIFFWriteScanline(*tif, tif_line, line) < 0){
            cerr << "Write problem!" << endl;
            exit(4);
        }
    }
}

double getRandomDouble(double min, double max){

    double f = (double) rand() / RAND_MAX;
    return min + f * (max - min);

}
*/
void write_line_tiff(TIFF* tif, double tif_line[], int line){

    uint32 width_band;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width_band);

    truncateArray(tif_line, width_band, 4);

    if (TIFFWriteScanline(tif, tif_line, line) < 0){
        cerr << "Write problem!" << endl;
        exit(4);
    }

};

void close_tifs(TIFF* tifs[], int quant_tifs){
    for(int i = 1; i < quant_tifs; i++)
        TIFFClose(tifs[i]);
};

void truncateArray(double* array, int size, int dec) {

    int pot = pow(10, dec);
    int arrayAux[size];

    for (int i = 0; i < size; i++){
        
        arrayAux[i] = int(floor(array[i] * pot));
        array[i] = double(arrayAux[i]) / pot;

    }

}

/*
The following definitions are from The art of computer programming by Knuth
*/

bool approximatelyEqual(double a, double b){
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}

bool essentiallyEqual(double a, double b){
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * EPS);
}

bool definitelyGreaterThan(double a, double b){
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}

bool definitelyLessThan(double a, double b){
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * EPS);
}