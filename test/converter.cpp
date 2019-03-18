#include "tiffio.h"
#include <string>
#include <math.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>

using namespace std;

const double EPS = 1e-7;
const double NaN = -sqrt(-1.0);

struct PixelReader{
	uint16 sampleFormat;
	uint8 byteSize;
	tdata_t buffer;

	PixelReader();
	PixelReader(uint16 _sampleFormat, uint8 _byteSize,tdata_t _buffer);

	double read_pixel(uint32 colunm);
};

PixelReader::PixelReader() {
	sampleFormat = 0;
	byteSize = 0;
	buffer = NULL;
};

PixelReader::PixelReader(uint16 _sampleFormat, uint8 _byteSize, tdata_t _buffer){
	sampleFormat = _sampleFormat;
	byteSize = _byteSize;
	buffer = _buffer;
};

double PixelReader::read_pixel(uint32 colunm){
	double ret = 0;
	switch(sampleFormat){
		case 1:
			{
				uint64 value = 0;
				memcpy(&value, buffer + (colunm * byteSize), byteSize);
				ret = value;
			}
			break;
		case 2:
			{
				int64 value = 0;
				memcpy(&value, buffer + (colunm * byteSize), byteSize);
				ret = value;
			}
			break;
		case 3:
			switch(byteSize){
				case 4:
					{
						float value = 0;
						memcpy(&value, buffer + (colunm * byteSize), byteSize);
						ret = value;
					}
					break;
				case 8:
					{
						double value = 0;
						memcpy(&value, buffer + (colunm * byteSize), byteSize);
						ret = value;
					}
					break;
				case 16:
					{
						long double value = 0;
						memcpy(&value, buffer + (colunm * byteSize), byteSize);
						ret = value;
					}
					break;
				default:
					cerr << "Unsupported operation!" << endl;
					exit(2);
			}
			break;
		default:
			cerr << "Unsupported operation!" << endl;
			exit(2);
	}
	return ret;
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

void converter_tif(string base_tif_path, string output_tif_path, double mask = 0.0){
    uint32 heigth_band, width_band;
    uint16 sample_band;

    TIFF *base = TIFFOpen(base_tif_path.c_str(), "rm");
    TIFFGetField(base, TIFFTAG_IMAGEWIDTH, &width_band);
    TIFFGetField(base, TIFFTAG_IMAGELENGTH, &heigth_band);
    TIFFGetField(base, TIFFTAG_SAMPLEFORMAT, &sample_band);

    TIFF *output = TIFFOpen(output_tif_path.c_str(), "w8m");
    setup(output, base);

    PixelReader pixel_read_band;
    tdata_t line_band;
    double output_line_band[width_band];

    unsigned short byte_size_band = TIFFScanlineSize(base) / width_band;
    line_band = _TIFFmalloc(TIFFScanlineSize(base));
    pixel_read_band = PixelReader(sample_band, byte_size_band, line_band);

    for(int line = 0; line < heigth_band; line ++){
        if(TIFFReadScanline(base, line_band, line) < 0){
            cerr << "Read problem" << endl;
            exit(3);
        }

        for(int col = 0; col < width_band; col ++){
            output_line_band[col] = pixel_read_band.read_pixel(col);
                
            if(fabs(output_line_band[col] - mask) < EPS)
                output_line_band[col] = NaN;
        }

        if (TIFFWriteScanline(output, output_line_band, line) < 0){
            cerr << "Write problem!" << endl;
            exit(4);
        }
    }

    TIFFClose(base);
    TIFFClose(output);
};

int main(int argc, char *argv[]){
    string base_tif_path = argv[1];
    string output_tif_path = argv[2];
    double mask = 0.0;

    if(argc == 4)
        mask = atof(argv[3]);

    converter_tif(base_tif_path, output_tif_path, mask);

    return 0;
}