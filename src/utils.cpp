#include "utils.h"

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
					exit(1);
			}
			break;
		default:
			cerr << "Unsupported operation!" << endl;
			exit(1);
	}
	return ret;
};

MTL::MTL(){
    this->number_sensor = 0;
    this->julian_day = 0;
    this->year = 0;
    this->sun_elevation = 0;
};

MTL::MTL(string metadata_path){
    map<string, string> mtl;

    ifstream in(metadata_path);
    if(!in.is_open() || !in) exit(1);

    string line;
    while(getline(in, line)){
        stringstream lineReader(line);
        string token;
        vector<string> nline;
        while(lineReader >> token)
            nline.push_back(token);
        
        mtl[nline[0]] = nline[2];
    }

    char julian_day[3];
    julian_day[0] = mtl["LANDSAT_SCENE_ID"][14];
    julian_day[1] = mtl["LANDSAT_SCENE_ID"][15];
    julian_day[2] = mtl["LANDSAT_SCENE_ID"][16];

    char year[4];
    year[0] = mtl["LANDSAT_SCENE_ID"][10];
    year[1] = mtl["LANDSAT_SCENE_ID"][11];
    year[2] = mtl["LANDSAT_SCENE_ID"][12];
    year[3] = mtl["LANDSAT_SCENE_ID"][13];

    this->number_sensor = atoi(new char(mtl["LANDSAT_SCENE_ID"][3]));
    this->julian_day = atoi(julian_day);
    this->year = atoi(year);
    this->sun_elevation = atof(mtl["SUN_ELEVATION"].c_str());
    this->distance_earth_sun = atof(mtl["EARTH_SUN_DISTANCE"].c_str());
};

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

    TIFFGetField(read_bands[8], TIFFTAG_IMAGEWIDTH, &width_band);
    TIFFGetField(read_bands[8], TIFFTAG_IMAGELENGTH, &heigth_band);

    for(int line = 0; line < heigth_band; line ++){
        for(int i = 1; i < 9; i++)
            if(TIFFReadScanline(read_bands[i], line_bands[i], line) < 0){
                cerr << "Read problem" << endl;
                exit(3);
            }

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
            if(TIFFWriteScanline(write_bands[i], line_write_bands[i], line) < 0){
                cerr << "Write problem" << endl;
                exit(4);
            }
        
    }

    return (((double)quant_pixels_valid)/(heigth_band*width_band - quant_pixel_zero)) >= 0.01;
};

int set_mask(int number_sensor){
    if(number_sensor != 8)
        return 672;
    else
        return 2720;
};