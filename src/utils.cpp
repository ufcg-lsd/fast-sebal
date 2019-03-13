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
					exit(7);
			}
			break;
		default:
			cerr << "Unsupported operation!" << endl;
			exit(7);
	}
	return ret;
};

MTL::MTL(){
    this->number_sensor = 0;
    this->julian_day = 0;
    this->year = 0;
    this->sun_elevation = 0;
    this->rad_add_10 = 0;
    this->rad_mult_10 = 0;
};

MTL::MTL(string metadata_path){
    map<string, string> mtl;

    ifstream in(metadata_path);
    if(!in.is_open() || !in){
        cerr << "Open metadata problem!" << endl;
        exit(2);
    }

    string line;
    while(getline(in, line)){
        stringstream lineReader(line);
        string token;
        vector<string> nline;
        while(lineReader >> token)
            nline.push_back(token);
        
        mtl[nline[0]] = nline[2];
    }

    in.close();

    char julian_day[3];
    julian_day[0] = mtl["LANDSAT_SCENE_ID"][14];
    julian_day[1] = mtl["LANDSAT_SCENE_ID"][15];
    julian_day[2] = mtl["LANDSAT_SCENE_ID"][16];

    char year[4];
    year[0] = mtl["LANDSAT_SCENE_ID"][10];
    year[1] = mtl["LANDSAT_SCENE_ID"][11];
    year[2] = mtl["LANDSAT_SCENE_ID"][12];
    year[3] = mtl["LANDSAT_SCENE_ID"][13];

    int hours = atoi(mtl["SCENE_CENTER_TIME"].substr(1, 2).c_str());
    int minutes = atoi(mtl["SCENE_CENTER_TIME"].substr(4, 2).c_str());

    this->number_sensor = atoi(new char(mtl["LANDSAT_SCENE_ID"][3]));
    this->julian_day = atoi(julian_day);
    this->year = atoi(year);
    this->sun_elevation = atof(mtl["SUN_ELEVATION"].c_str());
    this->distance_earth_sun = atof(mtl["EARTH_SUN_DISTANCE"].c_str());
    this->image_hour = (hours + minutes / 60.0) * 100;

    if(this->number_sensor == 8){
        rad_mult_10 = atof(mtl["RADIANCE_MULT_BAND_10"].c_str());
        rad_add_10 = atof(mtl["RADIANCE_ADD_BAND_10"].c_str());
    }
};

Sensor::Sensor(int number_sensor, int year){
    string sensor_path = capture_parameter_path(number_sensor, year);
    load_parameter_values(sensor_path);
}

string Sensor::capture_parameter_path(int number_sensor, int year){
    switch(number_sensor){
        case 8:
            return "src/parametros/LC.data";
            break;
        case 7:
            return "src/parametros/ETM.data";
            break;
        case 5:
            if(year < 1992)
                return "src/parametros/TM1.data";
            else
                return "src/parametros/TM2.data";
            break;
        default:
            cerr << "Sensor problem" << endl;
            exit(6);
    }
}

void Sensor::load_parameter_values(string sensor_path){
    ifstream in(sensor_path);
    if(!in.is_open() || !in) {
        cerr << "Open sensor parameters problem" << endl;
        exit(1);
    }

    string line, token;
    for(int i = 1; i < 8; i++){
        getline(in, line);
        stringstream lineReader(line);
        vector<string> nline;
        while(lineReader >> token)
            nline.push_back(token);

        this->parameters[i][this->GRESCALE] = atof(nline[0].c_str());
        this->parameters[i][this->BRESCALE] = atof(nline[1].c_str());
        this->parameters[i][this->ESUN] = atof(nline[2].c_str());
        this->parameters[i][this->WB] = atof(nline[3].c_str());
    }

    in.close();
};

Station::Station(){
    this->temperature_image = 0;
};

Station::Station(string station_data_path, double image_hour){
    ifstream in(station_data_path);
    if(!in.is_open() || !in){
        cerr << "Open station data problem!" << endl;
        exit(2);
    }

    string line;
    while(getline(in, line)){
        istringstream lineReader(line);
        vector<string> nline;
        string token;
        while(getline(lineReader, token, ';'))
            nline.push_back(token);

        if(nline.size())
            this->info.push_back(nline);
    }

    in.close();

    if(this->info.size() < 1){
        cerr << "Station data empty!" << endl;
        exit(12);
    }

    double diff = fabs(atof(this->info[0][2].c_str()) - image_hour);
    this->temperature_image = atof(this->info[0][6].c_str());
    this->v6 = atof(this->info[0][5].c_str());
    this->v7_max = atof(this->info[0][6].c_str());
    this->v7_min = atof(this->info[0][6].c_str());
    this->latitude = atof(this->info[0][3].c_str());
    this->longitude = atof(this->info[0][4].c_str());

    for(int i = 1; i < this->info.size(); i++){

        v7_max = max(v7_max, atof(this->info[i][6].c_str()));
        v7_min = min(v7_min, atof(this->info[i][6].c_str()));

        if(fabs(atof(this->info[i][2].c_str()) - image_hour) < diff){
            diff = fabs(atof(this->info[i][2].c_str()) - image_hour);
            this->temperature_image = atof(this->info[i][6].c_str());
            this->v6 = atof(this->info[i][5].c_str());
        }
    }
};

Candidate::Candidate(){
    this->ndvi = 0;
    this->temperature = 0;
    this->net_radiation = 0;
    this->soil_heat_flux = 0;
    this->ho = 0;
}

Candidate::Candidate(double ndvi, double temperature, double net_radiation, double soil_heat_flux, double ho){
    this->ndvi = ndvi;
    this->temperature = temperature;
    this->net_radiation = net_radiation;
    this->soil_heat_flux = soil_heat_flux;
    this->ho = ho;
}

void Candidate::setAerodynamicResistance(double u200, double A_ZOM, double B_ZOM, double VON_KARMAN){
    this->zom = exp(A_ZOM + B_ZOM * this->ndvi);
    this->ustar = (VON_KARMAN * u200)/log(200/this->zom);
    this->aerodynamic_resistance.push_back(log(2/0.1)/(this->ustar * VON_KARMAN));
}

bool compare_candidate_temperature(Candidate a, Candidate b){
    return a.temperature < b.temperature;
}

bool compare_candidate_ho(Candidate a, Candidate b){
    return a.ho < b.ho;
}

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

void write_line_tiff(TIFF* tif, double tif_line[], int line){
    if (TIFFWriteScanline(tif, tif_line, line) < 0){
        cerr << "Write problem!" << endl;
        exit(4);
    }
};

void close_tifs(TIFF* tifs[], int quant_tifs){
    for(int i = 1; i < quant_tifs; i++)
        TIFFClose(tifs[i]);
};