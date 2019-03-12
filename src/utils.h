#pragma once

#include "types.h"

struct PixelReader{
	uint16 sampleFormat;
	uint8 byteSize;
	tdata_t buffer;

	PixelReader();
	PixelReader(uint16 _sampleFormat, uint8 _byteSize,tdata_t _buffer);

	double read_pixel(uint32 colunm);
};

struct MTL{
    int number_sensor, julian_day, year;
    double sun_elevation, distance_earth_sun;
	double rad_mult_10, rad_add_10;
	double image_hour;

    MTL();
    MTL(string metadata_path);
};

struct Sensor{
	double parameters[8][4];
	const int GRESCALE = 0;
	const int BRESCALE = 1;
	const int ESUN = 2;
	const int WB = 3;

	Sensor();
	Sensor(int number_sensor, int year);
	string capture_parameter_path(int number_sensor, int year);
	void load_parameter_values(string sensor_path);
};

struct Station {
	vector< vector<string> > info;
	double temperature_image;

	Station();
	Station(string estation_data_path, double image_hour);
};

struct Candidate{

	Candidate();
};

bool analisy_shadow(TIFF* read_bands[], TIFF* write_bands[], int number_sensor);

int set_mask(int number_sensor);

void setup(TIFF* new_tif, TIFF* base_tif);

void check_open_tiff(TIFF* tif);

void read_line_tiff(TIFF* tif, double tif_line[], int line);

void read_line_tiff(TIFF* tif, tdata_t tif_line, int line);

void write_line_tiff(TIFF* tif, double tif_line[], int line);

void close_tifs(TIFF* tifs[], int quant_tifs);