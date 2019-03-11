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

struct Estacao {
	double latitude; //v4
	double longitude; //v5
	vector<double> hora; //n sei o que eh, mas parece ser hora v3
	vector<double> temperatura; //v7
	vector<double> v6; //n sei o que eh

	Estacao();
	Estacao(string dados_estacao_path);
}

bool analisy_shadow(TIFF* read_bands[], TIFF* write_bands[], int number_sensor);

int set_mask(int number_sensor);

void setup(TIFF* new_tif, TIFF* base_tif);

void close_tifs(TIFF* tifs[], int quant_tifs);