#pragma once

#include "types.h"
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

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

    MTL();
    MTL(string metadata_path);
};

bool analisy_shadow(TIFF* read_bands[], TIFF* write_bands[], int number_sensor);

int set_mask(int number_sensor);