#pragma once

#include "types.h"

const double VON_KARMAN = 0.41;
const double GRAVITY = 9.81;
const double RHO = 1.15; //Air density
const double SPECIFIC_HEAT_AIR = 1004;
const double GSC = 0.082; //Solar constant

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
	const int WIND_SPEED = 3;
	const double SURFACE_ROUGNESS = 0.024;
	const double A_ZOM = -3;
	const double B_ZOM = 6.47;
	const double INTERNALIZATION_FACTOR = 0.16;

	Station();
	Station(string estation_data_path, double image_hour);
};

struct Candidate{
	double ndvi, temperature, aerodynamic_resistance, ustar;
	double net_radiation, soil_heat_flux, ho, zom;

	Candidate();
	Candidate(double ndvi, double temperature, double net_radiation, double soil_heat_flux, double ho);
	void setAerodynamicResistance(double u200, double A_ZOM, double B_ZOM, double VON_KARMAN);
};

bool analisy_shadow(TIFF* read_bands[], TIFF* write_bands[], int number_sensor);

int set_mask(int number_sensor);

void setup(TIFF* new_tif, TIFF* base_tif);

void check_open_tiff(TIFF* tif);

void read_line_tiff(TIFF* tif, double tif_line[], int line);

void read_line_tiff(TIFF* tif, tdata_t tif_line, int line);

void write_line_tiff(TIFF* tif, double tif_line[], int line);

void close_tifs(TIFF* tifs[], int quant_tifs);