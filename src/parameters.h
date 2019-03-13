#pragma once

#include "types.h"

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
	double v6, v7_max, v7_min;
	double latitude, longitude;

	const int WIND_SPEED = 3;
	const double SURFACE_ROUGHNESS = 0.024;
	const double A_ZOM = -3;
	const double B_ZOM = 6.47;
	const double INTERNALIZATION_FACTOR = 0.16;

	Station();
	Station(string estation_data_path, double image_hour);
};