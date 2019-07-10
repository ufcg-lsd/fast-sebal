#pragma once

#include "types.h"

/**
 * @brief  Struct to hold some metadata informations.
 */
struct MTL{

    int number_sensor, julian_day, year;
    double sun_elevation, distance_earth_sun;
	double rad_mult_10, rad_add_10;
	double image_hour;

    /**
     * @brief  Empty constructor. Setting all attributes to 0.
     */
    MTL();

    /**
     * @brief  Constructor receiving the metadata path.
     * @param  metadata_path: Metadata file path.
     */
    MTL(string metadata_path);

};

/**
 * @brief  Satellite sensor struct.
 */
struct Sensor{
	double parameters[8][4];
	const int GRESCALE = 0;
	const int BRESCALE = 1;
	const int ESUN = 2;
	const int WB = 3;

	/**
	 * @brief  Empty constructor.
	 */
	Sensor();

	/**
	 * @brief  Constructor.
	 * @param  number_sensor: Number of the satellite sensor.
	 * @param  year: Year of image.
	 */
	Sensor(int number_sensor, int year);

	/**
	 * @brief  Get the path of the correct sensor parameters based on his number.
	 * @param  number_sensor: Number of the satellite sensor.
	 * @param  year: Year of image. 
	 * @retval Path to sensor parameter file.
	 */
	string capture_parameter_path(int number_sensor, int year);

	/**
	 * @brief  Loads the sensor parameters.
	 * @param  sensor_path: Path to sensor parameter file.
	 */
	void load_parameter_values(string sensor_path);
};

/**
 * @brief  Struct to hold the weather station data.
 */
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

	/**
	 * @brief  Empty constructor. Set temperature_image to 0.
	 */
	Station();

	/**
	 * @brief  Constructor.
	 * @param  station_data_path: Weather station data file.
	 * @param  image_hour: Image hour.
	 */
	Station(string station_data_path, double image_hour);
};