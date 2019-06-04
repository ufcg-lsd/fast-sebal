#pragma once

#include "types.h"
#include "pixel_reader.h"

bool analisy_shadow(TIFF* read_bands[], TIFF* write_bands[], int number_sensor);

int set_mask(int number_sensor);

void setup(TIFF* new_tif, TIFF* base_tif);

void check_open_tiff(TIFF* tif);

void read_line_tiff(TIFF* tif, double tif_line[], int line);

void read_line_tiff(TIFF* tif, tdata_t tif_line, int line);

double read_position_tiff(TIFF* tif, int col, int line);

void write_line_tiff(TIFF* tif, double tif_line[], int line);

void close_tifs(TIFF* tifs[], int quant_tifs);

void fill_tiff(TIFF** tif, double min, double max);

double getRandomDouble(double min, double max);

void truncateArray(double* array, int size, int dec);

/*
The following definitions are from The art of computer programming by Knuth
*/

bool approximatelyEqual(double a, double b);

bool essentiallyEqual(double a, double b);

bool definitelyGreaterThan(double a, double b);

bool definitelyLessThan(double a, double b);
