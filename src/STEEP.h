#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

void compute_H0(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]);

void filter_valid_values(double *target_line, double *target_values, int width_band, int *pos);

void get_quartiles(TIFF *target, double *v_quartile, int height_band, int width_band, double first_interval, double last_interval);

Candidate getHotPixelSTEPP(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);

Candidate getColdPixelSTEPP(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);
