#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

void quartile(TIFF** target, double* vQuartile, int height_band, int width_band);

void hoFunction(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]);

Candidate getHotPixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);

Candidate getColdPixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);