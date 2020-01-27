#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

void testLandCoverHomogeneity(TIFF* landCover, TIFF* mask);

void testHomogeneity(TIFF* ndvi, TIFF* surface_temperature, TIFF* albedo, TIFF* mask);

void testMorphological(TIFF* input, TIFF* output, int groupSize);

Candidate esaHotPixel(TIFF** mask, TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);

Candidate esaColdPixel(TIFF** mask, TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);