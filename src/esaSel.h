#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

void testLandCoverHomogeneity(TIFF* landCover, TIFF* mask);

void testHomogeneity(TIFF* ndvi, TIFF* surface_temperature, TIFF* albedo, TIFF* maskLC, TIFF* output);

void testMorphological(TIFF* input, TIFF* output, int groupSize);

void hoCalc(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]);

pair<Candidate, Candidate> esaPixelSelect(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, TIFF** landCover, int height_band, int width_band, string output_path);