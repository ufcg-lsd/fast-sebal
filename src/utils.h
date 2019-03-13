#pragma once

#include "types.h"

bool analisy_shadow(TIFF* read_bands[], TIFF* write_bands[], int number_sensor);

int set_mask(int number_sensor);

void setup(TIFF* new_tif, TIFF* base_tif);

void check_open_tiff(TIFF* tif);

void read_line_tiff(TIFF* tif, double tif_line[], int line);

void read_line_tiff(TIFF* tif, tdata_t tif_line, int line);

void write_line_tiff(TIFF* tif, double tif_line[], int line);

void close_tifs(TIFF* tifs[], int quant_tifs);