#pragma once

#include "types.h"
#include "utils.h"

void radiance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, int width_band, double radiance_line[][8]);
void reflectance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, double reflectance_line[][8]);
string tal_function(TIFF* raster_elevation, string output_path);
void albedo_function(double reflectance_line[][8], Sensor sensor, double tal_line[], int width_band, int number_sensor, int line, TIFF* albedo);
void short_wave_radiation_function();
void ndvi_function();
void savi_function();
void lai_function();
void evi_function();
void enb_emissivity_function();
void eo_emissivity_function();
void ea_emissivity_function();
void kelvin_surface_function();
void large_wave_radiation_surface_function();
void large_wave_radiation_atmosphere_function();
void radiation_net_function();
void soil_heat_flux_function();