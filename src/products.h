#pragma once

#include "types.h"
#include "utils.h"

string tal_function(TIFF* raster_elevation, string output_path);
void radiance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, int width_band, double radiance_line[][8]);
void reflectance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, double reflectance_line[][8]);
void albedo_function(double reflectance_line[][8], Sensor sensor, double tal_line[], int width_band, int number_sensor, int line, double albedo_line[]);
void ndvi_function(double reflectance_line[][8], int width_band, int line, TIFF *ndvi);
void lai_function(double reflectance_line[][8], int width_band, int line, TIFF *lai);
void evi_function();
void enb_emissivity_function();
void eo_emissivity_function();
void ea_emissivity_function();
void kelvin_surface_function();
void short_wave_radiation_function(double tal_line[], MTL mtl, int width_band, double short_wave_radiation_line[]);
void large_wave_radiation_surface_function();
void large_wave_radiation_atmosphere_function();
void net_radiation_function(double tal_line[], double albedo_line[], MTL mtl, int width_band, int line, TIFF *net_radiation);
void soil_heat_flux_function();