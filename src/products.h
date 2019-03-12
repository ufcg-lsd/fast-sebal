#pragma once

#include "types.h"
#include "utils.h"
#include <algorithm>
#include <math.h>

string tal_function(TIFF* raster_elevation, string output_path);
void radiance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, int width_band, double radiance_line[][8]);
void reflectance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, double reflectance_line[][8]);
void albedo_function(double reflectance_line[][8], Sensor sensor, double tal_line[], int width_band, int number_sensor, double albedo_line[]);
void ndvi_function(double reflectance_line[][8], int width_band, double ndvi_line[]);
void lai_function(double reflectance_line[][8], int width_band, double lai_line[]);
void evi_function(double reflectance_line[][8], int width_band, double evi_line[]);
void enb_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double enb_emissivity_line[]);
void eo_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double eo_emissivity_line[]);
void ea_emissivity_function(double tal_line[], int width_band, double ea_emissivity_line[]);
void surface_temperature_function(double radiance_line[][8], double enb_emissivity_line[], int number_sensor, int width_band, double surface_temperature_line[]);
void short_wave_radiation_function(double tal_line[], MTL mtl, int width_band, double short_wave_radiation_line[]);
void large_wave_radiation_surface_function(double eo_emissivity_line[], double surface_temperature_line[], int width_band, double large_wave_radiation_surface_line[]);
void large_wave_radiation_atmosphere_function(double ea_emissivity_line[], int width_band, double temperature, double large_wave_radiation_atmosphere_line[]);
void net_radiation_function(double short_wave_radiation_line[], double large_wave_radiation_surface_line[],
                            double large_wave_radiation_atmosphere_line[], double albedo_line[],
                            double eo_emissivity_line[], int width_band, double net_radiation_line[]);
void soil_heat_flux_function(double ndvi_line[], double surface_temperature_line[], double albedo_line[], double net_radiation_line[], int width_band, double soil_heat_flux[]);
void ho_function(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]);
Candidate select_hot_pixel(TIFF* ndvi, TIFF* surface_temperature, TIFF* net_radiation, TIFF* soil_heat, int heigth_band, int width_band);
Candidate select_cold_pixel(TIFF* ndvi, TIFF* surface_temperature, TIFF* net_radiation, TIFF* soil_heat, int heigth_band, int width_band);
void sensible_heat_flux_function(Candidate hot_pixel, Candidate cold_pixel);