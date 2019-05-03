#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

string tal_function(TIFF* raster_elevation, string output_path);
void radiance_function(TIFF* read_bands[], MTL mtl, Sensor sensor, int width_band, int line, double radiance_line[][8]);
void reflectance_function(TIFF* read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, int line, double reflectance_line[][8]);
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
Candidate select_hot_pixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** net_radiation, TIFF** soil_heat, int heigth_band, int width_band);
Candidate select_cold_pixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** net_radiation, TIFF** soil_heat, int heigth_band, int width_band);
void zom_fuction(double A_ZOM, double B_ZOM, double ndvi_line[], int width_band, double zom_line[]);
void ustar_fuction(double u200, double zom_line[], int width_band, double ustar_line[]);
void aerodynamic_resistence_fuction(double ustar_line[], int width_band, double aerodynamic_resistence_line[]);
void sensible_heat_flux_function(Candidate hot_pixel, Candidate cold_pixel, double u200, double zom_line[], double ustar_line[], double aerodynamic_resistence_line[], double surface_temperature_line[], int width_band, double sensible_heat_flux_line[]);
void latent_heat_flux_function(double net_radiation_line[], double soil_heat_flux_line[], double sensible_heat_flux_line[], int width_band, double latent_heat_flux[]);
void net_radiation_24h_function(double albedo_line[], double Ra24h, double Rs24h, int width_band, double net_radiation_24h_line[]);
void evapotranspiration_fraction_fuction(double latent_heat_flux_line[], double net_radiation_line[], double soil_heat_line[], int width_band, double evapotranspiration_fraction_line[]);
void sensible_heat_flux_24h_fuction(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double sensible_heat_flux_24h_line[]);
void latent_heat_flux_24h_function(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double latent_heat_flux_24h_line[]);
void evapotranspiration_24h_function(double latent_heat_flux_24h_line[], Station station, int width_band, double evapotranspiration_24h_line[]);