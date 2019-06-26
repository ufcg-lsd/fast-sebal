#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

/**
 * @brief  
 * @note   
 * @param  raster_elevation: 
 * @param  output_path: 
 * @retval 
 */
string tal_function(TIFF* raster_elevation, string output_path);

/**
 * @brief  
 * @note   
 * @param  read_bands[]: 
 * @param  mtl: 
 * @param  sensor: 
 * @param  width_band: 
 * @param  line: 
 * @param  radiance_line[][8]: 
 * @retval None
 */
void radiance_function(TIFF* read_bands[], MTL mtl, Sensor sensor, int width_band, int line, double radiance_line[][8]);

/**
 * @brief  
 * @note   
 * @param  read_bands[]: 
 * @param  mtl: 
 * @param  sensor: 
 * @param  radiance_line[][8]: 
 * @param  width_band: 
 * @param  line: 
 * @param  reflectance_line[][8]: 
 * @retval None
 */
void reflectance_function(TIFF* read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, int line, double reflectance_line[][8]);

/**
 * @brief  
 * @note   
 * @param  reflectance_line[][8]: 
 * @param  sensor: 
 * @param  tal_line[]: 
 * @param  width_band: 
 * @param  number_sensor: 
 * @param  albedo_line[]: 
 * @retval None
 */
void albedo_function(double reflectance_line[][8], Sensor sensor, double tal_line[], int width_band, int number_sensor, double albedo_line[]);

/**
 * @brief  
 * @note   
 * @param  reflectance_line[][8]: 
 * @param  width_band: 
 * @param  ndvi_line[]: 
 * @retval None
 */
void ndvi_function(double reflectance_line[][8], int width_band, double ndvi_line[]);

/**
 * @brief  
 * @note   
 * @param  reflectance_line[][8]: 
 * @param  width_band: 
 * @param  lai_line[]: 
 * @retval None
 */
void lai_function(double reflectance_line[][8], int width_band, double lai_line[]);

/**
 * @brief  
 * @note   
 * @param  reflectance_line[][8]: 
 * @param  width_band: 
 * @param  evi_line[]: 
 * @retval None
 */
void evi_function(double reflectance_line[][8], int width_band, double evi_line[]);

/**
 * @brief  
 * @note   
 * @param  lai_line[]: 
 * @param  ndvi_line[]: 
 * @param  width_band: 
 * @param  enb_emissivity_line[]: 
 * @retval None
 */
void enb_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double enb_emissivity_line[]);

/**
 * @brief  
 * @note   
 * @param  lai_line[]: 
 * @param  ndvi_line[]: 
 * @param  width_band: 
 * @param  eo_emissivity_line[]: 
 * @retval None
 */
void eo_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double eo_emissivity_line[]);

/**
 * @brief  
 * @note   
 * @param  tal_line[]: 
 * @param  width_band: 
 * @param  ea_emissivity_line[]: 
 * @retval None
 */
void ea_emissivity_function(double tal_line[], int width_band, double ea_emissivity_line[]);

/**
 * @brief  
 * @note   
 * @param  radiance_line[][8]: 
 * @param  enb_emissivity_line[]: 
 * @param  number_sensor: 
 * @param  width_band: 
 * @param  surface_temperature_line[]: 
 * @retval None
 */
void surface_temperature_function(double radiance_line[][8], double enb_emissivity_line[], int number_sensor, int width_band, double surface_temperature_line[]);

/**
 * @brief  
 * @note   
 * @param  tal_line[]: 
 * @param  mtl: 
 * @param  width_band: 
 * @param  short_wave_radiation_line[]: 
 * @retval None
 */
void short_wave_radiation_function(double tal_line[], MTL mtl, int width_band, double short_wave_radiation_line[]);

/**
 * @brief  
 * @note   
 * @param  eo_emissivity_line[]: 
 * @param  surface_temperature_line[]: 
 * @param  width_band: 
 * @param  large_wave_radiation_surface_line[]: 
 * @retval None
 */
void large_wave_radiation_surface_function(double eo_emissivity_line[], double surface_temperature_line[], int width_band, double large_wave_radiation_surface_line[]);

/**
 * @brief  
 * @note   
 * @param  ea_emissivity_line[]: 
 * @param  width_band: 
 * @param  temperature: 
 * @param  large_wave_radiation_atmosphere_line[]: 
 * @retval None
 */
void large_wave_radiation_atmosphere_function(double ea_emissivity_line[], int width_band, double temperature, double large_wave_radiation_atmosphere_line[]);

/**
 * @brief  
 * @note   
 * @param  short_wave_radiation_line[]: 
 * @param  large_wave_radiation_surface_line[]: 
 * @param  large_wave_radiation_atmosphere_line[]: 
 * @param  albedo_line[]: 
 * @param  eo_emissivity_line[]: 
 * @param  width_band: 
 * @param  net_radiation_line[]: 
 * @retval None
 */
void net_radiation_function(double short_wave_radiation_line[], double large_wave_radiation_surface_line[],
                            double large_wave_radiation_atmosphere_line[], double albedo_line[],
                            double eo_emissivity_line[], int width_band, double net_radiation_line[]);

/**
 * @brief  
 * @note   
 * @param  ndvi_line[]: 
 * @param  surface_temperature_line[]: 
 * @param  albedo_line[]: 
 * @param  net_radiation_line[]: 
 * @param  width_band: 
 * @param  soil_heat_flux[]: 
 * @retval None
 */
void soil_heat_flux_function(double ndvi_line[], double surface_temperature_line[], double albedo_line[], double net_radiation_line[], int width_band, double soil_heat_flux[]);

/**
 * @brief  
 * @note   
 * @param  net_radiation_line[]: 
 * @param  soil_heat_flux[]: 
 * @param  width_band: 
 * @param  ho_line[]: 
 * @retval None
 */
void ho_function(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]);

/**
 * @brief  
 * @note   
 * @param  ndvi: 
 * @param  surface_temperature: 
 * @param  net_radiation: 
 * @param  soil_heat: 
 * @param  height_band: 
 * @param  width_band: 
 * @retval 
 */
Candidate select_hot_pixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);

/**
 * @brief  
 * @note   
 * @param  ndvi: 
 * @param  surface_temperature: 
 * @param  net_radiation: 
 * @param  soil_heat: 
 * @param  height_band: 
 * @param  width_band: 
 * @retval 
 */
Candidate select_cold_pixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);

/**
 * @brief  
 * @note   
 * @param  A_ZOM: 
 * @param  B_ZOM: 
 * @param  ndvi_line[]: 
 * @param  width_band: 
 * @param  zom_line[]: 
 * @retval None
 */
void zom_fuction(double A_ZOM, double B_ZOM, double ndvi_line[], int width_band, double zom_line[]);

/**
 * @brief  
 * @note   
 * @param  u200: 
 * @param  zom_line[]: 
 * @param  width_band: 
 * @param  ustar_line[]: 
 * @retval None
 */
void ustar_fuction(double u200, double zom_line[], int width_band, double ustar_line[]);

/**
 * @brief  
 * @note   
 * @param  ustar_line[]: 
 * @param  width_band: 
 * @param  aerodynamic_resistance_line[]: 
 * @retval None
 */
void aerodynamic_resistance_fuction(double ustar_line[], int width_band, double aerodynamic_resistance_line[]);

/**
 * @brief  
 * @note   
 * @param  net_radiation_line[]: 
 * @param  soil_heat_flux_line[]: 
 * @param  sensible_heat_flux_line[]: 
 * @param  width_band: 
 * @param  latent_heat_flux[]: 
 * @retval None
 */
void latent_heat_flux_function(double net_radiation_line[], double soil_heat_flux_line[], double sensible_heat_flux_line[], int width_band, double latent_heat_flux[]);

/**
 * @brief  
 * @note   
 * @param  albedo_line[]: 
 * @param  Ra24h: 
 * @param  Rs24h: 
 * @param  width_band: 
 * @param  net_radiation_24h_line[]: 
 * @retval None
 */
void net_radiation_24h_function(double albedo_line[], double Ra24h, double Rs24h, int width_band, double net_radiation_24h_line[]);

/**
 * @brief  
 * @note   
 * @param  latent_heat_flux_line[]: 
 * @param  net_radiation_line[]: 
 * @param  soil_heat_line[]: 
 * @param  width_band: 
 * @param  evapotranspiration_fraction_line[]: 
 * @retval None
 */
void evapotranspiration_fraction_fuction(double latent_heat_flux_line[], double net_radiation_line[], double soil_heat_line[], int width_band, double evapotranspiration_fraction_line[]);

/**
 * @brief  
 * @note   
 * @param  evapotranspiration_fraction_line[]: 
 * @param  net_radiation_24h_line[]: 
 * @param  width_band: 
 * @param  sensible_heat_flux_24h_line[]: 
 * @retval None
 */
void sensible_heat_flux_24h_fuction(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double sensible_heat_flux_24h_line[]);

/**
 * @brief  
 * @note   
 * @param  evapotranspiration_fraction_line[]: 
 * @param  net_radiation_24h_line[]: 
 * @param  width_band: 
 * @param  latent_heat_flux_24h_line[]: 
 * @retval None
 */
void latent_heat_flux_24h_function(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double latent_heat_flux_24h_line[]);

/**
 * @brief  
 * @note   
 * @param  latent_heat_flux_24h_line[]: 
 * @param  station: 
 * @param  width_band: 
 * @param  evapotranspiration_24h_line[]: 
 * @retval None
 */
void evapotranspiration_24h_function(double latent_heat_flux_24h_line[], Station station, int width_band, double evapotranspiration_24h_line[]);