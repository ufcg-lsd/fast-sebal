#pragma once

#include "types.h"
#include "utils.h"
#include "parameters.h"
#include "pixel_reader.h"
#include "candidate.h"

/**
 * @brief  The spectral radiance for each band is computed.
 * @param  read_bands[]: Satellite bands.
 * @param  mtl: MTL struct.
 * @param  sensor: Sensor struct.
 * @param  width_band: Band width.
 * @param  line: Line to be calculated.
 * @param  radiance_line[][8]: Auxiliary array for save the calculated value of radiance for each band.
 */
void radiance_function(TIFF* read_bands[], MTL mtl, Sensor sensor, int width_band, int line, double radiance_line[][8], double noData);

/**
 * @brief  The reflectivity for each band (ρλ) is computed (model F02). 
 * @param  read_bands[]: Satellite bands.
 * @param  mtl: MTL struct.
 * @param  sensor: Sensor struct.
 * @param  radiance_line[][8]: Radiance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  line: Line to be calculated.
 * @param  reflectance_line[][8]: Auxiliary array for save the calculated value of reflectance for each band.
 */
void reflectance_function(TIFF* read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, int line, double reflectance_line[][8], double noData);

/**
 * @brief  The surface albedo is computed.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  sensor: Sensor struct.
 * @param  tal_line[]: Array containing the specified line from the tal TIFF.
 * @param  width_band: Band width.
 * @param  number_sensor: Number of the satellite sensor.
 * @param  albedo_line[]: Auxiliary array for save the calculated value of albedo for the line.
 */
void albedo_function(double reflectance_line[][8], Sensor sensor, double tal_line[], int width_band, int number_sensor, double albedo_line[]);

/**
 * @brief  Normalized Difference Vegetation Index (NDVI) is computed.
 * @note   Values for NDVI range between -1 and 1.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  ndvi_line[]: Auxiliary array for save the calculated value of NDVI for the line.
 */
void ndvi_function(double reflectance_line[][8], int width_band, double ndvi_line[]);

/**
 * @brief  Leaf Area Index (LAI) is computed.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  lai_line[]: Auxiliary array for save the calculated value of LAI for the line.
 */
void lai_function(double reflectance_line[][8], int width_band, double lai_line[]);

/**
 * @brief  Enhanced Vegetation Index (EVI) is computed.
 * @param  reflectance_line[][8]: Reflectance for the specific line for each band.
 * @param  width_band: Band width.
 * @param  evi_line[]: Auxiliary array for save the calculated value of EVI for the line.
 */
void evi_function(double reflectance_line[][8], int width_band, double evi_line[]);

/**
 * @brief  Calculates emissivity representing surface behavior for thermal emission in the relatively narrow band 6 of Landsat (10.4 to 12.5 µm),
           expressed as enb.
 * @param  lai_line[]: Array containing the specified line from the LAI computation.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  enb_emissivity_line[]: Auxiliary array for save the calculated value of Enb for the line.
 */
void enb_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double enb_emissivity_line[]);

/**
 * @brief  Calculates emissivity representing surface behavior for thermal emission in the broad thermal spectrum (6 to 14 µm), expressed as eο.  
 * @param  lai_line[]: Array containing the specified line from the LAI computation.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  eo_emissivity_line[]: Auxiliary array for save the calculated value of Eo for the line.
 */
void eo_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double eo_emissivity_line[]);

/**
 * @brief  Calculates the atmospheric emissivity (ea).
 * @param  tal_line[]: Array containing the specified line from the tal computation.
 * @param  width_band: Band width.
 * @param  ea_emissivity_line[]: Auxiliary array for save the calculated value of Ea for the line.
 */
void ea_emissivity_function(double tal_line[], int width_band, double ea_emissivity_line[]);

/**
 * @brief  The surface temperature (TS) is computed.
 * @param  radiance_line[][8]: Radiance for the specific line for each band.
 * @param  enb_emissivity_line[]: Array containing the specified line from the Enb computation.
 * @param  number_sensor: Number of the satellite sensor.
 * @param  width_band: Band width.
 * @param  surface_temperature_line[]: Auxiliary array for save the calculated value of TS for the line.
 */
void surface_temperature_function(double radiance_line[][8], double enb_emissivity_line[], int number_sensor, int width_band, double surface_temperature_line[]);

/**
 * @brief  Computes Short Wave Radiation (Rs).
 * @param  tal_line[]: Array containing the specified line from the tal computation.
 * @param  mtl: MTL Struct.
 * @param  width_band: Band width.
 * @param  short_wave_radiation_line[]: Auxiliary array for save the calculated value of Rs for the line.
 */
void short_wave_radiation_function(double tal_line[], MTL mtl, int width_band, double short_wave_radiation_line[]);

/**
 * @brief  Computes Large Wave Radiation from Surface (RLSup)
 * @param  eo_emissivity_line[]: Array containing the specified line from the Eo computation.
 * @param  surface_temperature_line[]: Array containing the specified line from the TS computation.
 * @param  width_band: Band width.
 * @param  large_wave_radiation_surface_line[]: Auxiliary array for save the calculated value of RLSup for the line.
 */
void large_wave_radiation_surface_function(double eo_emissivity_line[], double surface_temperature_line[], int width_band, double large_wave_radiation_surface_line[]);

/**
 * @brief  Computes Large Wave Radiation from Atmosphere (RLatm)
 * @param  ea_emissivity_line[]: Array containing the specified line from the Ea computation.
 * @param  width_band: Band width.
 * @param  temperature: Near surface air temperature in Kelvin.
 * @param  large_wave_radiation_atmosphere_line[]: Auxiliary array for save the calculated value of RLatm for the line.
 */
void large_wave_radiation_atmosphere_function(double ea_emissivity_line[], int width_band, double temperature, double large_wave_radiation_atmosphere_line[]);

/**
 * @brief  The net surface radiation flux (Rn) is computed.
 * @param  short_wave_radiation_line[]: Array containing the specified line from the Rs computation.
 * @param  large_wave_radiation_surface_line[]: Array containing the specified line from the RLSup computation.
 * @param  large_wave_radiation_atmosphere_line[]: Array containing the specified line from the RLatm computation.
 * @param  albedo_line[]: Array containing the specified line from the albedo computation.
 * @param  eo_emissivity_line[]: Array containing the specified line from the Eo computation.
 * @param  width_band: Band width.
 * @param  net_radiation_line[]: Auxiliary array for save the calculated value of Rn for the line.
 */
void net_radiation_function(double short_wave_radiation_line[], double large_wave_radiation_surface_line[],
                            double large_wave_radiation_atmosphere_line[], double albedo_line[],
                            double eo_emissivity_line[], int width_band, double net_radiation_line[]);

/**
 * @brief  Computes the Soil heat flux (G).    
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  surface_temperature_line[]: Array containing the specified line from the TS computation.
 * @param  albedo_line[]: Array containing the specified line from the albedo computation.
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  width_band: Band width.
 * @param  soil_heat_flux[]: Auxiliary array for save the calculated value of G for the line.
 */
void soil_heat_flux_function(double ndvi_line[], double surface_temperature_line[], double albedo_line[], double net_radiation_line[], int width_band, double soil_heat_flux[]);

/**
 * @brief  Computes the HO.
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  soil_heat_flux[]: Array containing the specified line from the G computation.
 * @param  width_band: Band width.
 * @param  ho_line[]: Auxiliary array for save the calculated value of HO for the line.
 */
void ho_function(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]);

/**
 * @brief  Select the hot pixel.
 * @param  ndvi: NDVI TIFF.
 * @param  surface_temperature: TS TIFF.
 * @param  net_radiation: Rn TIFF.
 * @param  soil_heat: G TIFF.
 * @param  height_band: Band height.
 * @param  width_band: Band width.
 * @retval Candidate struct containing the hot pixel.
 */
Candidate select_hot_pixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);

/**
 * @brief  Select the cold pixel.
 * @param  ndvi: NDVI TIFF.
 * @param  surface_temperature: TS TIFF.
 * @param  net_radiation: Rn TIFF.
 * @param  soil_heat: G TIFF.
 * @param  height_band: Band height.
 * @param  width_band: Band width.
 * @retval Candidate struct containing the cold pixel.
 */
Candidate select_cold_pixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band);

/**
 * @brief  Computes the momentum roughness length (zom).
 * @param  A_ZOM: Correlation constant a.
 * @param  B_ZOM: Correlation constant b.
 * @param  ndvi_line[]: Array containing the specified line from the NDVI computation.
 * @param  width_band: Band width.
 * @param  zom_line[]: Auxiliary array for save the calculated value of zom for the line.
 */
void zom_fuction(double A_ZOM, double B_ZOM, double ndvi_line[], int width_band, double zom_line[]);

/**
 * @brief  The friction velocity (u*) is computed.
 * @param  u200: Wind speed at 200 m.
 * @param  zom_line[]: Array containing the specified line from the zom computation.
 * @param  width_band: Band width.
 * @param  ustar_line[]: Auxiliary array for save the calculated value of ustar for the line.
 */
void ustar_fuction(double u200, double zom_line[], int width_band, double ustar_line[]);

/**
 * @brief  Computes the aerodynamic resistance (Rah).   
 * @param  ustar_line[]: Array containing the specified line from the ustar computation.
 * @param  width_band: Band width.
 * @param  aerodynamic_resistance_line[]: Auxiliary array for save the calculated value of Rah for the line.
 */
void aerodynamic_resistance_fuction(double ustar_line[], int width_band, double aerodynamic_resistance_line[]);

/**
 * @brief  Computes Latent Heat Flux (LE).  
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  soil_heat_flux_line[]: Array containing the specified line from the G computation.
 * @param  sensible_heat_flux_line[]: Array containing the specified line from the H computation.
 * @param  width_band: Band width.
 * @param  latent_heat_flux[]: Auxiliary array for save the calculated value of LE for the line.
 */
void latent_heat_flux_function(double net_radiation_line[], double soil_heat_flux_line[], double sensible_heat_flux_line[], int width_band, double latent_heat_flux[]);

/**
 * @brief  Calculates the Net Radiation for 24 hours (Rn24h).
 * @param  albedo_line[]: Array containing the specified line from the albedo computation.
 * @param  Ra24h: Ra24h: Extraterrestrial Radiation defined as solar short wave radiation in the absence of an atmosphere (Ra24h).
 * @param  Rs24h: Short wave radiation incident in 24 hours (Rs24h).
 * @param  width_band: Band width.
 * @param  net_radiation_24h_line[]: Auxiliary array for save the calculated value of Rn24h for the line.
 */
void net_radiation_24h_function(double albedo_line[], double Ra24h, double Rs24h, int width_band, double net_radiation_24h_line[]);

/**
 * @brief  The Reference ET Fraction (EF) is computed.
 * @param  latent_heat_flux_line[]: Array containing the specified line from the LE computation.
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  soil_heat_line[]: Array containing the specified line from the G computation.
 * @param  width_band: Band width.
 * @param  evapotranspiration_fraction_line[]: Auxiliary array for save the calculated value of EF for the line.
 */
void evapotranspiration_fraction_fuction(double latent_heat_flux_line[], double net_radiation_line[], double soil_heat_line[], int width_band, double evapotranspiration_fraction_line[]);

/**
 * @brief  Computes Sensible Heat Flux for 24 hours (H24h).
 * @param  evapotranspiration_fraction_line[]: Array containing the specified line from the EF computation.
 * @param  net_radiation_24h_line[]: Array containing the specified line from the Rn24h computation.
 * @param  width_band: Band width.
 * @param  sensible_heat_flux_24h_line[]: Auxiliary array for save the calculated value of H24h for the line.
 */
void sensible_heat_flux_24h_fuction(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double sensible_heat_flux_24h_line[]);

/**
 * @brief  Calculates Latente Heat Flux for 24 hours (LE24h).
 * @param  evapotranspiration_fraction_line[]: Array containing the specified line from the EF computation.
 * @param  net_radiation_24h_line[]: Array containing the specified line from the Rn24h computation.
 * @param  width_band: Band width.
 * @param  latent_heat_flux_24h_line[]: Auxiliary array for save the calculated value of LE24h for the line.
 */
void latent_heat_flux_24h_function(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double latent_heat_flux_24h_line[]);

/**
 * @brief  Computes the Evapotranspiration for 24 hours (ET24h)
 * @param  latent_heat_flux_24h_line[]: Array containing the specified line from the LE24h computation.
 * @param  station: Station struct.
 * @param  width_band: Band width.
 * @param  evapotranspiration_24h_line[]: Auxiliary array for save the calculated value of ET24h for the line.
 */
void evapotranspiration_24h_function(double latent_heat_flux_24h_line[], Station station, int width_band, double evapotranspiration_24h_line[]);