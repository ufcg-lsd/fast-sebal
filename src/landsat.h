#pragma once

#include "types.h"
#include "utils.h"
#include "products.h"
#include "parameters.h"
#include "asebal.h"

/**
 * @brief  Struct to manage the products calculation.
 */
struct Landsat{

    // Path to TIFF files.
    string tal_path, output_path;
    string albedo_path, ndvi_path;
    string evi_path, lai_path;
    string soil_heat_path, surface_temperature_path;
    string net_radiation_path, evapotranspiration_fraction_path;
    string evapotranspiration_24h_path;
    string zom_path, ustar_path, aerodynamic_resistance_path, sensible_heat_flux_path;
    string ustar_tif0_path, ustar_tif1_path, aerodynamic_resistance_tif0_path, aerodynamic_resistance_tif1_path;
    string latent_heat_flux_path, net_radiation_24h_path, latent_heat_flux_24h_path;

    /**
     * @brief  Empty constructor.
     */
    Landsat();

    /**
     * @brief  Constructor of the struct.
     * @param  tal_path: Path to tal TIFF.
     * @param  output_path: Output path where TIFF should be saved.
     */
    Landsat(string tal_path, string output_path);

	/**
	 * @brief  Calculates the partials products (e. g. Albedo, NDVI, Rn, G) of the SEBAL execution.
	 * @param  read_bands[]: Satellite images as TIFFs.
	 * @param  mtl: MTL struct.
	 * @param  station: Station struct.
	 * @param  sensor: Sensor struct.
	 */
	void process_partial_products(TIFF* read_bands[], MTL mtl, Station station, Sensor sensor);

    /**
     * @brief  Process the final products (e. g. Evapotranspiration 24 hours) of the SEBAL execution.
     * @param  station: Station struct.
	 * @param  mtl: MTL struct.
     */
    void process_final_products(Station station, MTL mtl);

    /**
     * @brief  Initializes TIFFs of the partial execution products as writable. Doing their setup based upon the tal TIFF characteristics.
     * @param  **tal: Tal TIFF.
     * @param  **albedo: Albedo TIFF.
     * @param  **ndvi: NDVI TIFF.
     * @param  **evi: EVI TIFF.
     * @param  **lai: LAI TIFF.
     * @param  **soil_heat: Soil heat flux TIFF.
     * @param  **surface_temperature: Surface temperature TIFF.
     * @param  **net_radiation: Net radiation TIFF.
     */
    void create_tiffs(TIFF **tal, TIFF **albedo, TIFF **ndvi, TIFF **evi, TIFF **lai, TIFF **soil_heat, TIFF **surface_temperature, TIFF **net_radiation);
    
    /**
     * @brief  Open the partial products TIFF as readble TIFFs and create the final products TIFF (evapotranspiration_fraction and evapotranspiration_24h) as writable. For use them at the final phase.
     * @param  **albedo: Albedo TIFF.
     * @param  **ndvi: NDVI TIFF. 
     * @param  **soil_heat: Soil heat flux TIFF.
     * @param  **surface_temperature: Surface temperature TIFF.
     * @param  **net_radiation: Net radiation TIFF.
     * @param  **evapotranspiration_fraction: Evapotranspiration fraction TIFF.
     * @param  **evapotranspiration_24h: Evapotranspiration 24 hours TIFF.
     */
    void open_tiffs(TIFF **albedo, TIFF **ndvi, TIFF **soil_heat, TIFF **surface_temperature, TIFF **net_radiation, TIFF **evapotranspiration_fraction, TIFF **evapotranspiration_24h);
    
    /**
     * @brief  Writes values from an array to a specific line in a TIFF. Doing this for each respective array and TIFF at the vectors parameters passed.
     * @note:  The positions both vectors should be corresponding arrays and TIFFs.
     * @param  products_line: Vector containing the arrays of a line to be written on a respective TIFF.
     * @param  products: Vector containing the respective TIFF for each array.
     * @param  line: Number of the line that should be written.
     */
    void save_tiffs(vector<double*> products_line, vector<TIFF*> products, int line);
};