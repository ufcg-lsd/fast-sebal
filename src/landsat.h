#pragma once

#include "types.h"
#include "utils.h"
#include "products.h"
#include <vector>

struct Landsat{
    string tal_path, output_path;
    string albedo_path, ndvi_path;
    string evi_path, lai_path;
    string soil_heat_path, surface_temperature_path;
    string net_radiation_path;

    Landsat();
    Landsat(string tal_path, string output_path);

	void process_parcial_products(TIFF* read_bands[], MTL mtl, Station station, Sensor sensor);
    void process_final_products(vector<TIFF*> parcial_products);
    void create_tiffs(TIFF *tal, TIFF *albedo, TIFF *ndvi, TIFF *evi, TIFF *lai, TIFF *soil_heat, TIFF *surface_temperature, TIFF *net_radiation);
    void open_tiffs(TIFF *albedo, TIFF *ndvi, TIFF *evi, TIFF *lai, TIFF *soil_heat, TIFF *surface_temperature, TIFF *net_radiation);
    void save_tiffs(vector<double*> products_line, vector<TIFF*> products, int line);
};