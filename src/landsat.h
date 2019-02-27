#pragma once

#include "types.h"
#include "utils.h"
#include "products.h"

struct Landsat{
    string tal_path, output_path;
    string albedo_path, ndvi_path;
    string evi_path, lai_path;
    string soil_heat_path, surface_temperature_path;
    string net_radiation_path;

    Landsat();
    Landsat(string tal_path, string output_path);

	void process(TIFF* read_bands[], MTL mtl, Sensor sensor);
};