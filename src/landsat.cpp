#include "landsat.h"

Landsat::Landsat(){
};

Landsat::Landsat(string tal_path, string output_path){
    this->tal_path = tal_path;
    this->output_path = output_path;
    this->albedo_path = output_path + "/albedo.tif";
    this->ndvi_path = output_path + "/ndvi.tif";
    this->evi_path = output_path + "/evi.tif";
    this->lai_path = output_path + "/lai.tif";
    this->soil_heat_path = output_path + "/G.tif";
    this->surface_temperature_path = output_path + "/TS.tif";
    this->net_radiation_path = output_path + "/Rn.tif";
};

void Landsat::process_parcial_products(TIFF* read_bands[], MTL mtl, Station station, Sensor sensor){
    uint32 heigth_band, width_band;

    PixelReader pixel_read_bands[8];
    tdata_t line_bands[8];

    for(int i = 1; i < 8; i++){
        uint16 sample_band;
        
        TIFFGetField(read_bands[i], TIFFTAG_SAMPLEFORMAT, &sample_band);
        TIFFGetField(read_bands[i], TIFFTAG_IMAGEWIDTH, &width_band);

        unsigned short byte_size_band = TIFFScanlineSize(read_bands[i]) / width_band;

        line_bands[i] = _TIFFmalloc(TIFFScanlineSize(read_bands[i]));

        pixel_read_bands[i] = PixelReader(sample_band, byte_size_band, line_bands[i]);
    }

    TIFFGetField(read_bands[1], TIFFTAG_IMAGEWIDTH, &width_band);
    TIFFGetField(read_bands[1], TIFFTAG_IMAGELENGTH, &heigth_band);

    TIFF *tal = TIFFOpen(this->tal_path.c_str(), "rm");
    check_open_tiff(tal);

    double tal_line[width_band];

    TIFF *albedo, *ndvi, *evi, *lai, *soil_heat, *surface_temperature, *net_radiation;
    create_tiffs(tal, albedo, ndvi, evi, lai, soil_heat, surface_temperature, net_radiation);

    //Declare array with product information
    double albedo_line[width_band], ndvi_line[width_band];
    double evi_line[width_band], lai_line[width_band];
    double soil_heat_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band];

    //Declare auxiliars arrays
    double radiance_line[width_band][8];
    double reflectance_line[width_band][8];

    //Declare arrays of auxiliars products
    double eo_emissivity_line[width_band], ea_emissivity_line[width_band], enb_emissivity_line[width_band];
    double large_wave_radiation_atmosphere_line[width_band], large_wave_radiation_surface_line[width_band];
    double short_wave_radiation_line[width_band];

    for(int line = 0; line < heigth_band; line ++){
        radiance_function(pixel_read_bands, mtl, sensor, width_band, radiance_line);
        reflectance_function(pixel_read_bands, mtl, sensor, radiance_line, width_band, reflectance_line);

        read_line_tiff(tal, tal_line, line);

        albedo_function(reflectance_line, sensor, tal_line, width_band, mtl.number_sensor, albedo_line);
        short_wave_radiation_function(tal_line, mtl, width_band, short_wave_radiation_line);
        ndvi_function(reflectance_line, width_band, ndvi_line);
        lai_function(reflectance_line, width_band, lai_line);
        evi_function(reflectance_line, width_band, evi_line);
        enb_emissivity_function(lai_line, ndvi_line, width_band, enb_emissivity_line);
        eo_emissivity_function(lai_line, ndvi_line, width_band, eo_emissivity_line);
        surface_temperature_function(radiance_line, enb_emissivity_line, mtl.number_sensor, width_band, surface_temperature_line);
        large_wave_radiation_surface_function(eo_emissivity_line, surface_temperature_line, width_band, large_wave_radiation_surface_line);
        ea_emissivity_function(tal_line, width_band, ea_emissivity_line);
        large_wave_radiation_atmosphere_function(ea_emissivity_line, width_band, station.temperature_image, large_wave_radiation_atmosphere_line);
        net_radiation_function(short_wave_radiation_line, large_wave_radiation_surface_line, large_wave_radiation_atmosphere_line, albedo_line, eo_emissivity_line, width_band, net_radiation_line);
        soil_heat_flux_function(ndvi_line, surface_temperature_line, albedo_line, net_radiation_line, width_band, soil_heat_line);

        save_tiffs(vector<double*> {albedo_line, ndvi_line, evi_line, lai_line, soil_heat_line, surface_temperature_line, net_radiation_line},
                    vector<TIFF*> {albedo, ndvi, evi, lai, soil_heat, surface_temperature, net_radiation}, line);
    }

    TIFFClose(albedo);
    TIFFClose(ndvi);
    TIFFClose(evi);
    TIFFClose(lai);
    TIFFClose(soil_heat);
    TIFFClose(surface_temperature);
    TIFFClose(net_radiation);
    TIFFClose(tal);
};

void Landsat::process_final_products(Station station){
    TIFF *albedo, *ndvi, *evi, *lai, *soil_heat, *surface_temperature, *net_radiation;
    open_tiffs(albedo, ndvi, evi, lai, soil_heat, surface_temperature, net_radiation);

    uint32 heigth_band, width_band;
    TIFFGetField(albedo, TIFFTAG_IMAGELENGTH, &heigth_band);
    TIFFGetField(albedo, TIFFTAG_IMAGEWIDTH, &width_band);

    Candidate hot_pixel = select_hot_pixel(ndvi, surface_temperature, net_radiation, soil_heat, heigth_band, width_band);
    Candidate cold_pixel = select_cold_pixel(ndvi, surface_temperature, net_radiation, soil_heat, heigth_band, width_band);

    double sensible_heat_flux_line[band_width];

    for(int line = 0; line < heigth_band; line++){

        /*  TODO
            station.v6, n sei como acessar, n sei como ta funcionado o station
        */
        double ustar_station = (VON_KARMAN * station.v6)/(log(station.WIND_SPEED/station.SURFACE_ROUGHNESS));
        double u200 = ustar_station/(VON_KARMAN * log(200 / station.SURFACE_ROUGHNESS));

        hot_pixel.setAerodynamicResistence(u200, station.A_ZOM, station.B_ZOM, VON_KARMAN);
        cold_pixel.setAerodynamicResistence(u200, station.A_ZOM, station.B_ZOM, VON_KARMAN);

        ustar_fuction(); //TODO
        aerodynamic_resistence_fuction(); //TODO
        sensible_heat_flux_function(hot_pixel, cold_pixel, width_band);

    }

    
};

void Landsat::create_tiffs(TIFF *tal, TIFF *albedo, TIFF *ndvi, TIFF *evi, TIFF *lai, TIFF *soil_heat, TIFF *surface_temperature, TIFF *net_radiation){
    albedo = TIFFOpen(albedo_path.c_str(), "w8m");
    setup(albedo, tal);

    ndvi = TIFFOpen(ndvi_path.c_str(), "w8m");
    setup(ndvi, tal);

    evi = TIFFOpen(evi_path.c_str(), "w8m");
    setup(evi, tal);

    lai = TIFFOpen(lai_path.c_str(), "w8m");
    setup(lai, tal);

    soil_heat = TIFFOpen(soil_heat_path.c_str(), "w8m");
    setup(soil_heat, tal);

    surface_temperature = TIFFOpen(surface_temperature_path.c_str(), "w8m");
    setup(surface_temperature, tal);

    net_radiation = TIFFOpen(net_radiation_path.c_str(), "w8m");
    setup(net_radiation, tal);
};

void Landsat::open_tiffs(TIFF *albedo, TIFF *ndvi, TIFF *evi, TIFF *lai, TIFF *soil_heat, TIFF *surface_temperature, TIFF *net_radiation){

    albedo = TIFFOpen(albedo_path.c_str(), "rm");
    ndvi = TIFFOpen(ndvi_path.c_str(), "rm");
    evi = TIFFOpen(evi_path.c_str(), "rm");
    lai = TIFFOpen(lai_path.c_str(), "rm");
    soil_heat = TIFFOpen(soil_heat_path.c_str(), "rm");
    surface_temperature = TIFFOpen(surface_temperature_path.c_str(), "rm");
    net_radiation = TIFFOpen(net_radiation_path.c_str(), "rm");

}

void Landsat::save_tiffs(vector<double*> products_line, vector<TIFF*> products, int line){

    for (int i = 0; i < 7; i++){
        write_line_tiff(products[i], products_line[i], line);
    }
};