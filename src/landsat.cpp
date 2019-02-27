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

void Landsat::process(TIFF* read_bands[], MTL mtl, Sensor sensor){
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

    TIFF* tal = TIFFOpen(this->tal_path.c_str(), "rm");
    double tal_line[width_band];
    double radiance_line[width_band][8];
    double reflectance_line[width_band][8];

    TIFF* albedo = TIFFOpen(albedo_path.c_str(), "w8m");
    setup(albedo, tal);

    TIFF* ndvi = TIFFOpen(ndvi_path.c_str(), "w8m");
    setup(ndvi, tal);

    TIFF* evi = TIFFOpen(evi_path.c_str(), "w8m");
    setup(evi, tal);

    TIFF* lai = TIFFOpen(lai_path.c_str(), "w8m");
    setup(lai, tal);

    TIFF* soil_heat = TIFFOpen(soil_heat_path.c_str(), "w8m");
    setup(soil_heat, tal);

    TIFF* surface_temperature = TIFFOpen(surface_temperature_path.c_str(), "w8m");
    setup(surface_temperature, tal);

    TIFF* net_radiation = TIFFOpen(net_radiation_path.c_str(), "w8m");
    setup(net_radiation, tal);

    for(int line = 0; line < heigth_band; line ++){
        radiance_function(pixel_read_bands, mtl, sensor, width_band, radiance_line);
        reflectance_function(pixel_read_bands, mtl, sensor, radiance_line, width_band, reflectance_line);

        if(TIFFReadScanline(tal, tal_line, line) < 0){
            cerr << "Read problem" << endl;
            exit(3);
        }

        albedo_function(reflectance_line, sensor, tal_line, width_band, mtl.number_sensor, line, albedo);
        
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