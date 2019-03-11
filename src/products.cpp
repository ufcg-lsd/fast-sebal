#include "products.h"

string tal_function(TIFF *raster_elevation, string output_path){
    uint32 heigth_band, width_band;
    uint16 sample_band;

    TIFFGetField(raster_elevation, TIFFTAG_IMAGEWIDTH, &width_band);
    TIFFGetField(raster_elevation, TIFFTAG_IMAGELENGTH, &heigth_band);
    TIFFGetField(raster_elevation, TIFFTAG_SAMPLEFORMAT, &sample_band);

    string tal_path = output_path + "/tal.tif";
    TIFF *tal = TIFFOpen(tal_path.c_str(), "w8m");
    setup(tal, raster_elevation);

    PixelReader pixel_read_band;
    tdata_t line_band;
    double tal_line_band[width_band];

    unsigned short byte_size_band = TIFFScanlineSize(raster_elevation) / width_band;
    line_band = _TIFFmalloc(TIFFScanlineSize(raster_elevation));
    pixel_read_band = PixelReader(sample_band, byte_size_band, line_band);

    for (int line = 0; line < heigth_band; line++){
        if (TIFFReadScanline(raster_elevation, line_band, line) < 0){
            cerr << "Read problem" << endl;
            exit(3);
        }

        for (int col = 0; col < width_band; col++){
            double pixel_read = pixel_read_band.read_pixel(col);

            if (fabs(pixel_read - 0) <= EPS)
                tal_line_band[col] = NaN;
            else
                tal_line_band[col] = 0.75 + 2 * 1e-5 * pixel_read;
        }

        if (TIFFWriteScanline(tal, tal_line_band, line) < 0){
            cerr << "Write problem" << endl;
            exit(4);
        }
    }

    TIFFClose(tal);

    return tal_path;
}

void radiance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, int width_band, double radiance_line[][8]){
    if (mtl.number_sensor == 8){
        for (int col = 0; col < width_band; col++){
            if (isnan(pixel_read_bands[7].read_pixel(col))){
                radiance_line[col][7] = NaN;
                continue;
            }
            radiance_line[col][7] = pixel_read_bands[7].read_pixel(col) * mtl.rad_mult_10 + mtl.rad_add_10;
        }
    }
    else{
        for (int i = 1; i < 8; i++){
            for (int col = 0; col < width_band; col++){
                if (isnan(pixel_read_bands[i].read_pixel(col))){
                    radiance_line[col][i] = NaN;
                    continue;
                }
                radiance_line[col][i] = pixel_read_bands[i].read_pixel(col) * sensor.parameters[i][sensor.GRESCALE] + sensor.parameters[i][sensor.BRESCALE];
            }
        }
    }
};

void reflectance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, double reflectance_line[][8]){
    double costheta = sin(mtl.sun_elevation * PI / 180);

    for (int i = 1; i < 8; i++){
        for (int col = 0; col < width_band; col++){
            if (isnan(pixel_read_bands[i].read_pixel(col))){
                reflectance_line[col][i] = NaN;
                continue;
            }

            if (mtl.number_sensor == 8){
                reflectance_line[col][i] = (pixel_read_bands[i].read_pixel(col) * sensor.parameters[i][sensor.GRESCALE] + sensor.parameters[i][sensor.BRESCALE]) / costheta;
            }
            else{
                reflectance_line[col][i] = (PI * radiance_line[col][i] * mtl.distance_earth_sun * mtl.distance_earth_sun) /
                                           (sensor.parameters[i][sensor.ESUN] * costheta);
            }
        }
    }
};

void albedo_function(double reflectance_line[][8], Sensor sensor, double tal_line[], int width_band, int number_sensor, int line, double albedo_line[]){
    int final_tif_calc = number_sensor == 8 ? 6 : 7;
    double albedo_line[width_band];

    for (int col = 0; col < width_band; col++){
        albedo_line[col] = reflectance_line[col][1] * sensor.parameters[1][sensor.WB] +
                            reflectance_line[col][2] * sensor.parameters[2][sensor.WB] +
                            reflectance_line[col][3] * sensor.parameters[3][sensor.WB] +
                            reflectance_line[col][4] * sensor.parameters[4][sensor.WB] +
                            reflectance_line[col][5] * sensor.parameters[5][sensor.WB] +
                            reflectance_line[col][final_tif_calc] * sensor.parameters[final_tif_calc][sensor.WB];
        albedo_line[col] = (albedo_line[col] - 0.03) / (tal_line[col] * tal_line[col]);
    }

};//

void ndvi_function(double reflectance_line[][8], int width_band, int line, TIFF *ndvi){
    double ndvi_line[width_band];

    for (int col = 0; col < width_band; col++){
        ndvi_line[col] = (reflectance_line[col][4] - reflectance_line[col][3]) /
                         (reflectance_line[col][4] + reflectance_line[col][3]);
    }

    if (TIFFWriteScanline(ndvi, ndvi_line, line) < 0){
        cerr << "Write problem in ndvi tif" << endl;
        exit(4);
    }
};//

void lai_function(double reflectance_line[][8], int width_band, int line, TIFF *lai){
    double savi_line[width_band];
    double lai_line[width_band];

    double L = 0.05;

    for (int col = 0; col < width_band; col++){

        savi_line[col] = ((1 + L) * (reflectance_line[col][4] - reflectance_line[col][3])) /
                         (L + (reflectance_line[col][4] + reflectance_line[col][3]));
        
        if (!isnan(savi_line[col]) && savi_line[col] > 0.687){
            lai_line[col] = 6;
        }
        else if (!isnan(savi_line[col]) && savi_line[col] < 0.1){
            lai_line[col] = 0;
        }
        else if (!isnan(savi_line[col])){
            lai_line[col] = -log((0.69 - savi_line[col]) / 0.59) / 0.91;
        }
    }

    if (TIFFWriteScanline(lai, lai_line, line) < 0){
        cerr << "Write problem in lai tif" << endl;
        exit(4);
    }
};//

void evi_function(double reflectance_line[][8], int width_band, int line, TIFF *evi){
    double evi_line[width_band];

    for (int col = 0; col < width_band; col++){
        evi_line[col] = 2.5 * ((reflectance_line[col][4] - reflectance_line[col][3])/
                        (reflectance_line[col][4] + 6 * reflectance_line[col][3] - 7.5 * reflectance_line[col][1] + 1));
    }

    if (TIFFWriteScanline(evi, evi_line, line) < 0){
        cerr << "Write problem in evi tif" << endl;
        exit(4);
    }
};//

void enb_emissivity_function(){

}; //usa LAI //so eh usada no TS

void eo_emissivity_function(){

}; //usa LAI

void ea_emissivity_function(){

}; //usa tal //so eh usada no RLat

void kelvin_surface_temperature_function(){

}; //usa enb e radiance

void short_wave_radiation_function(double tal_line[], MTL mtl, int width_band, double short_wave_radiation_line[]){
    double costheta = sin(mtl.sun_elevation * PI / 180);

    for (int col = 0; col < width_band; col++){
        short_wave_radiation_line[col] = (1367 * costheta * tal_line[col]) /
                                        (mtl.distance_earth_sun * mtl.distance_earth_sun);
    }
}; //Rs //usa costheta, tal e d_sun_eart

void large_wave_radiation_surface_function(){

}; //RLsup //usa eo e TS

void large_wave_radiation_atmosphere_function(){

}; //RLat //usa ea e .station

void net_radiation_function(double tal_line[], double albedo_line[], MTL mtl, int width_band, int line, TIFF *net_radiation){
    double net_radiation_line[width_band];
    double short_wave_radiation_line[width_band];
    double large_wave_radiation_surface_line[width_band];
    double large_wave_radiation_atmosphere_line[width_band];
    double eo_emissivity_line[width_band];

    short_wave_radiation_function(tal_line, mtl, width_band, short_wave_radiation_line);
    large_wave_radiation_atmosphere_function();
    large_wave_radiation_surface_function();
    eo_emissivity_function();

    for (int col = 0; col < width_band; col++){
        net_radiation_line[col] = short_wave_radiation_line[col] - (short_wave_radiation_line[col] * albedo_line[col]) +
                                large_wave_radiation_atmosphere_line[col] - large_wave_radiation_surface_line[col] -
                                (1 - eo_emissivity_line[col]) * large_wave_radiation_atmosphere_line[col];
        if(net_radiation_line[col] < 0)
            net_radiation_line[col] = 0;
    }

    if (TIFFWriteScanline(net_radiation, net_radiation_line, line) < 0){
        cerr << "Write problem in evi tif" << endl;
        exit(4);
    }

}; //Rn //usa rs, alb, RLat, RLsup, eo

void soil_heat_flux_function(){

}; //G //usa ndvi, ts, alb e rn