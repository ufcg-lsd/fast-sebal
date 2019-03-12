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
        read_line_tiff(raster_elevation, line_band, line);

        for (int col = 0; col < width_band; col++){
            double pixel_read = pixel_read_band.read_pixel(col);

            if (fabs(pixel_read - 0) <= EPS)
                tal_line_band[col] = NaN;
            else
                tal_line_band[col] = 0.75 + 2 * 1e-5 * pixel_read;
        }

        write_line_tiff(tal, tal_line_band, line);
    }

    TIFFClose(tal);

    return tal_path;
}; //tal

void radiance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, int width_band, double radiance_line[][8]){
    if (mtl.number_sensor == 8){
        for (int col = 0; col < width_band; col++){
            double pixel_value = pixel_read_bands[7].read_pixel(col);
            radiance_line[col][7] = pixel_value * mtl.rad_mult_10 + mtl.rad_add_10;
        }
    }
    else{
        for (int i = 1; i < 8; i++){
            for (int col = 0; col < width_band; col++){
                double pixel_value = pixel_read_bands[7].read_pixel(col);
                radiance_line[col][i] = pixel_value * sensor.parameters[i][sensor.GRESCALE] + sensor.parameters[i][sensor.BRESCALE];
            }
        }
    }
}; //rad

void reflectance_function(PixelReader pixel_read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, double reflectance_line[][8]){
    double costheta = sin(mtl.sun_elevation * PI / 180);

    for (int i = 1; i < 8; i++){
        for (int col = 0; col < width_band; col++){
            if (mtl.number_sensor == 8){
                double pixel_value = pixel_read_bands[i].read_pixel(col);
                reflectance_line[col][i] = (pixel_value * sensor.parameters[i][sensor.GRESCALE] + sensor.parameters[i][sensor.BRESCALE]) / costheta;
            }else{
                reflectance_line[col][i] = (PI * radiance_line[col][i] * mtl.distance_earth_sun * mtl.distance_earth_sun) /
                                           (sensor.parameters[i][sensor.ESUN] * costheta);
            }
        }
    }
}; //ref

void albedo_function(double reflectance_line[][8], Sensor sensor, double tal_line[], int width_band, int number_sensor, double albedo_line[]){
    int final_tif_calc = number_sensor == 8 ? 6 : 7;

    for (int col = 0; col < width_band; col++){
        albedo_line[col] = reflectance_line[col][1] * sensor.parameters[1][sensor.WB] +
                            reflectance_line[col][2] * sensor.parameters[2][sensor.WB] +
                            reflectance_line[col][3] * sensor.parameters[3][sensor.WB] +
                            reflectance_line[col][4] * sensor.parameters[4][sensor.WB] +
                            reflectance_line[col][5] * sensor.parameters[5][sensor.WB] +
                            reflectance_line[col][final_tif_calc] * sensor.parameters[final_tif_calc][sensor.WB];
        albedo_line[col] = (albedo_line[col] - 0.03) / (tal_line[col] * tal_line[col]);
    }

}; //alb

void ndvi_function(double reflectance_line[][8], int width_band, double ndvi_line[]){

    for (int col = 0; col < width_band; col++){
        ndvi_line[col] = (reflectance_line[col][4] - reflectance_line[col][3]) /
                         (reflectance_line[col][4] + reflectance_line[col][3]);
    }

}; //ndvi

void lai_function(double reflectance_line[][8], int width_band, double lai_line[]){
    double savi_line[width_band];
    double L = 0.05;

    for (int col = 0; col < width_band; col++){

        savi_line[col] = ((1 + L) * (reflectance_line[col][4] - reflectance_line[col][3])) /
                         (L + (reflectance_line[col][4] + reflectance_line[col][3]));
        
        if (!isnan(savi_line[col]) && savi_line[col] > 0.687)
            lai_line[col] = 6;
        else if (!isnan(savi_line[col]) && savi_line[col] < 0.1)
            lai_line[col] = 0;
        else
            lai_line[col] = -log((0.69 - savi_line[col]) / 0.59) / 0.91;
        
    }

}; //lai

void evi_function(double reflectance_line[][8], int width_band, double evi_line[]){

    for (int col = 0; col < width_band; col++){
        evi_line[col] = 2.5 * ((reflectance_line[col][4] - reflectance_line[col][3])/
                        (reflectance_line[col][4] + 6 * reflectance_line[col][3] - 7.5 * reflectance_line[col][1] + 1));
    }

}; //evi

void enb_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double enb_emissivity_line[]){

    for(int col = 0; col < width_band; col++){
        if(ndvi_line[col] < 0 || lai_line[col] > 2.99)
            enb_emissivity_line[col] = 0.98;
        else
            enb_emissivity_line[col] = 0.97 + 0.0033 * lai_line[col];
    }

    
}; //enb

void eo_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double eo_emissivity_line[]){

    for(int col = 0; col < width_band; col++){
        if(ndvi_line[col] < 0 || lai_line[col] > 2.99)
            eo_emissivity_line[col] = 0.98;
        else
            eo_emissivity_line[col] = 0.95 + 0.01 * lai_line[col];
    }

}; //eo

void ea_emissivity_function(double tal_line[], int width_band, double ea_emissivity_line[]){

    for (int col = 0; col < width_band; col++)
        ea_emissivity_line[col] = 0.85 * pow((-1 * log(tal_line[col])), 0.09);

}; //ea

void surface_temperature_function(double radiance_line[][8], double enb_emissivity_line[], int number_sensor, int width_band, double surface_temperature_line[]){
    double k1, k2;

    switch(number_sensor){
        case 5:
            k1 = 607.76;
            k2 = 1282.71;           

            break;
        case 7:
            k1 = 666.09;
            k2 = 1260.56;
        
            break;
        case 8:
            k1 = 774.8853;
            k2 = 1321.0789;

            break;
        default:
            cerr << "Sensor problem!";
            exit(6);
    }

    int radiance_number = (number_sensor == 5)? 6: 7;

    for(int col = 0; col < width_band; col++) 
        surface_temperature_line[col] = k2 / (log( (enb_emissivity_line[col] * k1 / radiance_line[col][radiance_number]) + 1));
    

}; //Ts

void short_wave_radiation_function(double tal_line[], MTL mtl, int width_band, double short_wave_radiation_line[]){
    double costheta = sin(mtl.sun_elevation * PI / 180);

    for (int col = 0; col < width_band; col++){
        short_wave_radiation_line[col] = (1367 * costheta * tal_line[col]) /
                                        (mtl.distance_earth_sun * mtl.distance_earth_sun);
    }
}; //Rs

void large_wave_radiation_surface_function(double eo_emissivity_line[], double surface_temperature_line[], int width_band, double large_wave_radiation_surface_line[]){

    for(int col = 0; col < width_band; col++){
        double temperature_pixel = surface_temperature_line[col];
        double surface_temperature_pow_4 = temperature_pixel * temperature_pixel * temperature_pixel * temperature_pixel;
        large_wave_radiation_surface_line[col] = eo_emissivity_line[col] * 5.67 * 1e-8 * surface_temperature_pow_4;
    }

}; //RLsup

void large_wave_radiation_atmosphere_function(double ea_emissivity_line[], int width_band, double temperature, double large_wave_radiation_atmosphere_line[]){

    double temperature_kelvin = temperature + 273.15;
    double temperature_kelvin_pow_4 = temperature_kelvin * temperature_kelvin * temperature_kelvin * temperature_kelvin;

    for(int col = 0; col < width_band; col++)
        large_wave_radiation_atmosphere_line[col] = ea_emissivity_line[col] * 5.67 * 1e-8 * temperature_kelvin_pow_4;

}; //RLatm

void net_radiation_function(double short_wave_radiation_line[], double large_wave_radiation_surface_line[],
                            double large_wave_radiation_atmosphere_line[], double albedo_line[],
                            double eo_emissivity_line[], int width_band, double net_radiation_line[]){

    for (int col = 0; col < width_band; col++){
        net_radiation_line[col] = short_wave_radiation_line[col] - (short_wave_radiation_line[col] * albedo_line[col]) +
                                large_wave_radiation_atmosphere_line[col] - large_wave_radiation_surface_line[col] -
                                (1 - eo_emissivity_line[col]) * large_wave_radiation_atmosphere_line[col];

        if(net_radiation_line[col] < 0)
            net_radiation_line[col] = 0;
    }

}; //Rn

void soil_heat_flux_function(double ndvi_line[], double surface_temperature_line[], double albedo_line[], double net_radiation_line[], int width_band, double soil_heat_line[]){

    for(int col = 0; col < width_band; col ++){
        if(ndvi_line[col] >= 0){
            double ndvi_pixel_pow_4 = ndvi_line[col] * ndvi_line[col] * ndvi_line[col] * ndvi_line[col];
            soil_heat_line[col] = (surface_temperature_line[col] - 273.15) * (0.0038 + 0.0074 * albedo_line[col]) *
                                (1 - 0.98 * ndvi_pixel_pow_4) * net_radiation_line[col];
        }else
            soil_heat_line[col] = 0.5 * net_radiation_line[col];
        
        if(soil_heat_line[col] < 0)
            soil_heat_line[col] = 0;

    }

}; //G

void ho_fuction(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]){

    for(int col = 0; col < width_band; col++)
        ho_line[col] = net_radiation_line[col] - soil_heat_flux[col];

}; //HO

Candidate select_hot_pixel(TIFF* ndvi, TIFF* surface_temperature, TIFF* net_radiation, TIFF* soil_heat, int heigth_band, int width_band){

    double ndvi_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band], soil_heat_line[width_band];
    double ho_line[width_band];

    vector<Candidate> pre_candidates;

    for(int line = 0; line < heigth_band; line ++){

        read_line_tiff(net_radiation, net_radiation_line, line);
        read_line_tiff(soil_heat, soil_heat_line, line);

        ho_fuction(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(ndvi, ndvi_line, line);
        read_line_tiff(surface_temperature, surface_temperature_line, line);

        for(int col = 0; col < width_band; col ++){
            if(!isnan(ndvi_line[col]) && ndvi_line[col] > 0.15 && ndvi_line[col] < 0.20 && surface_temperature_line[col] > 273.16){
                pre_candidates.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col]));
            }
        }
    }
    
    sort(pre_candidates.begin(), pre_candidates.end()); // ARRUMAR //ORDENAR POR TEMPERATURA
    int pos = round(0.95 * pre_candidates.size());
    double surface_temperature_hot_pixel = pre_candidates[pos].temperature;

    vector<Candidate> candidates;
    for(Candidate c : pre_candidates){
        if(c.temperature == surface_temperature_hot_pixel)
            candidates.push_back(c);
    }

    if(candidates.size() > 1){
        
        vector< pair<double, int> > ndvi_hot;

        sort(ho_c_hot.begin(), ho_c_hot.end());
        int posmin = round(0.25 * ho_c_hot.size()), posmax = round(0.75 * ho_c_hot.size());
        double ho_c_hot_min = ho_c_hot[posmin], ho_c_hot_max = ho_c_hot[posmax];

        for(int col = 0; col < width_band; col++){
            if(surface_temperature_line[col] == surface_temperature_hot_pixel && ho_line[col] > ho_c_hot_min && ho_line[col] < ho_c_hot_max){
                ndvi_hot.push_back(make_pair(ndvi_line[col], col));
            }
        }

        if(ndvi_hot.size() == 0) return;

        sort(ndvi_hot.begin(), ndvi_hot.end());
        hot_pixel_candidates.push_back(Candidate(line, ndvi_hot[0].second));
    }

    /*  TODO
        Provalmente tem como otimizar isso aqui, tlvz algumas das iterações seja
        desnecessária. Mas vou almoçar agr.
    */

};

Candidate select_cold_pixel(TIFF* ndvi, TIFF* surface_temperature, TIFF* net_radiation, TIFF* soil_heat, int heigth_band, int width_band){

};