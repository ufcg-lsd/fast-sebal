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

void radiance_function(TIFF* read_bands[], MTL mtl, Sensor sensor, int width_band, int line, double radiance_line[][8]){
    double line_band[width_band];

    if (mtl.number_sensor == 8){
        read_line_tiff(read_bands[7], line_band, line);
        for (int col = 0; col < width_band; col++)
            radiance_line[col][7] = line_band[col] * mtl.rad_mult_10 + mtl.rad_add_10;
    }
    else{
        for (int i = 1; i < 8; i++){
            read_line_tiff(read_bands[i], line_band, line);
            for (int col = 0; col < width_band; col++)
                radiance_line[col][i] = line_band[col] * sensor.parameters[i][sensor.GRESCALE] + sensor.parameters[i][sensor.BRESCALE];
        }
    }
}; //rad

void reflectance_function(TIFF* read_bands[], MTL mtl, Sensor sensor, double radiance_line[][8], int width_band, int line, double reflectance_line[][8]){
    double costheta = sin(mtl.sun_elevation * PI / 180);
    double line_band[width_band];

    for (int i = 1; i < 8; i++){
        read_line_tiff(read_bands[i], line_band, line);
        for (int col = 0; col < width_band; col++){
            if (mtl.number_sensor == 8)
                reflectance_line[col][i] = (line_band[col] * sensor.parameters[i][sensor.GRESCALE] + sensor.parameters[i][sensor.BRESCALE]) / costheta;
            else
                reflectance_line[col][i] = (PI * radiance_line[col][i] * mtl.distance_earth_sun * mtl.distance_earth_sun) /
                                           (sensor.parameters[i][sensor.ESUN] * costheta);
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
                                    soil_heat_line[col],
                                    ho_line[col]));
            }
        }
    }
    
    sort(pre_candidates.begin(), pre_candidates.end(), compare_candidate_temperature);
    int pos = round(0.95 * pre_candidates.size());
    double surface_temperature_hot_pixel = pre_candidates[pos].temperature;

    vector<Candidate> candidates;
    for(Candidate c : pre_candidates){
        if(c.temperature == surface_temperature_hot_pixel)
            candidates.push_back(c);
    }

    Candidate choosen;
    if(candidates.size() == 1){
        choosen = candidates[0];
    } else {
        vector< pair<double, int> > ndvi_hot;
        //TODO Verificar uso de floor ou round
        sort(candidates.begin(), candidates.end(), compare_candidate_ho); 
        int posmin = floor(0.25 * candidates.size()), posmax = floor(0.75 * candidates.size());
       
        choosen = candidates[posmin+1];
        for(int i = posmin+2; i < posmax; i++)
            if(candidates[i].ndvi < choosen.ndvi) choosen = candidates[i];
    }

    return choosen;
};

Candidate select_cold_pixel(TIFF* ndvi, TIFF* surface_temperature, TIFF* net_radiation, TIFF* soil_heat, int heigth_band, int width_band){
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
            if(!isnan(ndvi_line[col]) && !isnan(ho_line[col]) && ndvi_line[col] < 0 && surface_temperature_line[col] > 273.16){
                pre_candidates.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col],
                                    ho_line[col]));
            }
        }
    }
    
    sort(pre_candidates.begin(), pre_candidates.end(), compare_candidate_temperature);
    int pos = round(0.5 * pre_candidates.size());
    double surface_temperature_hot_pixel = pre_candidates[pos].temperature;

    vector<Candidate> candidates;
    for(Candidate c : pre_candidates){
        if(c.temperature == surface_temperature_hot_pixel)
            candidates.push_back(c);
    }

    Candidate choosen;
    if(candidates.size() == 1){
        choosen = candidates[0];
    } else {
        vector< pair<double, int> > ndvi_hot;
        //TODO Verificar uso de floor ou round
        sort(candidates.begin(), candidates.end(), compare_candidate_ho); 
        int posmin = floor(0.25 * candidates.size()), posmax = floor(0.75 * candidates.size());
       
        choosen = candidates[posmin+1];
        for(int i = posmin+2; i < posmax; i++)
            if(candidates[i].ndvi > choosen.ndvi) choosen = candidates[i];
    }

    return choosen;
};

void zom_fuction(double A_ZOM, double B_ZOM, double ndvi_line[], int width_band, double zom_line[]){

    for(int col = 0; col < width_band; col++)
        zom_line[col] = A_ZOM + B_ZOM * ndvi_line[col];

}; //zom

void ustar_fuction(double u200, double zom_line[], int width_band, double ustar_line[]){

    for(int col = 0; col < width_band; col++)
        ustar_line[col] = (VON_KARMAN * u200)/log(200/zom_line[col]);

}; //ustar

void aerodynamic_resistence_fuction(double ustar_line[], int width_band, double aerodynamic_resistence_line[]){

    for(int col = 0; col < width_band; col++)
        aerodynamic_resistence_line[col] = log(2/0.1)/(ustar_line[col] * VON_KARMAN);

}; //rah

void sensible_heat_flux_function(Candidate hot_pixel, Candidate cold_pixel, double u200, double zom_line[], double ustar_line[], double aerodynamic_resistence_line[], double surface_temperature_line[], int width_band, double sensible_heat_flux_line[]){
    double H_hot = hot_pixel.net_radiation - hot_pixel.soil_heat_flux;
    double rah_hot0;

    //LINE PIXEL CALCULATION
    double L[width_band];
    double y_01_line[width_band], y_2_line[width_band], x_200_line[width_band];
    double psi_01_line[width_band], psi_2_line[width_band], psi_200_line[width_band];
    
    for(unsigned i = 0; i < hot_pixel.aerodynamic_resistance.size(); i++){
        rah_hot0 = hot_pixel.aerodynamic_resistance[i];

        double dt_hot = (H_hot * rah_hot0)/(RHO * SPECIFIC_HEAT_AIR);
        double b = dt_hot/(hot_pixel.temperature - cold_pixel.temperature);
        double a = -b * (cold_pixel.temperature - 273.15);

        for(int col = 0; col < width_band; col++){
            sensible_heat_flux_line[col] = (RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_line[col] - 273.15)))/aerodynamic_resistence_line[col];
            double ustar_pow_3 = ustar_line[col] * ustar_line[col] * ustar_line[col];
            L[col] = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * surface_temperature_line[col])/(VON_KARMAN * GRAVITY * sensible_heat_flux_line[col]));
            
            y_01_line[col] = pow((1 - 16*0.1/L[col]), 0.25);
            y_2_line[col] = pow((1 - 16*2/L[col]), 0.25);
            x_200_line[col] = pow((1 - 16*200/L[col]), 0.25);

            if(L[col] > 0) psi_01_line[col] = -5 * (0.1/L[col]);
            else psi_01_line[col] = 2 * log((1 + y_01_line[col]*y_01_line[col])/2);

            if(L[col] > 0) psi_2_line[col] = -5 * (2/L[col]);
            else psi_2_line[col] = 2 * log((1 + y_2_line[col]*y_2_line[col])/2);

            if(L > 0) psi_200_line[col] = -5 * (2/L[col]);
            else psi_200_line[col] = 2 * log((1 + x_200_line[col])/2) + log((1 + x_200_line[col]*x_200_line[col])/2) - 2 * atan(x_200_line[col]) + 0.5 * acos(-1);

            ustar_line[col] = (VON_KARMAN * u200) / (log(200/zom_line[col]) - psi_200_line[col]);
            aerodynamic_resistence_line[col] = (log(2/0.1) - psi_2_line[col] + psi_01_line[col])/(ustar_line[col] * VON_KARMAN);
        }
        
    }

    double dt_hot = (H_hot * rah_hot0)/(RHO * SPECIFIC_HEAT_AIR);
    double b = dt_hot/(hot_pixel.temperature - cold_pixel.temperature);
    double a = -b * (cold_pixel.temperature - 273.15);

    for(int col = 0; col < width_band; col ++)
        sensible_heat_flux_line[col] = (RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_line[col] - 273.15)))/aerodynamic_resistence_line[col];

}; //H

void latent_heat_flux_function(double net_radiation_line[], double soil_heat_flux_line[], double sensible_heat_flux_line[], int width_band, double latent_heat_flux[]){

    for(int col = 0; col < width_band; col++)
        latent_heat_flux[col] = max(net_radiation_line[col] - soil_heat_flux_line[col] - sensible_heat_flux_line[col], 0.0);

}; //LE

void net_radiation_24h_function(double albedo_line[], double Ra24h, double Rs24h, int width_band, double net_radiation_24h_line[]){
    int FL = 110;

    for(int col = 0; col < width_band; col++)
        net_radiation_24h_line[col] = (1 - albedo_line[col])*Rs24h - FL * Rs24h/Ra24h;

}; //Rn24h_dB

void evapotranspiration_fraction_fuction(double latent_heat_flux_line[], double net_radiation_line[], double soil_heat_line[], int width_band, double evapotranspiration_fraction_line[]){

    for(int col = 0; col < width_band; col++)
        evapotranspiration_fraction_line[col] = latent_heat_flux_line[col]/(net_radiation_line[col] - soil_heat_line[col]);

}; //EF

void sensible_heat_flux_24h_fuction(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double sensible_heat_flux_24h_line[]){

    for(int col = 0; col < width_band; col++)
        sensible_heat_flux_24h_line[col] = (1 - evapotranspiration_fraction_line[col]) * net_radiation_24h_line[col];
        
}; //H24h_dB

void latent_heat_flux_24h_function(double evapotranspiration_fraction_line[], double net_radiation_24h_line[], int width_band, double latent_heat_flux_24h_line[]){

    for(int col = 0; col < width_band; col++)
        latent_heat_flux_24h_line[col] = evapotranspiration_fraction_line[col] * net_radiation_24h_line[col];
        
}; //LE24h_db

void evapotranspiration_24h_function(double latent_heat_flux_24h_line[], Station station, int width_band, double evapotranspiration_24h_line[]){

    for(int col = 0; col < width_band; col++)
        evapotranspiration_24h_line[col] = (latent_heat_flux_24h_line[col] * 86400)/((2.501 - 0.00236 * (station.v7_max + station.v7_min) / 2) * 1e+6);

}; //ET24h_dB