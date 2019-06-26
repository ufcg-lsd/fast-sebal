#include "products.h"

string tal_function(TIFF *raster_elevation, string output_path){
    uint32 height_band, width_band;
    uint16 sample_band;

    TIFFGetField(raster_elevation, TIFFTAG_IMAGEWIDTH, &width_band);
    TIFFGetField(raster_elevation, TIFFTAG_IMAGELENGTH, &height_band);
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

    for (int line = 0; line < height_band; line++){
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
        for (int col = 0; col < width_band; col++) {
            radiance_line[col][7] = line_band[col] * mtl.rad_mult_10 + mtl.rad_add_10;
        }
    }
    else{
        for (int i = 1; i < 8; i++){
            read_line_tiff(read_bands[i], line_band, line);
            for (int col = 0; col < width_band; col++)
                radiance_line[col][i] = min(line_band[col] * sensor.parameters[i][sensor.GRESCALE] + sensor.parameters[i][sensor.BRESCALE], 0.0);
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
        
        if (!isnan(savi_line[col]) && definitelyGreaterThan(savi_line[col], 0.687))
            lai_line[col] = 6;
        else if (!isnan(savi_line[col]) && definitelyLessThan(savi_line[col], 0.1))
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
        if(definitelyLessThan(ndvi_line[col], 0) || definitelyGreaterThan(lai_line[col], 2.99))
            enb_emissivity_line[col] = 0.98;
        else
            enb_emissivity_line[col] = 0.97 + 0.0033 * lai_line[col];
    }

    
}; //enb

void eo_emissivity_function(double lai_line[], double ndvi_line[], int width_band, double eo_emissivity_line[]){

    for(int col = 0; col < width_band; col++){
        if(definitelyLessThan(ndvi_line[col], 0) || definitelyGreaterThan(lai_line[col], 2.99))
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

        if(definitelyLessThan(net_radiation_line[col], 0))
            net_radiation_line[col] = 0;
    }

}; //Rn

void soil_heat_flux_function(double ndvi_line[], double surface_temperature_line[], double albedo_line[], double net_radiation_line[], int width_band, double soil_heat_line[]){

    for(int col = 0; col < width_band; col ++){
        if(essentiallyEqual(ndvi_line[col], 0) || definitelyGreaterThan(ndvi_line[col], 0)){
            double ndvi_pixel_pow_4 = ndvi_line[col] * ndvi_line[col] * ndvi_line[col] * ndvi_line[col];
            soil_heat_line[col] = (surface_temperature_line[col] - 273.15) * (0.0038 + 0.0074 * albedo_line[col]) *
                                (1 - 0.98 * ndvi_pixel_pow_4) * net_radiation_line[col];
        }else
            soil_heat_line[col] = 0.5 * net_radiation_line[col];
        
        if(definitelyLessThan(soil_heat_line[col], 0))
            soil_heat_line[col] = 0;

    }

}; //G

void ho_function(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]){

    for(int col = 0; col < width_band; col++)
        ho_line[col] = net_radiation_line[col] - soil_heat_flux[col];

}; //HO

Candidate select_hot_pixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band){
    
    //Auxiliary arrays
    double ndvi_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band], soil_heat_line[width_band];
    double ho_line[width_band];

    //Contains the candidates with NDVI between 0.15 and 0.20, which surface temperature is greater than 273.16
    vector<Candidate> pre_candidates;

    for(int line = 0; line < height_band; line ++){
        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        ho_function(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);

        for(int col = 0; col < width_band; col ++){
            if(!isnan(ndvi_line[col]) && definitelyGreaterThan(ndvi_line[col], 0.15) && definitelyLessThan(ndvi_line[col], 0.20) && definitelyGreaterThan(surface_temperature_line[col], 273.16)){
                pre_candidates.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col],
                                    ho_line[col],
                                    line, col));
            }
        }

    }

    //Sort the candidates by their temperatures and choose the surface temperature of the hot pixel
    sort(pre_candidates.begin(), pre_candidates.end(), compare_candidate_temperature);
    int pos = floor(0.95 * pre_candidates.size());
    double surfaceTempHot = pre_candidates[pos].temperature;

    //Select only the ones with temperature equals the surface temperature of the hot pixel
    vector<double> ho_candidates;
    Candidate lastHOCandidate;
    for(Candidate c : pre_candidates){
        if(essentiallyEqual(c.temperature, surfaceTempHot)){
            ho_candidates.push_back(c.ho);
            lastHOCandidate = c;
        }
    }

    if(ho_candidates.size() == 1){
        return lastHOCandidate;
    }

    //Select the limits of HOs
    sort(ho_candidates.begin(), ho_candidates.end());
    double HO_min = ho_candidates[floor(0.25 * ho_candidates.size())];
    double HO_max = ho_candidates[floor(0.75 * ho_candidates.size())];

    //Contains the final candidates which HO is in (HO_min, HO_max) and surface temperature is greater than 273.16
    vector<Candidate> final_candidates;

    for(int line = 0; line < height_band; line ++){

        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        ho_function(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);

        for(int col = 0; col < width_band; col ++){
            if(definitelyGreaterThan(ho_line[col], HO_min) && definitelyLessThan(ho_line[col], HO_max) && essentiallyEqual(surface_temperature_line[col], surfaceTempHot)){
                final_candidates.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col],
                                    ho_line[col],
                                    line, col));
            }
        }

    }
    
    //Calculate the coefficient of variation, after the extract
    for(int i = 0; i < final_candidates.size(); i++){
        final_candidates[i].extract_coefficient_variation(*ndvi);
    }

    //Choose as candidate the pixel with the minor CV
    Candidate choosen = final_candidates[0];

    for(int i = 1; i < final_candidates.size(); i++){
        
        if(definitelyLessThan(final_candidates[i].coefficient_variation, choosen.coefficient_variation))
            choosen = final_candidates[i];

    }

    return choosen;
}

Candidate select_cold_pixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band){
    
    //Auxiliary arrays
    double ndvi_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band], soil_heat_line[width_band];
    double ho_line[width_band];

    //Contains the candidates with NDVI less than 0, which surface temperature is greater than 273.16
    vector<Candidate> pre_candidates;

    for(int line = 0; line < height_band; line ++){

        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        ho_function(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);

        for(int col = 0; col < width_band; col ++){
            if(!isnan(ndvi_line[col]) && !isnan(ho_line[col]) && definitelyLessThan(ndvi_line[col], 0) && definitelyGreaterThan(surface_temperature_line[col], 273.16)){
                pre_candidates.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col],
                                    ho_line[col],
                                    line, col));
            }
        }

    }

    //Sort the candidates by their temperatures and choose the surface temperature of the hot pixel
    sort(pre_candidates.begin(), pre_candidates.end(), compare_candidate_temperature);
    int pos = floor(0.5 * pre_candidates.size());
    double surfaceTempCold = pre_candidates[pos].temperature;

    //Select only the ones with temperature equals the surface temperature of the Cold pixel
    vector<double> ho_candidates;
    Candidate lastHOCandidate;
    for(Candidate c : pre_candidates){
        if(essentiallyEqual(c.temperature, surfaceTempCold)){
            ho_candidates.push_back(c.ho);
            lastHOCandidate = c;
        }
    }

    if(ho_candidates.size() == 1){
        return lastHOCandidate;
    }

    //Select the limits of HOs
    sort(ho_candidates.begin(), ho_candidates.end());
    double HO_min = ho_candidates[floor(0.25 * ho_candidates.size())];
    double HO_max = ho_candidates[floor(0.75 * ho_candidates.size())];

    //Contains the final candidates which HO is in (HO_min, HO_max) and surface temperature is greater than 273.16
    vector<Candidate> final_candidates;

    for(int line = 0; line < height_band; line ++){

        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        ho_function(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);

        for(int col = 0; col < width_band; col ++){
            if(definitelyGreaterThan(ho_line[col], HO_min) && definitelyLessThan(ho_line[col], HO_max) && essentiallyEqual(surface_temperature_line[col], surfaceTempCold)){
                final_candidates.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col],
                                    ho_line[col],
                                    line, col));
            }
        }

    }
    
    //Calculate the coefficient of variation, after the extract
    for(int i = 0; i < final_candidates.size(); i++){
        final_candidates[i].extract_negative_neighbour(*ndvi);
    }

    //Choose as candidate the pixel with the minor CV
    Candidate choosen = final_candidates[0];
    for(int i = 1; i < final_candidates.size(); i++){
        if(final_candidates[i].negative_neighbour > choosen.negative_neighbour)
            choosen = final_candidates[i];
    }

    return choosen;
}

void zom_fuction(double A_ZOM, double B_ZOM, double ndvi_line[], int width_band, double zom_line[]){

    for(int col = 0; col < width_band; col++)
        zom_line[col] = exp(A_ZOM + B_ZOM * ndvi_line[col]);

}; //zom

void ustar_fuction(double u200, double zom_line[], int width_band, double ustar_line[]){

    for(int col = 0; col < width_band; col++)
        ustar_line[col] = (VON_KARMAN * u200)/log(200/zom_line[col]);

}; //ustar

void aerodynamic_resistance_fuction(double ustar_line[], int width_band, double aerodynamic_resistance_line[]){

    for(int col = 0; col < width_band; col++)
        aerodynamic_resistance_line[col] = log(20)/(ustar_line[col] * VON_KARMAN);

}; //rah

void latent_heat_flux_function(double net_radiation_line[], double soil_heat_flux_line[], double sensible_heat_flux_line[], int width_band, double latent_heat_flux[]){

    for(int col = 0; col < width_band; col++)
        latent_heat_flux[col] = net_radiation_line[col] - soil_heat_flux_line[col] - sensible_heat_flux_line[col];

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
