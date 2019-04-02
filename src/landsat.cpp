#include "landsat.h"

Landsat::Landsat(){
};

Landsat::Landsat(string tal_path, string output_path){
    this->tal_path = tal_path;
    this->output_path = output_path;
    this->albedo_path = output_path + "/alb.tif";
    this->ndvi_path = output_path + "/NDVI.tif";
    this->evi_path = output_path + "/EVI.tif";
    this->lai_path = output_path + "/LAI.tif";
    this->soil_heat_path = output_path + "/G.tif";
    this->surface_temperature_path = output_path + "/TS.tif";
    this->net_radiation_path = output_path + "/Rn.tif";
    this->evapotranspiration_fraction_path = output_path + "/EF.tif";
    this->evapotranspiration_24h_path = output_path + "/ET24h.tif";
};

void Landsat::process_parcial_products(TIFF* read_bands[], MTL mtl, Station station, Sensor sensor){
    uint32 heigth_band, width_band;

    TIFFGetField(read_bands[1], TIFFTAG_IMAGEWIDTH, &width_band);
    TIFFGetField(read_bands[1], TIFFTAG_IMAGELENGTH, &heigth_band);

    TIFF *tal = TIFFOpen(this->tal_path.c_str(), "rm");
    check_open_tiff(tal);

    double tal_line[width_band];

    TIFF *albedo, *ndvi, *evi, *lai, *soil_heat, *surface_temperature, *net_radiation;
    create_tiffs(&tal, &albedo, &ndvi, &evi, &lai, &soil_heat, &surface_temperature, &net_radiation);

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
        radiance_function(read_bands, mtl, sensor, width_band, line, radiance_line);
        reflectance_function(read_bands, mtl, sensor, radiance_line, width_band, line, reflectance_line);

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

void Landsat::process_final_products(Station station, MTL mtl){
    TIFF *albedo, *ndvi, *soil_heat, *surface_temperature, *net_radiation;
    TIFF *evapotranspiration_fraction, *evapotranspiration_24h;

    open_tiffs(&albedo, &ndvi, &soil_heat, &surface_temperature, &net_radiation, &evapotranspiration_fraction, &evapotranspiration_24h);

    uint32 heigth_band, width_band;
    TIFFGetField(albedo, TIFFTAG_IMAGELENGTH, &heigth_band);
    TIFFGetField(albedo, TIFFTAG_IMAGEWIDTH, &width_band);

    Candidate hot_pixel = select_hot_pixel(&ndvi, &surface_temperature, &net_radiation, &soil_heat, heigth_band, width_band);
    Candidate cold_pixel = select_cold_pixel(&ndvi, &surface_temperature, &net_radiation, &soil_heat, heigth_band, width_band);
/*
    To run without selecting pixels.

    Candidate(double ndvi, double temperature, double net_radiation, double soil_heat_flux, double ho, int line, int col);

    Candidate hot_pixel = Candidate(0.227679967882490234375, 308.387664794921875, 476.35150146484375, 101.90081024169921875, 476.35150146484375 - 101.90081024169921875, 0, 0);
    Candidate cold_pixel = Candidate(-0.1008398681879043579102, 297.303009033203125, 766.78015136718750, 383.39007568359375, 766.780105136718750 - 383.39007568359375, 0, 0);
*/   
    double sensible_heat_flux_line[width_band];
    double zom_line[width_band];
    double ustar_line[width_band];
    double aerodynamic_resistence_line[width_band];
    double latent_heat_flux_line[width_band];

    double ustar_station = (VON_KARMAN * station.v6)/(log(station.WIND_SPEED/station.SURFACE_ROUGHNESS));
    double u200 = (ustar_station/VON_KARMAN) * log(200 / station.SURFACE_ROUGHNESS);

    hot_pixel.setAerodynamicResistance(u200, station.A_ZOM, station.B_ZOM, VON_KARMAN);
    cold_pixel.setAerodynamicResistance(u200, station.A_ZOM, station.B_ZOM, VON_KARMAN);
    hot_pixel.toString();
    cold_pixel.toString();

    double H_hot = hot_pixel.net_radiation - hot_pixel.soil_heat_flux;
    double value_pixel_rah = hot_pixel.aerodynamic_resistance[0];
    double rah_hot0;

    //HOT PIXEL AERODYNAMIC CALCULATION
    int i = 0;
    do{
        rah_hot0 = hot_pixel.aerodynamic_resistance[i];

        double dt_hot = (H_hot * rah_hot0)/(RHO * SPECIFIC_HEAT_AIR);
        double b = dt_hot/(hot_pixel.temperature - cold_pixel.temperature);
        double a = -b * (cold_pixel.temperature - 273.15);

        double sensible_heat_flux = (RHO * SPECIFIC_HEAT_AIR * (a + b * (hot_pixel.temperature - 273.15)))/rah_hot0;
        double ustar_pow_3 = hot_pixel.ustar * hot_pixel.ustar * hot_pixel.ustar;
        double L = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * hot_pixel.temperature)/(VON_KARMAN * GRAVITY * sensible_heat_flux));
        double y_01 = pow((1 - 16*0.1/L), 0.25);
        double y_2 = pow((1 - 16*2/L), 0.25);
        double x_200 = pow((1 - 16*200/L), 0.25);

        double psi_01;
        if(!isnan(L) && L > 0) psi_01 = -5 * (0.1 / L);
        else psi_01 = 2 * log((1 + y_01 * y_01)/2);

        double psi_2;
        if(!isnan(L) && L > 0) psi_2 = -5 * (2/L);
        else psi_2 = 2 * log((1 + y_2 * y_2)/2);

        double psi_200;
        if(!isnan(L) && L > 0) psi_200 = -5 * (2/L);
        else psi_200 = 2 * log((1 + x_200)/2) + log((1 + x_200*x_200)/2) - 2 * atan(x_200) + 0.5 * PI;

        hot_pixel.ustar = (VON_KARMAN * u200) / (log(200/hot_pixel.zom) - psi_200);
        hot_pixel.aerodynamic_resistance.push_back((log(20) - psi_2 + psi_01)/(hot_pixel.ustar * VON_KARMAN));
        i++;
    } while(fabs(1 - rah_hot0/hot_pixel.aerodynamic_resistance[i]) >= 0.05);

    printf("HOT PIXEL\n");
    hot_pixel.toString();

    printf("COLD PIXEL\n");
    cold_pixel.toString();

    //Parcial products
    double ndvi_line[width_band], surface_temperature_line[width_band];
    double soil_heat_line[width_band], net_radiation_line[width_band];
    double albedo_line[width_band];

    //Outhers products
    double net_radiation_24h_line[width_band];
    double evapotranspiration_fraction_line[width_band];
    double sensible_heat_flux_24h_line[width_band];
    double latent_heat_flux_24h_line[width_band];
    double evapotranspiration_24h_line[width_band];

    //Upscalling temporal
    double dr = (1 / mtl.distance_earth_sun) * (1 / mtl.distance_earth_sun);
    double sigma = 0.409*sin(((2*PI/365)*mtl.julian_day)-1.39);
    double phi = (PI/180) * station.latitude;
    double omegas = acos(-tan(phi)*tan(sigma));
    double Ra24h = (((24*60/PI)*GSC*dr)*(omegas*sin(phi)*
                                sin(sigma)+cos(phi)*cos(sigma)*sin(omegas)))*(1000000/86400.0);

    //Short wave radiation incident in 24 hours (Rs24h)
    double Rs24h = station.INTERNALIZATION_FACTOR * sqrt(station.v7_max - station.v7_min) * Ra24h;

    for(int line = 0; line < heigth_band; line++){
        read_line_tiff(ndvi, ndvi_line, line);
        read_line_tiff(surface_temperature, surface_temperature_line, line);
        read_line_tiff(net_radiation, net_radiation_line, line);
        read_line_tiff(soil_heat, soil_heat_line, line);
        read_line_tiff(albedo, albedo_line, line);

        zom_fuction(station.A_ZOM, station.B_ZOM, ndvi_line, width_band, zom_line); 
        ustar_fuction(u200, zom_line, width_band, ustar_line);
        aerodynamic_resistence_fuction(ustar_line, width_band, aerodynamic_resistence_line);
        sensible_heat_flux_function(hot_pixel, cold_pixel, u200, zom_line, ustar_line, aerodynamic_resistence_line, surface_temperature_line, width_band, sensible_heat_flux_line);
        latent_heat_flux_function(net_radiation_line, soil_heat_line, sensible_heat_flux_line, width_band, latent_heat_flux_line);
        net_radiation_24h_function(albedo_line, Ra24h, Rs24h, width_band, net_radiation_24h_line);
        evapotranspiration_fraction_fuction(latent_heat_flux_line, net_radiation_line, soil_heat_line, width_band, evapotranspiration_fraction_line);
        sensible_heat_flux_24h_fuction(evapotranspiration_fraction_line, net_radiation_24h_line, width_band, sensible_heat_flux_24h_line);
        latent_heat_flux_24h_function(evapotranspiration_fraction_line, net_radiation_24h_line, width_band, latent_heat_flux_24h_line);
        evapotranspiration_24h_function(latent_heat_flux_24h_line, station, width_band, evapotranspiration_24h_line);

        save_tiffs(vector<double*> {evapotranspiration_fraction_line, evapotranspiration_24h_line}, 
               vector<TIFF*> {evapotranspiration_fraction, evapotranspiration_24h}, line);
    }
    
    TIFFClose(albedo);
    TIFFClose(ndvi);
    TIFFClose(soil_heat);
    TIFFClose(surface_temperature);
    TIFFClose(net_radiation);
    TIFFClose(evapotranspiration_fraction);
    TIFFClose(evapotranspiration_24h);

};

void Landsat::create_tiffs(TIFF **tal, TIFF **albedo, TIFF **ndvi, TIFF **evi, TIFF **lai, TIFF **soil_heat, TIFF **surface_temperature, TIFF **net_radiation){
    *albedo = TIFFOpen(albedo_path.c_str(), "w8m");
    setup(*albedo, *tal);

    *ndvi = TIFFOpen(ndvi_path.c_str(), "w8m");
    setup(*ndvi, *tal);

    *evi = TIFFOpen(evi_path.c_str(), "w8m");
    setup(*evi, *tal);

    *lai = TIFFOpen(lai_path.c_str(), "w8m");
    setup(*lai, *tal);

    *soil_heat = TIFFOpen(soil_heat_path.c_str(), "w8m");
    setup(*soil_heat, *tal);

    *surface_temperature = TIFFOpen(surface_temperature_path.c_str(), "w8m");
    setup(*surface_temperature, *tal);

    *net_radiation = TIFFOpen(net_radiation_path.c_str(), "w8m");
    setup(*net_radiation, *tal);
};

void Landsat::open_tiffs(TIFF **albedo, TIFF **ndvi, TIFF **soil_heat, TIFF **surface_temperature, TIFF **net_radiation, TIFF **evapotranspiration_fraction, TIFF **evapotranspiration_24h){

    *albedo = TIFFOpen(albedo_path.c_str(), "rm");
    *ndvi = TIFFOpen(ndvi_path.c_str(), "rm");
    *soil_heat = TIFFOpen(soil_heat_path.c_str(), "rm");
    *surface_temperature = TIFFOpen(surface_temperature_path.c_str(), "rm");
    *net_radiation = TIFFOpen(net_radiation_path.c_str(), "rm");
    
    *evapotranspiration_fraction = TIFFOpen(evapotranspiration_fraction_path.c_str(), "w8m");
    setup(*evapotranspiration_fraction, *albedo);

    *evapotranspiration_24h = TIFFOpen(evapotranspiration_24h_path.c_str(), "w8m");
    setup(*evapotranspiration_24h, *albedo);

}

void Landsat::save_tiffs(vector<double*> products_line, vector<TIFF*> products, int line){

    for (unsigned i = 0; i < products.size(); i++)
        write_line_tiff(products[i], products_line[i], line);

};
