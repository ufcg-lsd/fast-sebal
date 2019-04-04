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
    this->zom_path = output_path + "/zom.tif";
    this->ustar_path = output_path + "/ustar.tif";
    this->aerodynamic_resistence_path = output_path + "/Rah.tif";
    this->sensible_heat_flux_path = output_path + "/H.tif";
    this->ustar_tif0_path = output_path + "/ustar_tif0.tif";
    this->ustar_tif1_path = output_path + "/ustar_tif1.tif";
    this->aerodynamic_resistence_tif0_path = output_path + "/Rah_tif0.tif";
    this->aerodynamic_resistence_tif1_path = output_path + "/Rah_tif1.tif";
    this->latent_heat_flux_path = output_path + "/LatentHF.tif";
    this->net_radiation_24h_path = output_path + "/Rn24h.tif";
    this->latent_heat_flux_24h_path = output_path + "/LatentHF24h.tif";
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

    //hot_pixel.setAerodynamicResistance(u200, station.A_ZOM, station.B_ZOM, VON_KARMAN);
    //cold_pixel.setAerodynamicResistance(u200, station.A_ZOM, station.B_ZOM, VON_KARMAN);
    hot_pixel.toString();
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
    
    //FIXME: auxiliar products TIFFs
    TIFF *zom, *ustar, *aerodynamic_resistence;
    zom = TIFFOpen(zom_path.c_str(), "w8m");
    setup(zom, albedo);

    ustar = TIFFOpen(ustar_path.c_str(), "w8m");
    setup(ustar, albedo);

    aerodynamic_resistence = TIFFOpen(aerodynamic_resistence_path.c_str(), "w8m");
    setup(aerodynamic_resistence, albedo);

    for(int line = 0; line < heigth_band; line++){
        read_line_tiff(ndvi, ndvi_line, line);
        read_line_tiff(surface_temperature, surface_temperature_line, line);
        read_line_tiff(net_radiation, net_radiation_line, line);
        read_line_tiff(soil_heat, soil_heat_line, line);
        read_line_tiff(albedo, albedo_line, line);

        zom_fuction(station.A_ZOM, station.B_ZOM, ndvi_line, width_band, zom_line); 
        ustar_fuction(u200, zom_line, width_band, ustar_line);
        aerodynamic_resistence_fuction(ustar_line, width_band, aerodynamic_resistence_line);

        save_tiffs(vector<double*> {zom_line, ustar_line, aerodynamic_resistence_line}, 
               vector<TIFF*> {zom, ustar, aerodynamic_resistence}, line);
    }

    //Initial zom, ustar and aerodynamic_resistence are calculated and saved.
    //Continuing the sebal calculation

    //Extract the hot pixel aerodynamic_resistance
    hot_pixel.aerodynamic_resistance.push_back(read_position_tiff(aerodynamic_resistence, hot_pixel.col, hot_pixel.line));
    double H_hot = hot_pixel.net_radiation - hot_pixel.soil_heat_flux;

    TIFFClose(ndvi);
    TIFFClose(zom);
    TIFFClose(ustar);
    TIFFClose(aerodynamic_resistence);

    TIFF *ustar_tif0, *ustar_tif1, *aerodynamic_resistence_tif0, *aerodynamic_resistence_tif1, *sensible_heat_flux;
    zom = TIFFOpen(zom_path.c_str(), "rm"); //It's not modified into the rah cycle

    //It's only written into the rah cycle
    sensible_heat_flux = TIFFOpen(sensible_heat_flux_path.c_str(), "w8m");
    setup(sensible_heat_flux, albedo);

    int i = 0;
    bool Erro = true;
    
    //Auxiliar arrays calculation      
    double L[width_band];
    double y_01_line[width_band], y_2_line[width_band], x_200_line[width_band];
    double psi_01_line[width_band], psi_2_line[width_band], psi_200_line[width_band];

    double rah_hot0;
    double rah_hot;

    while(Erro) {
        cout << i << endl;
        if(i%2) {
            //Since ustar is both write and read into the rah cycle, two TIFF will be needed
            ustar_tif0 = TIFFOpen(ustar_tif1_path.c_str(), "rm");
            ustar_tif1 = TIFFOpen(ustar_path.c_str(), "w8m");
            setup(ustar_tif1, albedo);

            //Since ustar is both write and read into the rah cycle, two TIFF will be needed
            aerodynamic_resistence_tif0 = TIFFOpen(aerodynamic_resistence_tif1_path.c_str(), "rm");
            aerodynamic_resistence_tif1 = TIFFOpen(aerodynamic_resistence_path.c_str(), "w8m");
            setup(aerodynamic_resistence_tif1, albedo);
        } else {
            //Since ustar is both write and read into the rah cycle, two TIFF will be needed
            ustar_tif0 = TIFFOpen(ustar_path.c_str(), "rm");
            ustar_tif1 = TIFFOpen(ustar_tif1_path.c_str(), "w8m");
            setup(ustar_tif1, albedo);

            //Since ustar is both write and read into the rah cycle, two TIFF will be needed
            aerodynamic_resistence_tif0 = TIFFOpen(aerodynamic_resistence_path.c_str(), "rm");
            aerodynamic_resistence_tif1 = TIFFOpen(aerodynamic_resistence_tif1_path.c_str(), "w8m");
            setup(aerodynamic_resistence_tif1, albedo);
        }

        for(int line = 0; line < heigth_band; line++){

            //Reading data needed
            read_line_tiff(surface_temperature, surface_temperature_line, line);
            read_line_tiff(zom, zom_line, line);
            read_line_tiff(ustar_tif0, ustar_line, line);
            read_line_tiff(aerodynamic_resistence_tif0, aerodynamic_resistence_line, line);

            double dt_hot = H_hot * rah_hot0 / (RHO * SPECIFIC_HEAT_AIR);
            double b = dt_hot/(hot_pixel.temperature - cold_pixel.temperature);
            double a = -b * (cold_pixel.temperature - 273.15);

            for(int col = 0; col < width_band; col++) {
                sensible_heat_flux_line[col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_line[col] - 273.15))/aerodynamic_resistence_line[col];
                cout << sensible_heat_flux_line[col] << endl;
                double ustar_pow_3 = ustar_line[col] * ustar_line[col] * ustar_line[col];
                L[col] = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * surface_temperature_line[col])/(VON_KARMAN * GRAVITY * sensible_heat_flux_line[col]));

                y_01_line[col] = pow((1 - (16*0.1)/L[col]), 0.25);
                y_2_line[col] = pow((1 - (16*2)/L[col]), 0.25);
                x_200_line[col] = pow((1 - (16*200)/L[col]), 0.25);

                if(!isnan(L[col]) && L[col] > 0) psi_01_line[col] = -5 * (0.1/L[col]);
                else psi_01_line[col] = 2 * log((1 + y_01_line[col]*y_01_line[col])/2);

                if(!isnan(L[col]) && L[col] > 0) psi_2_line[col] = -5 * (2/L[col]);
                else psi_2_line[col] = 2 * log((1 + y_2_line[col]*y_2_line[col])/2);

                if(!isnan(L[col]) && L > 0) psi_200_line[col] = -5 * (2/L[col]);
                else psi_200_line[col] = 2 * log((1 + x_200_line[col])/2) + log((1 + x_200_line[col]*x_200_line[col])/2) - 2 * atan(x_200_line[col]) + 0.5 * PI;

                ustar_line[col] = (VON_KARMAN * u200) / (log(200/zom_line[col]) - psi_200_line[col]);
                aerodynamic_resistence_line[col] = (log(2/0.1) - psi_2_line[col] + psi_01_line[col])/(ustar_line[col] * VON_KARMAN);

                if(line == hot_pixel.line && col == hot_pixel.col) {
                    cout << "L " << L[col] << endl;
                    cout << "y_line " << y_01_line[col] << " " << y_2_line[col] << " " << x_200_line[col] << endl;
                    cout << "psi-line " << psi_01_line[col] << " " << psi_2_line[col] << " " << psi_200_line[col] << endl;
                    cout << ustar_line[col] << endl;
                    cout << aerodynamic_resistence_line[col] << endl;

                    rah_hot = aerodynamic_resistence_line[col];
                    hot_pixel.aerodynamic_resistance.push_back(rah_hot);
                }
            
            }

            //Saving new ustar e rah
            save_tiffs(vector<double*> {ustar_line, aerodynamic_resistence_line}, 
                    vector<TIFF*> {ustar_tif1, aerodynamic_resistence_tif1}, line);

        }

        TIFFClose(ustar_tif0);
        TIFFClose(ustar_tif1);
        TIFFClose(aerodynamic_resistence_tif0);
        TIFFClose(aerodynamic_resistence_tif1);

        cout << fabs(1 - rah_hot0/rah_hot) << endl;
        Erro = (fabs(1 - rah_hot0/rah_hot) >= 0.05);
        i++;

    }

    TIFFClose(zom);
    TIFFClose(ustar);
    TIFFClose(aerodynamic_resistence);
    TIFFClose(surface_temperature);

    if(i%2) {
        printf("Rah_after is aerodynamic_resistence_tif1_path\n");
        aerodynamic_resistence_tif1 = TIFFOpen(aerodynamic_resistence_tif1_path.c_str(), "rm");
    } else {
        printf("Rah_after is aerodynamic_resistence_path\n");
        aerodynamic_resistence_tif1 = TIFFOpen(aerodynamic_resistence_path.c_str(), "rm");
    }

    double dt_hot = H_hot * rah_hot / (RHO * SPECIFIC_HEAT_AIR);
    double b = dt_hot/(hot_pixel.temperature - cold_pixel.temperature);
    double a = -b * (cold_pixel.temperature - 273.15);

    for(int line = 0; line < heigth_band; line++){
        read_line_tiff(aerodynamic_resistence_tif1, aerodynamic_resistence_line, line);
            
        for(int col = 0; col < width_band; col++) {
            sensible_heat_flux_line[col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_line[col] - 273.15))/aerodynamic_resistence_line[col];
        }

        save_tiffs(vector<double*> {sensible_heat_flux_line}, 
                    vector<TIFF*> {sensible_heat_flux}, line);
    }

    TIFFClose(sensible_heat_flux);

    //End of Rah correction

    //Continuing to the final products

    TIFF *latent_heat_flux, *net_radiation_24h, *latent_heat_flux_24h;

    latent_heat_flux = TIFFOpen(latent_heat_flux_path.c_str(), "w8m");
    setup(latent_heat_flux, albedo);

    net_radiation_24h = TIFFOpen(net_radiation_24h_path.c_str(), "w8m");
    setup(net_radiation_24h, albedo);

    latent_heat_flux_24h = TIFFOpen(latent_heat_flux_24h_path.c_str(), "w8m");
    setup(latent_heat_flux_24h, albedo);

    sensible_heat_flux = TIFFOpen(sensible_heat_flux_path.c_str(), "rm");

    for(int line = 0; line < heigth_band; line++){
        read_line_tiff(net_radiation, net_radiation_line, line);
        read_line_tiff(soil_heat, soil_heat_line, line);
        read_line_tiff(albedo, albedo_line, line);
        read_line_tiff(sensible_heat_flux, sensible_heat_flux_line, line);

        latent_heat_flux_function(net_radiation_line, soil_heat_line, sensible_heat_flux_line, width_band, latent_heat_flux_line);
        net_radiation_24h_function(albedo_line, Ra24h, Rs24h, width_band, net_radiation_24h_line);
        evapotranspiration_fraction_fuction(latent_heat_flux_line, net_radiation_line, soil_heat_line, width_band, evapotranspiration_fraction_line);
        sensible_heat_flux_24h_fuction(evapotranspiration_fraction_line, net_radiation_24h_line, width_band, sensible_heat_flux_24h_line);
        latent_heat_flux_24h_function(evapotranspiration_fraction_line, net_radiation_24h_line, width_band, latent_heat_flux_24h_line);
        evapotranspiration_24h_function(latent_heat_flux_24h_line, station, width_band, evapotranspiration_24h_line);
        
        save_tiffs(vector<double*> {latent_heat_flux_line, net_radiation_24h_line, latent_heat_flux_24h_line, evapotranspiration_fraction_line, evapotranspiration_24h_line}, 
               vector<TIFF*> {latent_heat_flux, net_radiation_24h, latent_heat_flux_24h, evapotranspiration_fraction, evapotranspiration_24h}, line);
    
    }
    
    TIFFClose(albedo);
    TIFFClose(soil_heat);
    TIFFClose(net_radiation);
    TIFFClose(sensible_heat_flux);
    TIFFClose(latent_heat_flux);
    TIFFClose(net_radiation_24h);
    TIFFClose(latent_heat_flux_24h);
    TIFFClose(evapotranspiration_fraction);
    TIFFClose(evapotranspiration_24h);

    /*
    TIFFClose(albedo);
    TIFFClose(ndvi);
    TIFFClose(soil_heat);
    TIFFClose(surface_temperature);
    TIFFClose(net_radiation);
    TIFFClose(evapotranspiration_fraction);
    TIFFClose(evapotranspiration_24h);
    */

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
