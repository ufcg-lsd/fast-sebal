#include "utils.h"
#include "landsat.h"
#include "products.h"

using namespace std;

/*
    Complete execution

    arg 01 - band 01 path
    arg 02 - band 02 path
    arg 03 - band 03 path
    arg 04 - band 04 path
    arg 05 - band 05 path
    arg 06 - band 06 path
    arg 07 - band 07 path
    arg 08 - band bqa path
    arg 09 - metadata path
    arg 10 - raster elevation path
    arg 11 - station data path
    arg 12 - output path

*/

/*
int main(int argc, char *argv[]){
    string output_path = argv[12];

    string metadata_path = argv[9];
    MTL mtl = MTL(metadata_path);

    string station_data_path = argv[11];
    Station station = Station(station_data_path, mtl.image_hour);
    
    Sensor sensor = Sensor(mtl.number_sensor, mtl.year);

    TIFF *read_bands[9], *no_shadow_bands[9];
    for(int i = 1; i < 9; i++){
        string path_tiff_base = argv[i];
        read_bands[i] = TIFFOpen(path_tiff_base.c_str(), "rm");

        string path_tiff_write = path_tiff_base.substr(0, path_tiff_base.size() - 4) + "_no_shadow.tif";
        no_shadow_bands[i] = TIFFOpen(path_tiff_write.c_str(), "w8m");

        check_open_tiff(read_bands[i]);

        //Check write bands (ERROR)

        setup(no_shadow_bands[i], read_bands[i]);
    }

    string raster_elevation_path = argv[10];
    TIFF* raster_elevation = TIFFOpen(raster_elevation_path.c_str(), "rm");
    check_open_tiff(raster_elevation);

    string tal_path = tal_function(raster_elevation, output_path);
    TIFFClose(raster_elevation);

    if(!analisy_shadow(read_bands, no_shadow_bands, mtl.number_sensor)){
        close_tifs(read_bands, 9);
        close_tifs(no_shadow_bands, 9);
        cerr << "Shadow problem" << endl;
        exit(5);
    }

    close_tifs(read_bands, 9);
    close_tifs(no_shadow_bands, 9);

    for(int i = 1; i < 8; i++){
        string path_tiff_base = argv[i];
        path_tiff_base = path_tiff_base.substr(0, path_tiff_base.size() - 4) + "_no_shadow.tif";

        read_bands[i] = TIFFOpen(path_tiff_base.c_str(), "rm");
        check_open_tiff(read_bands[i]);
    }

    Landsat landsat = Landsat(tal_path, output_path);
    landsat.process_parcial_products(read_bands, mtl, station, sensor);
    landsat.process_final_products(station, mtl);

    close_tifs(read_bands, 8);
    return 0;
}


/*
    Execution test (after resample function)

    arg 01 - band 01 remsampled path
    arg 02 - band 02 remsampled path
    arg 03 - band 03 remsampled path
    arg 04 - band 04 remsampled path
    arg 05 - band 05 remsampled path
    arg 06 - band 06 remsampled path
    arg 07 - band 07 remsampled path
    arg 08 - metadata path
    arg 09 - tal path
    arg 10 - station data path
    arg 11 - output path

    Flag
    -dist=xxxx - value distance between sun and earth

    ./run input/B2_converted.tif input/B3_converted.tif input/B4_converted.tif input/B5_converted.tif input/B6_converted.tif input/B7_converted.tif input/B10_converted.tif input/MTL.txt tal_converted.tif input/station.csv results -dist=0.98330

*/
int main(int argc, char *argv[]){
    string output_path = argv[11];

    string metadata_path = argv[8];
    MTL mtl = MTL(metadata_path);
    
    string station_data_path = argv[10];
    Station station = Station(station_data_path, mtl.image_hour);
    
    Sensor sensor = Sensor(mtl.number_sensor, mtl.year);

    if(argc == 13){
        string dist_flag = argv[12];
        if(dist_flag.substr(0, 6) == "-dist=")
            mtl.distance_earth_sun = atof(dist_flag.substr(6, dist_flag.size()).c_str());
    }
    
    string tal_path = argv[9];

    TIFF *bands_resampled[8];
    for(int i = 1; i < 8; i++){
        string path_tiff_base = argv[i];
        bands_resampled[i] = TIFFOpen(path_tiff_base.c_str(), "rm");
        check_open_tiff(bands_resampled[i]);
    }

    Landsat landsat = Landsat(tal_path, output_path);
    landsat.process_parcial_products(bands_resampled, mtl, station, sensor);
    landsat.process_final_products(station, mtl);

    close_tifs(bands_resampled, 8);
    return 0;
}

/*
./run input/MTL.txt input/station.csv -dist=0.98330

void print_tiff(TIFF* tif) {

    uint32 heigth_band, width_band;
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &heigth_band);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width_band);

    double tif_line[width_band];

    int cont = 0;

    for(int line = 0; line < heigth_band; line++){
        read_line_tiff(tif, tif_line, line);
        for(int col = 0; col < width_band; col ++){
           printf("%.4lf", tif_line[col]);
           cont++;

           cout << (cont%7 ? " " : "\n");
        }

    }

    cout << endl;

}

void save_tiffs(vector<double*> products_line, vector<TIFF*> products, int line){

    for (unsigned i = 0; i < products.size(); i++)
        write_line_tiff(products[i], products_line[i], line);

};


double getRandomDouble(double min, double max){

    double f = (double) rand() / RAND_MAX;
    return min + f * (max - min);

}

void fill_tiff(TIFF* tif, double min, double max){
    
    uint32 heigth_band, width_band;
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &heigth_band);
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width_band);

    double tif_line[width_band];

    for(int line = 0; line < heigth_band; line++){

        for(int col = 0; col < width_band; col ++){
           tif_line[col] = getRandomDouble(min, max);
        }
        
        if (TIFFWriteScanline(tif, tif_line, line) < 0){
            cerr << "Write problem!" << endl;
            exit(4);
        }
    }
}

void setup(TIFF *tif) {

    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH     , 10);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH    , 10);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE  , 64);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT   , 3);
    TIFFSetField(tif, TIFFTAG_COMPRESSION    , 1);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC    , 1);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP   , 1);
    TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT , 1);
    TIFFSetField(tif, TIFFTAG_XRESOLUTION    , 1);
    TIFFSetField(tif, TIFFTAG_YRESOLUTION    , 1);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG   , PLANARCONFIG_CONTIG);

}

int main(int argc, char *argv[]){
    
    srand (time(NULL));

    string metadata_path = argv[1];
    MTL mtl = MTL(metadata_path);
    
    string station_data_path = argv[2];
    Station station = Station(station_data_path, mtl.image_hour);
    
    Sensor sensor = Sensor(mtl.number_sensor, mtl.year);

    if(argc == 4){
        string dist_flag = argv[3];
        if(dist_flag.substr(0, 6) == "-dist=")
            mtl.distance_earth_sun = atof(dist_flag.substr(6, dist_flag.size()).c_str());
    }

    //TIFF Fake data
    //TIFF *ndvi = TIFFOpen("meuNDVI.tif", "w8m");
    //setup(ndvi);
    //fill_tiff(ndvi, 0.13, 0.40);
    //TIFFClose(ndvi);
    TIFF * ndvi = TIFFOpen("Testes/Teste06/meuNDVI.tif", "rm");

    //TIFF *surface_temperature = TIFFOpen("meuTS.tif", "w8m");
    //setup(surface_temperature);
    //fill_tiff(surface_temperature, 290, 315);
    //TIFFClose(surface_temperature);
    TIFF *surface_temperature = TIFFOpen("Testes/Teste06/meuTS.tif", "rm");

    //TIFF *net_radiation = TIFFOpen("meuRn.tif", "w8m");
    //setup(net_radiation);
    //fill_tiff(net_radiation, 400, 655);
    //TIFFClose(net_radiation);
    TIFF *net_radiation = TIFFOpen("Testes/Teste06/meuRn.tif", "rm");

    //TIFF *soil_heat = TIFFOpen("meuG.tif", "w8m");
    //setup(soil_heat);
    //fill_tiff(soil_heat, 70, 110);
    //TIFFClose(soil_heat);
    TIFF *soil_heat = TIFFOpen("Testes/Teste06/meuG.tif", "rm");

    //TIFF *albedo = TIFFOpen("meuAlbedo.tif", "w8m");
    //setup(albedo);
    //fill_tiff(albedo, 0.15, 0.45);
    //TIFFClose(albedo);
    TIFF *albedo = TIFFOpen("Testes/Teste06/meuAlbedo.tif", "rm");

    Candidate hot_pixel = Candidate();
    hot_pixel.col = 6;
    hot_pixel.line = 4;
    hot_pixel.ndvi = read_position_tiff(ndvi, hot_pixel.col, hot_pixel.line);
    hot_pixel.soil_heat_flux = read_position_tiff(soil_heat, hot_pixel.col, hot_pixel.line);
    hot_pixel.temperature = read_position_tiff(surface_temperature, hot_pixel.col, hot_pixel.line);
    hot_pixel.net_radiation = read_position_tiff(net_radiation, hot_pixel.col, hot_pixel.line);
    
    Candidate cold_pixel = Candidate();
    cold_pixel.col = 8;
    cold_pixel.line = 2;
    cold_pixel.ndvi = read_position_tiff(ndvi, cold_pixel.col, cold_pixel.line);
    cold_pixel.soil_heat_flux = read_position_tiff(soil_heat, cold_pixel.col, cold_pixel.line);
    cold_pixel.temperature = read_position_tiff(surface_temperature, cold_pixel.col, cold_pixel.line);
    cold_pixel.net_radiation = read_position_tiff(net_radiation, cold_pixel.col, cold_pixel.line);

    cout << "HOT" << endl;
    hot_pixel.toString();

    cout << "COLD" << endl;
    cold_pixel.toString();

    uint32 heigth_band, width_band;
    TIFFGetField(albedo, TIFFTAG_IMAGELENGTH, &heigth_band);
    TIFFGetField(albedo, TIFFTAG_IMAGEWIDTH, &width_band);

    double sensible_heat_flux_line[width_band];
    double zom_line[width_band];
    double ustar_line[width_band];
    double aerodynamic_resistence_line[width_band];
    double latent_heat_flux_line[width_band];

    double ustar_station = (VON_KARMAN * station.v6)/(log(station.WIND_SPEED/station.SURFACE_ROUGHNESS));
    double u200 = (ustar_station/VON_KARMAN) * log(200 / station.SURFACE_ROUGHNESS);

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
    zom = TIFFOpen("meuZom.tif", "w8m");
    setup(zom);

    ustar = TIFFOpen("meuUstar.tif", "w8m");
    setup(ustar);

    aerodynamic_resistence = TIFFOpen("meuRah.tif", "w8m");
    setup(aerodynamic_resistence);

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
    
    TIFFClose(ndvi);
    TIFFClose(zom);
    TIFFClose(ustar);
    TIFFClose(aerodynamic_resistence);

    aerodynamic_resistence = TIFFOpen("meuRah.tif", "rm");

    //Extract the hot pixel aerodynamic_resistance
    hot_pixel.aerodynamic_resistance.push_back(read_position_tiff(aerodynamic_resistence, hot_pixel.col, hot_pixel.line));
    double H_hot = hot_pixel.net_radiation - hot_pixel.soil_heat_flux;

    TIFFClose(aerodynamic_resistence);

    TIFF *ustar_tif0, *ustar_tif1, *aerodynamic_resistence_tif0, *aerodynamic_resistence_tif1, *sensible_heat_flux;
    zom = TIFFOpen("meuZom.tif", "rm"); //It's not modified into the rah cycle

    cout << "Zom" << endl;
    print_tiff(zom);

    //It's only written into the rah cycle
    sensible_heat_flux = TIFFOpen("meuH.tif", "w8m");
    setup(sensible_heat_flux);

    int i = 0;
    bool Erro = true;
    
    //Auxiliar arrays calculation      
    double L[width_band];
    double y_01_line[width_band], y_2_line[width_band], x_200_line[width_band];
    double psi_01_line[width_band], psi_2_line[width_band], psi_200_line[width_band];
    double ustar_read_line[width_band], ustar_write_line[width_band], aerodynamic_resistence_read_line[width_band], aerodynamic_resistence_write_line[width_band];

    double rah_hot0;
    double rah_hot;

    //Auxiliar TIFFS
    TIFF *Ltif, *y01, *y2, *x200, *psi01, *psi2, *psi200;

    while(Erro) {
        cout << "Loop " << i << endl;
        rah_hot0 = hot_pixel.aerodynamic_resistance[i];

        if(i%2) {
            //Since ustar is both write and read into the rah cycle, two TIFF will be needed
            ustar_tif0 = TIFFOpen("meuUstar_2.tif", "rm");
            ustar_tif1 = TIFFOpen("meuUstar.tif", "w8m");
            setup(ustar_tif1);
            
            //Since ustar is both write and read into the rah cycle, two TIFF will be needed
            aerodynamic_resistence_tif0 = TIFFOpen("meuRah_2.tif", "rm");
            aerodynamic_resistence_tif1 = TIFFOpen("meuRah.tif", "w8m");
            setup(aerodynamic_resistence_tif1);
            
        } else {
            //Since ustar is both write and read into the rah cycle, two TIFF will be needed
            ustar_tif0 = TIFFOpen("meuUstar.tif", "rm");
            ustar_tif1 = TIFFOpen("meuUstar_2.tif", "w8m");
            setup(ustar_tif1);

            //Since ustar is both write and read into the rah cycle, two TIFF will be needed
            aerodynamic_resistence_tif0 = TIFFOpen("meuRah.tif", "rm");
            aerodynamic_resistence_tif1 = TIFFOpen("meuRah_2.tif", "w8m");
            setup(aerodynamic_resistence_tif1);

        }

        Ltif = TIFFOpen("meuL.tif", "w8m");
        setup(Ltif);

        y01 = TIFFOpen("meuy01.tif", "w8m");
        setup(y01);

        y2 = TIFFOpen("meuy2.tif", "w8m");
        setup(y2);

        x200 = TIFFOpen("meux200.tif", "w8m");
        setup(x200);

        psi01 = TIFFOpen("meupsi01.tif", "w8m");
        setup(psi01);

        psi2 = TIFFOpen("meupsi2.tif", "w8m");
        setup(psi2);

        psi200 = TIFFOpen("meupsi200.tif", "w8m");
        setup(psi200);

        cout << "Ustar before loop" << endl;
        print_tiff(ustar_tif0);
        cout << endl;

        cout << "Rah before loop" << endl;
        print_tiff(aerodynamic_resistence_tif0);
        cout << endl;
        
        for(int line = 0; line < heigth_band; line++){
            
            //Reading data needed
            read_line_tiff(surface_temperature, surface_temperature_line, line);
            read_line_tiff(zom, zom_line, line);
            read_line_tiff(ustar_tif0, ustar_read_line, line);
            read_line_tiff(aerodynamic_resistence_tif0, aerodynamic_resistence_read_line, line);

            double dt_hot = H_hot * rah_hot0 / (RHO * SPECIFIC_HEAT_AIR);
            double b = dt_hot/(hot_pixel.temperature - cold_pixel.temperature);
            double a = -b * (cold_pixel.temperature - 273.15);

            for(int col = 0; col < width_band; col++) {
                sensible_heat_flux_line[col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_line[col] - 273.15))/aerodynamic_resistence_read_line[col];
                
                double ustar_pow_3 = ustar_read_line[col] * ustar_read_line[col] * ustar_read_line[col];
                L[col] = -1 * ((RHO * SPECIFIC_HEAT_AIR * ustar_pow_3 * surface_temperature_line[col])/(VON_KARMAN * GRAVITY * sensible_heat_flux_line[col]));

                y_01_line[col] = pow((1 - (16*0.1)/L[col]), 0.25);
                y_2_line[col] = pow((1 - (16*2)/L[col]), 0.25);
                x_200_line[col] = pow((1 - (16*200)/L[col]), 0.25);

                if(!isnan(L[col]) && L[col] > 0) psi_01_line[col] = -5 * (0.1/L[col]);
                else psi_01_line[col] = 2 * log((1 + y_01_line[col]*y_01_line[col])/2);

                if(!isnan(L[col]) && L[col] > 0) psi_2_line[col] = -5 * (2/L[col]);
                else psi_2_line[col] = 2 * log((1 + y_2_line[col]*y_2_line[col])/2);

                if(!isnan(L[col]) && L[col] > 0) psi_200_line[col] = -5 * (2/L[col]);
                else psi_200_line[col] = 2 * log((1 + x_200_line[col])/2) + log((1 + x_200_line[col]*x_200_line[col])/2) - 2 * atan(x_200_line[col]) + 0.5 * PI;

                ustar_write_line[col] = (VON_KARMAN * u200) / (log(200/zom_line[col]) - psi_200_line[col]);
                aerodynamic_resistence_write_line[col] = (log(2/0.1) - psi_2_line[col] + psi_01_line[col])/(ustar_write_line[col] * VON_KARMAN);

                if(line == hot_pixel.line && col == hot_pixel.col) {
                    rah_hot = aerodynamic_resistence_write_line[col];
                    hot_pixel.aerodynamic_resistance.push_back(rah_hot);
                }
            
            }

            //Saving new ustar e rah
            save_tiffs(vector<double*> {ustar_write_line, aerodynamic_resistence_write_line, L, y_01_line, y_2_line, x_200_line, psi_01_line, psi_2_line, psi_200_line}, 
                    vector<TIFF*> {ustar_tif1, aerodynamic_resistence_tif1, Ltif, y01, y2, x200, psi01, psi2, psi200}, line);

        }

        TIFFClose(Ltif);
        TIFFClose(y01);
        TIFFClose(y2);
        TIFFClose(x200);
        TIFFClose(psi01);
        TIFFClose(psi2);
        TIFFClose(psi200);

        Ltif = TIFFOpen("meuL.tif", "rm");
        y01 = TIFFOpen("meuy01.tif", "rm");
        y2 = TIFFOpen("meuy2.tif", "rm");
        x200 = TIFFOpen("meux200.tif", "rm");
        psi01 = TIFFOpen("meupsi01.tif", "rm");
        psi2 = TIFFOpen("meupsi2.tif", "rm");
        psi200 = TIFFOpen("meupsi200.tif", "rm");

        cout << "L" << endl;
        print_tiff(Ltif);
        cout << endl;

        cout << "y_0.1" << endl;
        print_tiff(y01);
        cout << endl;

        cout << "y_2" << endl;
        print_tiff(y2);
        cout << endl;

        cout << "x_200" << endl;
        print_tiff(x200);
        cout << endl;

        cout << "psi_0.1" << endl;
        print_tiff(psi01);
        cout << endl;

        cout << "psi_2" << endl;
        print_tiff(psi2);
        cout << endl;

        cout << "psi_200" << endl;
        print_tiff(psi200);
        cout << endl;

        TIFFClose(Ltif);
        TIFFClose(y01);
        TIFFClose(y2);
        TIFFClose(x200);
        TIFFClose(psi01);
        TIFFClose(psi2);
        TIFFClose(psi200);

        TIFFClose(ustar_tif0);
        TIFFClose(ustar_tif1);
        TIFFClose(aerodynamic_resistence_tif0);
        TIFFClose(aerodynamic_resistence_tif1);
        
        cout << fabs(1 - rah_hot0/rah_hot) << endl;
        Erro = (fabs(1 - rah_hot0/rah_hot) >= 0.05);
        i++;

    }
    
    TIFFClose(zom);
    
    if(i%2) {
        printf("Rah_after is meuRah_2\n");
        aerodynamic_resistence_tif0 = TIFFOpen("meuRah_2.tif", "rm");
    } else {
        printf("Rah_after is meuRah\n");
        aerodynamic_resistence_tif0 = TIFFOpen("meuRah.tif", "rm");
    }
    cout << "Depois de ler o tif final..." << endl;
    double dt_hot = H_hot * rah_hot0 / (RHO * SPECIFIC_HEAT_AIR);
    double b = dt_hot/(hot_pixel.temperature - cold_pixel.temperature);
    double a = -b * (cold_pixel.temperature - 273.15);

    cout << "H Hot" << H_hot << endl;
    cout << "Rah hot" << rah_hot0 << endl;
    cout << "RHO " << RHO << endl;
    cout << "CP " << SPECIFIC_HEAT_AIR;
    cout << "DT HOT " << dt_hot << endl;
    cout << "B " << b << endl;
    cout << "A " << a << endl;

    for(int line = 0; line < heigth_band; line++){
        read_line_tiff(aerodynamic_resistence_tif0, aerodynamic_resistence_line, line);
        read_line_tiff(surface_temperature, surface_temperature_line, line);
        read_line_tiff(net_radiation, net_radiation_line, line);
        read_line_tiff(soil_heat, soil_heat_line, line);
            
        for(int col = 0; col < width_band; col++) {
            sensible_heat_flux_line[col] = RHO * SPECIFIC_HEAT_AIR * (a + b * (surface_temperature_line[col] - 273.15))/aerodynamic_resistence_line[col];

            if (!isnan(sensible_heat_flux_line[col]) && sensible_heat_flux_line[col] > (net_radiation_line[col] - soil_heat_line[col])) {
                sensible_heat_flux_line[col] = net_radiation_line[col] - soil_heat_line[col];
            }
            

        }

        save_tiffs(vector<double*> {sensible_heat_flux_line}, 
                    vector<TIFF*> {sensible_heat_flux}, line);
    }

    TIFFClose(surface_temperature);
    TIFFClose(aerodynamic_resistence_tif0);
    TIFFClose(sensible_heat_flux);
    
    sensible_heat_flux = TIFFOpen("meuH.tif", "rm");
    print_tiff(sensible_heat_flux);

    TIFF *latent_heat_flux, *net_radiation_24h, *latent_heat_flux_24h, *evapotranspiration_fraction, *evapotranspiration_24h;

    latent_heat_flux = TIFFOpen("meuLatentHF.tif", "w8m");
    setup(latent_heat_flux, albedo);

    net_radiation_24h = TIFFOpen("meuRn24h.tif", "w8m");
    setup(net_radiation_24h, albedo);

    latent_heat_flux_24h = TIFFOpen("meuLatentHF24h.tif", "w8m");
    setup(latent_heat_flux_24h, albedo);

    evapotranspiration_fraction = TIFFOpen("meuEF.tif", "w8m");
    setup(evapotranspiration_fraction, albedo);

    evapotranspiration_24h = TIFFOpen("meuET24h.tif", "w8m");
    setup(evapotranspiration_24h, albedo);

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

    return 0;
}
*/