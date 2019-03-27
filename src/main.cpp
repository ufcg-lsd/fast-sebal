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
*/

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

    /*  DEBUG
    printf("Image hour: %.10lf\n", mtl.image_hour);
    printf("Julian day: %d\n", mtl.julian_day);
    printf("Number sensor: %d\n", mtl.number_sensor);
    printf("Rad 10 add: %.10lf\n", mtl.rad_add_10);
    printf("Rad 10 mult: %.10lf\n", mtl.rad_mult_10);
    printf("Sun elevation: %.10lf\n", mtl.sun_elevation);
    printf("Costheta: %.10f\n", sin(mtl.sun_elevation * PI / 180));
    printf("Year: %d\n", mtl.year);
    */
    
    string station_data_path = argv[10];
    Station station = Station(station_data_path, mtl.image_hour);

    /*  DEGUG
    printf("V6 station: %.10lf\n", station.v6);
    printf("V7 max station: %.10lf\n", station.v7_max);
    printf("V7 min station: %.10lf\n", station.v7_min);
    printf("Temperature image: %.10lf\n", station.temperature_image);
    printf("Latitude: %.10lf\n", station.latitude);
    printf("Longitude: %.10lf\n", station.longitude);
   

    double ustar_station = (VON_KARMAN * station.v6)/(log(station.WIND_SPEED/station.SURFACE_ROUGHNESS));
    double u200 = (ustar_station/VON_KARMAN) * log(200 / station.SURFACE_ROUGHNESS);

    printf("Ustar station: %.10lf\n", ustar_station);
    printf("u200: %.10lf\n", u200);

    double ndvi_hot_pixel = 0.1930204209;
    double zom_hot_pixel = exp(station.A_ZOM + station.B_ZOM * ndvi_hot_pixel);
    double ustar_hot_pixel = (VON_KARMAN * u200)/log(200/zom_hot_pixel);

    printf("Zom hot pixel: %.10lf\n", zom_hot_pixel);
    printf("Ustar hot pixel: %.10lf\n", ustar_hot_pixel);

    double ndvi_cold_pixel = -0.0010276547;
    double zom_cold_pixel = exp(station.A_ZOM + station.B_ZOM * ndvi_cold_pixel);
    double ustar_cold_pixel = (VON_KARMAN * u200)/log(200/zom_cold_pixel);

    printf("Zom cold pixel: %.10lf\n", zom_cold_pixel);
    printf("Ustar cold pixel: %.10lf\n", ustar_cold_pixel);

    return 0;
     */
    
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
double getRandomDouble(double min, double max){

    double f = (double) rand() / RAND_MAX;
    return min + f * (max - min);

}

void fill_tiff(TIFF** tif, double min, double max){
    
    double tif_line[30];

    for(int line = 0; line < 30; line++){

        for(int col = 0; col < 30; col ++){
           tif_line[col] = getRandomDouble(min, max);
        }

        if (TIFFWriteScanline(*tif, tif_line, line) < 0){
            cerr << "Write problem!" << endl;
            exit(4);
        }
    }
}

int main(){
    srand (time(NULL));

    TIFF *ndvi, *soil_heat, *surface_temperature, *net_radiation;

    ndvi = TIFFOpen("/local/workspace/saps/fast-sebal/testeS/meuNDVI.tif", "w8m");
    TIFFSetField(ndvi, TIFFTAG_IMAGEWIDTH     , 30); 
    TIFFSetField(ndvi, TIFFTAG_IMAGELENGTH    , 30);
    TIFFSetField(ndvi, TIFFTAG_BITSPERSAMPLE  , 64);
    TIFFSetField(ndvi, TIFFTAG_SAMPLEFORMAT   , 3);
    TIFFSetField(ndvi, TIFFTAG_COMPRESSION    , 1);
    TIFFSetField(ndvi, TIFFTAG_PHOTOMETRIC    , 1);
    TIFFSetField(ndvi, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(ndvi, TIFFTAG_ROWSPERSTRIP   , 1);
    TIFFSetField(ndvi, TIFFTAG_RESOLUTIONUNIT , 1);
    TIFFSetField(ndvi, TIFFTAG_XRESOLUTION    , 1);
    TIFFSetField(ndvi, TIFFTAG_YRESOLUTION    , 1);
    TIFFSetField(ndvi, TIFFTAG_PLANARCONFIG   , PLANARCONFIG_CONTIG);

    soil_heat = TIFFOpen("/local/workspace/saps/fast-sebal/testeS/meuG.tif", "w8m");
    TIFFSetField(soil_heat, TIFFTAG_IMAGEWIDTH     , 30); 
    TIFFSetField(soil_heat, TIFFTAG_IMAGELENGTH    , 30);
    TIFFSetField(soil_heat, TIFFTAG_BITSPERSAMPLE  , 64);
    TIFFSetField(soil_heat, TIFFTAG_SAMPLEFORMAT   , 3);
    TIFFSetField(soil_heat, TIFFTAG_COMPRESSION    , 1);
    TIFFSetField(soil_heat, TIFFTAG_PHOTOMETRIC    , 1);
    TIFFSetField(soil_heat, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(soil_heat, TIFFTAG_ROWSPERSTRIP   , 1);
    TIFFSetField(soil_heat, TIFFTAG_RESOLUTIONUNIT , 1);
    TIFFSetField(soil_heat, TIFFTAG_XRESOLUTION    , 1);
    TIFFSetField(soil_heat, TIFFTAG_YRESOLUTION    , 1);
    TIFFSetField(soil_heat, TIFFTAG_PLANARCONFIG   , PLANARCONFIG_CONTIG);

    surface_temperature = TIFFOpen("/local/workspace/saps/fast-sebal/testeS/meuTS.tif", "w8m");
    TIFFSetField(surface_temperature, TIFFTAG_IMAGEWIDTH     , 30); 
    TIFFSetField(surface_temperature, TIFFTAG_IMAGELENGTH    , 30);
    TIFFSetField(surface_temperature, TIFFTAG_BITSPERSAMPLE  , 64);
    TIFFSetField(surface_temperature, TIFFTAG_SAMPLEFORMAT   , 3);
    TIFFSetField(surface_temperature, TIFFTAG_COMPRESSION    , 1);
    TIFFSetField(surface_temperature, TIFFTAG_PHOTOMETRIC    , 1);
    TIFFSetField(surface_temperature, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(surface_temperature, TIFFTAG_ROWSPERSTRIP   , 1);
    TIFFSetField(surface_temperature, TIFFTAG_RESOLUTIONUNIT , 1);
    TIFFSetField(surface_temperature, TIFFTAG_XRESOLUTION    , 1);
    TIFFSetField(surface_temperature, TIFFTAG_YRESOLUTION    , 1);
    TIFFSetField(surface_temperature, TIFFTAG_PLANARCONFIG   , PLANARCONFIG_CONTIG);

    net_radiation = TIFFOpen("/local/workspace/saps/fast-sebal/testeS/meuRn.tif", "w8m");
    TIFFSetField(net_radiation, TIFFTAG_IMAGEWIDTH     , 30); 
    TIFFSetField(net_radiation, TIFFTAG_IMAGELENGTH    , 30);
    TIFFSetField(net_radiation, TIFFTAG_BITSPERSAMPLE  , 64);
    TIFFSetField(net_radiation, TIFFTAG_SAMPLEFORMAT   , 3);
    TIFFSetField(net_radiation, TIFFTAG_COMPRESSION    , 1);
    TIFFSetField(net_radiation, TIFFTAG_PHOTOMETRIC    , 1);
    TIFFSetField(net_radiation, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(net_radiation, TIFFTAG_ROWSPERSTRIP   , 1);
    TIFFSetField(net_radiation, TIFFTAG_RESOLUTIONUNIT , 1);
    TIFFSetField(net_radiation, TIFFTAG_XRESOLUTION    , 1);
    TIFFSetField(net_radiation, TIFFTAG_YRESOLUTION    , 1);
    TIFFSetField(net_radiation, TIFFTAG_PLANARCONFIG   , PLANARCONFIG_CONTIG);

    fill_tiff(&ndvi, -0.2, 0.8);
    fill_tiff(&soil_heat, 20.0, 380.0);
    fill_tiff(&surface_temperature, 280.0, 315.0);
    fill_tiff(&net_radiation, 90.0, 760.0);

    TIFFClose(ndvi);
    TIFFClose(soil_heat);
    TIFFClose(surface_temperature);
    TIFFClose(net_radiation);

    ndvi = TIFFOpen("/local/workspace/saps/fast-sebal/testeS/meuNDVI.tif", "rm");
    soil_heat = TIFFOpen("/local/workspace/saps/fast-sebal/testeS/meuG.tif", "rm");
    surface_temperature = TIFFOpen("/local/workspace/saps/fast-sebal/testeS/meuTS.tif", "rm");
    net_radiation = TIFFOpen("/local/workspace/saps/fast-sebal/testeS/meuRn.tif", "rm");

    Candidate hot = select_hot_pixel(&ndvi, &surface_temperature, &net_radiation, &soil_heat, 30, 30);
    Candidate cold = select_cold_pixel(&ndvi, &surface_temperature, &net_radiation, &soil_heat, 30, 30);
    hot.toString();
    cold.toString();

    return 0;
}*/