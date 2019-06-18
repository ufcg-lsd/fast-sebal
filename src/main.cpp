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