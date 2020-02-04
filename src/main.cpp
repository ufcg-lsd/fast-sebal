#include "utils.h"
#include "landsat.h"
#include "products.h"

using namespace std;

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

    string landCoverPath = (argc >= 13) ? argv[12] : "";
    printf("PATH: %s\n", landCoverPath.c_str());

    int method = 0;
    if(argc >= 14){
        string flag = argv[13];
        if(flag.substr(0, 6) == "-meth=")
            method = flag[6] - '0';
    }

    double noData = NaN;
    if(argc >= 15){
        string noData_flag = argv[14];
        if(noData_flag.substr(0,5) == "-nan=")
            noData = atof(noData_flag.substr(5, noData_flag.size()).c_str());
    }
    
    string tal_path = argv[9];

    TIFF *bands_resampled[8];
    for(int i = 1; i < 8; i++){
        string path_tiff_base = argv[i];
        bands_resampled[i] = TIFFOpen(path_tiff_base.c_str(), "rm");
        check_open_tiff(bands_resampled[i]);
    }

    //Timing
    chrono::steady_clock::time_point begin, end;
    chrono::duration< double, micro > time_span_us;

    Landsat landsat = Landsat(tal_path, output_path, method, noData, landCoverPath);
    //printf("PHASE 1 - START, %d\n", int(time(NULL)));
    begin = chrono::steady_clock::now();
    landsat.process_partial_products(bands_resampled, mtl, station, sensor);
    end = chrono::steady_clock::now();
    time_span_us = chrono::duration_cast< chrono::duration<double, micro> >(end - begin);
    printf("PHASE 1 - DURATION, %.5f\n", time_span_us);

    begin = chrono::steady_clock::now();
    //printf("PHASE 2 - START, %d\n", int(time(NULL)));
    landsat.process_final_products(station, mtl);
    close_tiffs(bands_resampled, 8);
    end = chrono::steady_clock::now();
    time_span_us = chrono::duration_cast< chrono::duration<double, micro> >(end - begin);
    printf("PHASE 2 - DURATION, %.5f\n", time_span_us);
    return 0;
}