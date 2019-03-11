#include "utils.h"
#include "landsat.h"
#include "products.h"

using namespace std;

int main(int argc, char *argv[]){
    string output_path = argv[12];

    string metadata_path = argv[9];
    MTL mtl = MTL(metadata_path);

    string dados_estacao_path = argv[11];
    Estacao estacao = Estacao(dados_estacao_path);
    
    Sensor sensor = Sensor(mtl.number_sensor, mtl.year);

    TIFF *read_bands[9], *no_shadow_bands[9];
    for(int i = 1; i < 9; i++){
        string path_tiff_base = argv[i];
        read_bands[i] = TIFFOpen(path_tiff_base.c_str(), "rm");

        string path_tiff_write = path_tiff_base.substr(0, path_tiff_base.size() - 4) + "_no_shadow.tif";
        no_shadow_bands[i] = TIFFOpen(path_tiff_write.c_str(), "w8m");

        if(!read_bands[i]){
            cerr << "Open problem" << endl;
            return 1;
        } 

        //Check write bands (ERROR)

        setup(no_shadow_bands[i], read_bands[i]);
    }

    string raster_elevation_path = argv[10];
    TIFF* raster_elevation = TIFFOpen(raster_elevation_path.c_str(), "rm");
    if(!raster_elevation){
        cerr << "Open raster elevation problem" << endl;
        return 1;
    }
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
        if(!read_bands[i]){
            cerr << "Open problem" << endl;
            return 1;
        }
    }

    Landsat landsat = Landsat(tal_path, output_path);
    landsat.process(read_bands, mtl, sensor);

    //calculate ET, ET24h

    close_tifs(read_bands, 8);
    return 0;
}