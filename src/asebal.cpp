#include "asebal.h"

void quartile(TIFF** target, double* vQuartile){

    uint32 height_band, width_band;
    TIFFGetField(target, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(target, TIFFTAG_IMAGEWIDTH, &width_band);

    const int SIZE = height_band * width_band;
    double target_line[width_band];
    double* target_values = (double*) malloc(sizeof(double) * SIZE);
    if(target_values == NULL) exit(15);
    int pos = 0;

    for(int line = 0; line < height_band; line++){
        read_line_tiff(target, target_line, line);

        for(int col = 0; col < width_band; col++){

            target_values[pos] = target_line[col];
            pos++;

        }

    }

    sort(target_values, target_values + SIZE);

    //First quartile
    vQuartile[0] = target_values[int(floor(0.25 * SIZE))];

    //Second quartile
    vQuartile[1] = target_values[int(floor(0.5 * SIZE))];

    //Third quartile
    vQuartile[2] = target_values[int(floor(0.75 * SIZE))];

    //Fourth quartile
    vQuartile[3] = target_values[SIZE-1];

}

Candidate getHotPixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band){

    double ndvi_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band], soil_heat_line[width_band];
    double ho_line[width_band], albedo_line[width_band];

    double* ndviQuartile = (double*) malloc(sizeof(double) * 4);
    double* tsQuartile = (double*) malloc(sizeof(double) * 4);
    double* albedoQuartile = (double*) malloc(sizeof(double) * 4);

    quartile(&ndvi, ndviQuartile);
    quartile(&surface_temperature, ndviQuartile);
    quartile(&albedo, ndviQuartile);

    //Creating first pixel group
    vector<Candidate> candidatesGroupI;

    for(int line = 0; line < height_band; line ++){
        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        ho_function(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);
        read_line_tiff(*albedo, albedo, line);

        for(int col = 0; col < width_band; col ++){

            bool albedoValid = !isnan(albedo_line[col]) && albedo_line[col] > albedoQuartile[2];
            bool ndviValid = !isnan(ndvi_line[col]) && ndvi_line[col] > 0.10 && ndvi_line[col] < ndviQuartile[0];
            bool tsValid = !isnan(surface_temperature_line[col]) && surface_temperature_line[col] > tsQuartile[2];

            if(albedoValid && ndviValid && tsValid){
                candidatesGroupI.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col],
                                    ho_line[col],
                                    line, col));
            }
        }

    }

    //Creating second pixel group, all values lower than the 3rd quartile are excluded
    sort(candidatesGroupI.begin(), candidatesGroupI.end(), compare_candidate_temperature);
    unsigned int pos = int(floor(candidatesGroupI.size() * 0.75));
    vector<Candidate> candidatesGroupII(candidatesGroupI.begin() + pos, candidatesGroupI.end());

    pos = int(floor(candidatesGroupII.size() * 0.5));
    Candidate hotPixel = candidatesGroupII[pos];
    return hotPixel;
}

Candidate getColdPixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band) {
    
    double ndvi_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band], soil_heat_line[width_band];
    double ho_line[width_band], albedo_line[width_band];

    double* ndviQuartile = (double*) malloc(sizeof(double) * 4);
    double* tsQuartile = (double*) malloc(sizeof(double) * 4);
    double* albedoQuartile = (double*) malloc(sizeof(double) * 4);

    quartile(&ndvi, ndviQuartile);
    quartile(&surface_temperature, ndviQuartile);
    quartile(&albedo, ndviQuartile);

    //Creating first pixel group
    vector<Candidate> candidatesGroupI;

    for(int line = 0; line < height_band; line ++){
        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        ho_function(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);
        read_line_tiff(*albedo, albedo, line);

        for(int col = 0; col < width_band; col ++){

            bool albedoValid = !isnan(albedo_line[col]) && albedo_line[col] < albedoQuartile[1];
            bool ndviValid = !isnan(ndvi_line[col]) && ndvi_line[col] >= ndviQuartile[3];
            bool tsValid = !isnan(surface_temperature_line[col]) && surface_temperature_line[col] < tsQuartile[0];

            if(albedoValid && ndviValid && tsValid){
                candidatesGroupI.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col],
                                    ho_line[col],
                                    line, col));
            }
        }

    }

    //Creating second pixel group, all pixels with values higher than the 1 st quartile are excluded
    sort(candidatesGroupI.begin(), candidatesGroupI.end(), compare_candidate_temperature);
    unsigned int pos = int(floor(candidatesGroupI.size() * 0.25));
    vector<Candidate> candidatesGroupII(candidatesGroupI.begin(), candidatesGroupI.begin() + pos);

    pos = int(floor(candidatesGroupII.size() * 0.5));
    Candidate coldPixel = candidatesGroupII[pos];
    return coldPixel;
}