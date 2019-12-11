#include "asebal.h"

void quartile(TIFF* target, double* vQuartile, int height_band, int width_band){

    const int SIZE = height_band * width_band;
    double target_line[width_band];
    double* target_values = (double*) malloc(sizeof(double) * SIZE);
    if(target_values == NULL) exit(15);
    int pos = 0;

    for(int line = 0; line < height_band; line++){
        read_line_tiff(target, target_line, line);

        for(int col = 0; col < width_band; col++){

            if(!isnan(target_line[col]) && !isinf(target_line[col])) {

                target_values[pos] = target_line[col];
                pos++;

            }

        }

    }

    sort(target_values, target_values + pos);

    //First quartile
    vQuartile[0] = target_values[int(floor(0.25 * pos))];

    //Second quartile
    vQuartile[1] = target_values[int(floor(0.5 * pos))];

    //Third quartile
    vQuartile[2] = target_values[int(floor(0.75 * pos))];

    //Fourth quartile
    vQuartile[3] = target_values[pos-1];

    printf("%.3f %.3f %.3f %.3f\n", vQuartile[0], vQuartile[1], vQuartile[2], vQuartile[3]);

    free(target_values);

}

void hoFunction(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]){

    for(int col = 0; col < width_band; col++)
        ho_line[col] = net_radiation_line[col] - soil_heat_flux[col];

};

Candidate getHotPixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band){

    double ndvi_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band], soil_heat_line[width_band];
    double ho_line[width_band], albedo_line[width_band];

    double* ndviQuartile = (double*) malloc(sizeof(double) * 4);
    double* tsQuartile = (double*) malloc(sizeof(double) * 4);
    double* albedoQuartile = (double*) malloc(sizeof(double) * 4);

    quartile(*ndvi, ndviQuartile, height_band, width_band);
    quartile(*surface_temperature, tsQuartile, height_band, width_band);
    quartile(*albedo, albedoQuartile, height_band, width_band);
    //Creating first pixel group
    vector<Candidate> candidatesGroupI;

    for(int line = 0; line < height_band; line ++){
        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        hoFunction(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);
        read_line_tiff(*albedo, albedo_line, line);

        for(int col = 0; col < width_band; col++){

            bool albedoValid = !isnan(albedo_line[col]) && albedo_line[col] > albedoQuartile[2];
            bool ndviValid = !isnan(ndvi_line[col]) && ndvi_line[col] < ndviQuartile[0]; //ndvi_line[col] > 0.10 && 
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

    free(ndviQuartile);
    free(tsQuartile);
    free(albedoQuartile);

    hotPixel.toString();

    return hotPixel;
}

Candidate getColdPixel(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, int height_band, int width_band) {
    
    double ndvi_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band], soil_heat_line[width_band];
    double ho_line[width_band], albedo_line[width_band];

    double* ndviQuartile = (double*) malloc(sizeof(double) * 4);
    double* tsQuartile = (double*) malloc(sizeof(double) * 4);
    double* albedoQuartile = (double*) malloc(sizeof(double) * 4);

    quartile(*ndvi, ndviQuartile, height_band, width_band);
    quartile(*surface_temperature, tsQuartile, height_band, width_band);
    quartile(*albedo, albedoQuartile, height_band, width_band);

    //Creating first pixel group
    vector<Candidate> candidatesGroupI;

    for(int line = 0; line < height_band; line ++){
        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        hoFunction(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);
        read_line_tiff(*albedo, albedo_line, line);

        for(int col = 0; col < width_band; col ++){

            bool albedoValid = !isnan(albedo_line[col]) && albedo_line[col] < albedoQuartile[1];
            bool ndviValid = !isnan(ndvi_line[col]) && ndvi_line[col] >= ndviQuartile[2]; //ndvi_line[col] >= ndviQuartile[3];
            bool tsValid = !isnan(surface_temperature_line[col]) && surface_temperature_line[col] < tsQuartile[0];

            if(albedoValid && ndviValid && tsValid){

                // printf("NDVI: %.10lf\n", ndvi_line[col]);
                // printf("TS: %.10lf\n", surface_temperature_line[col]);
                // printf("Rn: %.10lf\n", net_radiation_line[col]);
                // printf("G: %.10lf\n", soil_heat_line[col]);
                // printf("HO: %.10lf\n", ho_line[col]);
                // printf("Line: %d\n", line);
                // printf("Col: %d\n", col);

                candidatesGroupI.push_back(Candidate(ndvi_line[col],
                                    surface_temperature_line[col],
                                    net_radiation_line[col],
                                    soil_heat_line[col],
                                    ho_line[col],
                                    line, col));

                //candidatesGroupI[candidatesGroupI.size()-1].toString();
            }
        }

    }
    std::cout << "GROUP1: " << candidatesGroupI.size() << std::endl;
    //Creating second pixel group, all pixels with values higher than the 1 st quartile are excluded
    sort(candidatesGroupI.begin(), candidatesGroupI.end(), compare_candidate_temperature);

    // for(unsigned int i = 0; i < candidatesGroupI.size(); i++) {

    //     candidatesGroupI[i].toString();

    // }
    unsigned int pos = int(floor(candidatesGroupI.size() * 0.25));
    vector<Candidate> candidatesGroupII(candidatesGroupI.begin(), candidatesGroupI.begin() + pos);
    std::cout << "GROUP2: " << candidatesGroupII.size() << std::endl;
    // for(unsigned int i = 0; i < candidatesGroupII.size(); i++) {

    //     candidatesGroupII[i].toString();

    // }
    pos = int(floor(candidatesGroupII.size() * 0.5));
    Candidate coldPixel = candidatesGroupII[pos];

    free(ndviQuartile);
    free(tsQuartile);
    free(albedoQuartile);

    coldPixel.toString();

    return coldPixel;
}