#include "esaSel.h"

void testLandCoverHomogeneity(TIFF* landCover, TIFF* mask){

    uint32 height_band, width_band;
    TIFFGetField(landCover, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(landCover, TIFFTAG_IMAGEWIDTH, &width_band);

    for(int line = 0; line < height_band; line++) {

        // Create the respective line of the binary map of eligibles pixels
        double mask_line[width_band];

        for(int column = 0; column < width_band; column++) {

            double pixel_value;
            pixel_value = read_position_tiff(landCover, column, line);

            mask_line[column] = false;

            if(pixel_value == AGR) { //Verify if the pixel is an AGR pixel

                mask_line[column] = true;

                for(int i = -3; i <= 3 && mask_line[column]; i++){

                    for(int j = -3; j <= 3 && mask_line[column]; j++){

                        // Check if the neighbor is AGR too

                        if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band) {

                            pixel_value = read_position_tiff(landCover, column + i, line + j);
                            if(!isnan(pixel_value))
                                if(pixel_value != AGR)
                                    mask_line[column] = false;

                        }

                    }

                }

            }

        }

        write_line_tiff(mask, mask_line, line);

    }

}

void testHomogeneity(TIFF* ndvi, TIFF* surface_temperature, TIFF* albedo, TIFF* maskLC, TIFF* output){

    uint32 height_band, width_band;
    TIFFGetField(ndvi, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(ndvi, TIFFTAG_IMAGEWIDTH, &width_band);

    for(int line = 0; line < height_band; line++) {

        // Create the respective line of the binary map of eligibles pixels
        double mask_line[width_band];
        read_line_tiff(maskLC, mask_line, line);

        for(int column = 0; column < width_band; column++) {

            if(mask_line[column] == true) { //Verify if the pixel passed the land cover test

                vector<double> ndvi_neighbors;
                vector<double> ts_neighbors;
                vector<double> albedo_neighbors;

                double pixel_value;

                for(int i = -3; i <= 3; i++){

                    for(int j = -3; j <= 3; j++){

                        // Add for the NDVI, TS and Albedo the value of neighbors pixels into the respective vector

                        if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band) {

                            pixel_value = read_position_tiff(ndvi, column + i, line + j);
                            if(!isnan(pixel_value))
                                ndvi_neighbors.push_back(pixel_value);
                            
                            pixel_value = read_position_tiff(surface_temperature, column + i, line + j);
                            if(!isnan(pixel_value))
                                ts_neighbors.push_back(pixel_value);

                            pixel_value = read_position_tiff(albedo, column + i, line + j);
                            if(!isnan(pixel_value))
                                albedo_neighbors.push_back(pixel_value);

                        }

                    }

                }

                // Do the calculation of the dispersion measures from the NDVI, TS and Albedo

                double mean, sd, coefficient_variation;
                double meanNDVI, meanTS, meanAlb;
                double sdNDVI, sdTS, sdAlb;
                double cvNDVI, cvAlb;
                double sumNDVI = 0, sumTS = 0, sumAlb = 0;

                for(int i = 0; i < ndvi_neighbors.size(); i++) {

                    sumNDVI += ndvi_neighbors[i];
                    sumTS += ts_neighbors[i];
                    sumAlb += albedo_neighbors[i];

                }
                
                meanNDVI = sumNDVI / ndvi_neighbors.size();
                meanTS = sumTS / ts_neighbors.size();
                meanAlb = sumAlb / albedo_neighbors.size();

                sumNDVI = 0, sumTS = 0, sumAlb = 0;

                for(int i = 0; i < ndvi_neighbors.size(); i++) {
                
                    sumNDVI += (ndvi_neighbors[i] - mean) * (ndvi_neighbors[i] - mean);
                    sumTS += (ts_neighbors[i] - mean) * (ts_neighbors[i] - mean);
                    sumAlb += (albedo_neighbors[i] - mean) * (albedo_neighbors[i] - mean);

                }

                sdNDVI = sqrt(sumNDVI / ndvi_neighbors.size());
                sdTS = sqrt(sumTS / ts_neighbors.size());
                sdAlb = sqrt(sumAlb / albedo_neighbors.size());

                cvNDVI = sdNDVI / meanNDVI;
                cvAlb = sdAlb / meanAlb;


                // Check if the pixel is eligible
                mask_line[column] = (cvNDVI < 0.25) && (cvAlb < 0.25) && (sdTS < 1.5);

            }

        }

        write_line_tiff(output, mask_line, line);

    }

}

void testMorphological(TIFF* input, TIFF* output, int groupSize){

    // Read the entire TIFF to the memory
    // Create an same size matrix to serve as output

    uint32 height_band, width_band;
    TIFFGetField(input, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(input, TIFFTAG_IMAGEWIDTH, &width_band);

    int** inputM = (int **) malloc(height_band * sizeof(int *));
    int** outputM = (int **) malloc(height_band * sizeof(int *));

    for(int i = 0; i < height_band; i++){

        inputM[i] = (int *) malloc(width_band * sizeof(int));
        outputM[i] = (int *) malloc(width_band * sizeof(int));

    }

    for(int i = 0; i < height_band; i++){

        read_line_tiff(input, inputM[i], i);
        read_line_tiff(input, outputM[i], i);

    }

    // Apply the routine

    queue< pair<int, int> > fila;
    set< pair<int, int> > cont;

    for(int line = 0; line < height_band; line++) {

        for(int col = 0; col < width_band; col++) {

            if(inputM[line][col] == 1) {

                fila.push({line, col});
                cont.insert({line, col});
                inputM[line][col] = -1;

                while (!fila.empty()) {
                    
                    int i = fila.front().first;
                    int j = fila.front().second;
                    fila.pop();

                    if(j + 1 < width_band){
                        
                        if(inputM[i][j+1] == 1) {
                            fila.push({i, j+1});
                            cont.insert({i, j+1});
                            inputM[i][j+1] = -1;
                        }

                        if(i + 1 < height_band && inputM[i+1][j+1] == 1){
                            fila.push({i+1, j+1});
                            cont.insert({i+1, j+1});
                            inputM[i+1][j+1] = -1;
                        }

                        if(i > 0 && inputM[i-1][j+1] == 1){
                            fila.push({i-1, j+1});
                            cont.insert({i-1, j+1});
                            inputM[i-1][j+1] = -1;
                        }

                    }

                    if(j > 0){

                        if(inputM[i][j-1] == 1){
                            fila.push({i, j-1});
                            cont.insert({i, j-1});
                            inputM[i][j-1] = -1;
                        }

                        if(i + 1 < height_band && inputM[i+1][j-1] == 1){
                            fila.push({i+1, j-1});
                            cont.insert({i+1, j-1});
                            inputM[i+1][j-1] = -1;
                        }

                        if(i > 0 && inputM[i-1][j-1] == 1){
                            fila.push({i-1, j-1});
                            cont.insert({i-1, j-1});
                            inputM[i-1][j-1] = -1;
                        }

                    }

                    if(i + 1 < height_band && inputM[i+1][j] == 1){
                        fila.push({i+1, j});
                        cont.insert({i+1, j});
                        inputM[i+1][j] = -1;
                    }

                    if(i > 0 && inputM[i-1][j] == 1){
                        fila.push({i-1, j});
                        cont.insert({i-1, j});
                        inputM[i-1][j] = -1;
                    }

                }
                
                int group = cont.size();

                for(auto elem : cont) {

                    outputM[elem.first][elem.second] = (group >= groupSize);

                }

                cont.clear();

            } else if (inputM[line][col] == 0) {

                outputM[line][col] = 0;

            }

        }

    }

    // Write output TIFF

    for(int i = 0; i < height_band; i++){

        write_line_tiff(output, outputM[i], i);

    }

    for(int i = 0; i < height_band; i++){

        free(inputM[i]);
        free(outputM[i]);

    }

    free(inputM);
    free(outputM);

}

void hoCalc(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]){

    for(int col = 0; col < width_band; col++)
        ho_line[col] = net_radiation_line[col] - soil_heat_flux[col];

};

pair<Candidate, Candidate> esaPixelSelect(TIFF** ndvi, TIFF** surface_temperature, TIFF** albedo, TIFF** net_radiation, TIFF** soil_heat, TIFF** landCover, int height_band, int width_band, string output_path){

    //Testing land cover homogeneity
    TIFF* outputLC = TIFFOpen((output_path + "/outLC.tif").c_str(), "w8m");
    setup(outputLC, width_band, height_band, 32, 2);
    testLandCoverHomogeneity(*landCover, outputLC);
    TIFFClose(outputLC);

    //Testing the other params homogeneity
    outputLC = TIFFOpen((output_path + "/outLC.tif").c_str(), "r");
    TIFF* outputH = TIFFOpen((output_path + "/outH.tif").c_str(), "w8m");
    setup(outputH, width_band, height_band, 32, 2);
    testHomogeneity(*ndvi, *surface_temperature, *albedo, outputLC, outputH);
    TIFFClose(outputLC);
    TIFFClose(outputH);

    //Morphological test
    outputH = TIFFOpen((output_path + "/outH.tif").c_str(), "r");
    TIFF* outputAll = TIFFOpen((output_path + "/outAll.tif").c_str(), "w8m");
    setup(outputAll, width_band, height_band, 32, 2);
    testMorphological(outputH, outputAll, 50);
    TIFFClose(outputH);
    TIFFClose(outputAll);

    //Getting the candidates
    outputAll = TIFFOpen((output_path + "/outAll.tif").c_str(), "r");
    int all_condition[width_band];
    vector<Candidate> histTS;

    //Auxiliary arrays
    double ndvi_line[width_band], surface_temperature_line[width_band];
    double net_radiation_line[width_band], soil_heat_line[width_band];
    double ho_line[width_band];

    //Creating candidates array for TS and then for NDVI as a copy
    for(int line = 0; line < height_band; line++){

        read_line_tiff(outputAll, all_condition, line);

        read_line_tiff(*net_radiation, net_radiation_line, line);
        read_line_tiff(*soil_heat, soil_heat_line, line);

        hoCalc(net_radiation_line, soil_heat_line, width_band, ho_line);

        read_line_tiff(*ndvi, ndvi_line, line);
        read_line_tiff(*surface_temperature, surface_temperature_line, line);

        for(int col = 0; col < width_band; col++) {

            if(all_condition[col]){

                histTS.push_back(Candidate(ndvi_line[col],
                                           surface_temperature_line[col],
                                           net_radiation_line[col],
                                           soil_heat_line[col],
                                           ho_line[col],
                                           line, col));

            }

        }

    }

    vector<Candidate> histNDVI (histTS);

    sort(histTS.begin(), histTS.end(), compare_candidate_temperature);
    sort(histNDVI.begin(), histNDVI.end(), compare_candidate_ndvi);

    cout << "TS SIZE: " << histTS.size() << ", NDVI SIZE: " << histNDVI.size() << endl;


    // Select cold pixel
    int pixel_count = 0, n1 = 1, n2 = 1, ts_pos, ndvi_pos;
    vector<Candidate> coldPixels;
    while (pixel_count < 10 && !(n2 == 10 && n1 == 10)) {
        
        ts_pos = int(floor(n1/100.0 * histTS.size()));
        ndvi_pos = int(floor((100 - n2)/100.0 * histNDVI.size()));

        cout << "TS POS: " << ts_pos << ", NDVI POS: " << ndvi_pos << endl;

        for(int i = 0; i <= ts_pos; i++) {

            for(int j = histNDVI.size() - 1; j >= ndvi_pos; j--) {

                if(equals(histTS[i], histNDVI[j])) coldPixels.push_back(histTS[i]);

            }

        }

        if(n2 < 10) n2++;
        else if(n1 < 10) n1++;

    }

    //Select hot pixel
    pixel_count = 0, n1 = 1, n2 = 1;
    vector<Candidate> hotPixels;
    while (pixel_count < 10 && !(n2 == 10 && n1 == 10)) {
        
        ts_pos = int(floor((100 - n1)/100.0 * histTS.size()));
        ndvi_pos = int(floor(n2/100.0 * histNDVI.size()));

        cout << "TS POS: " << ts_pos << ", NDVI POS: " << ndvi_pos << endl;

        for(int i = 0; i <= ndvi_pos; i++) {

            for(int j = histTS.size() - 1; j >= ts_pos; j--) {

                if(equals(histTS[j], histNDVI[i])) coldPixels.push_back(histTS[j]);

            }

        }

        if(n2 < 10) n2++;
        else if(n1 < 10) n1++;

    }

    cout << "COLD SIZE: " << coldPixels.size() << ", HOT SIZE: " << hotPixels.size() << endl;

    sort(coldPixels.begin(), coldPixels.end(), compare_candidate_ndvi);
    sort(hotPixels.begin(), hotPixels.end(), compare_candidate_temperature);

    return {hotPixels[hotPixels.size() - 1], coldPixels[0]};
   
}
