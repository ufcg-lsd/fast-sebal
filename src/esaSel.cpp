#include "esaSel.h"

/**
 * @brief   Checks if the input value matches one of the agricultural land cover code.
 * @param   value: The input value.
 * @retval  TRUE if the values matches and FALSE otherwise.
 */
bool checkLandCode(int value){
    
    return (value == AGP) || (value == PAS) || (value == AGR) || (value == CAP) || (value == CSP) || (value == MAP);

}

/**
 * @brief   Tests the land cover homogeneity of a pixel.
 * @note    A the land cover of a pixel is homogeneous if every neighbour pixel inside a 7x7 window is also agricultural field.
 * @param   landCover: Land Cover TIFF.
 * @param   mask: Output binary TIFF, where pixels with 1 means that is a homogeneous pixel and 0 means the otherwise.
 */
void testLandCoverHomogeneity(TIFF* landCover, TIFF* mask){

    uint32 height_band, width_band;
    TIFFGetField(landCover, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(landCover, TIFFTAG_IMAGEWIDTH, &width_band);

    double** buffer = (double **) malloc(7 * sizeof(double *));

    for(int i = 0; i < 7; i++){
        buffer[i] = (double *) malloc(width_band * sizeof(double));
    }

    int relation[7] = {-1, -1, -1, -1, -1, -1, -1}, aux;

    for(int line = 0; line < height_band; line++) {

        // Create the respective line of the binary map of eligibles pixels
        int mask_line[width_band];

        for(int column = 0; column < width_band; column++) {

            int pixel_value;

            aux = line % 7;

            if(relation[aux] != line) {
                
                read_line_tiff(landCover, buffer[aux], line);
                relation[aux] = line;

            }

            pixel_value = buffer[aux][column];

            mask_line[column] = false;

            if(checkLandCode(pixel_value)) { //Verify if the pixel is an AGR pixel

                mask_line[column] = true;

                for(int i = -3; i <= 3 && mask_line[column]; i++){

                    for(int j = -3; j <= 3 && mask_line[column]; j++){

                        // Check if the neighbor is AGR too

                        if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band) {

                            aux = (line + j) % 7;

                            if(relation[aux] != (line + j)) {

                                read_line_tiff(landCover, buffer[aux], line + j);
                                relation[aux] = (line + j);

                            }

                            pixel_value = buffer[aux][column + i];

                            if(!isnan(pixel_value))
                                if(!checkLandCode(pixel_value))
                                    mask_line[column] = false;

                        }

                    }

                }

            }

        }

        write_line_tiff(mask, mask_line, line);

    }

    for(int i = 0; i < 7; i++){
        free(buffer[i]);
    }
    free(buffer);

}

/**
 * @brief   Tests the ndvi, surface_temperature and albedo homogeneity of a pixel.
 * @note    A pixel is homogeneous in these criteria if inside a 7x7 window the coefficient of variation of the albedo and ndvi is less or equal than 25%
 *          and the surface temperature has a standard deviation less or equal than 1.5 K.
 * @param   ndvi: NDVI TIFF.
 * @param   surface_temperature: TS TIFF.
 * @param   albedo: Albedo TIFF.
 * @param   maskLC: A binary TIFF conteining the data of the land cover homogeneity.
 * @param   output: A binary TIFF, where pixels with 1 means that is a homogeneous pixel in land cover, ndvi, surface temperature and albedo, and 0 means otherwise.
 */
void testHomogeneity(TIFF* ndvi, TIFF* surface_temperature, TIFF* albedo, TIFF* maskLC, TIFF* output){

    uint32 height_band, width_band;
    TIFFGetField(ndvi, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(ndvi, TIFFTAG_IMAGEWIDTH, &width_band);

    double **bufferTS = (double **) malloc(10 * sizeof(double *));
    double **bufferNDVI = (double **) malloc(10 * sizeof(double *));
    double **bufferAlb = (double **) malloc(10 * sizeof(double *));

    for(int i = 0; i < 7; i++){
        bufferTS[i] = (double *) malloc(width_band * sizeof(double));
        bufferNDVI[i] = (double *) malloc(width_band * sizeof(double));
        bufferAlb[i] = (double *) malloc(width_band * sizeof(double));
    }

    int relation[7] = {-1, -1, -1, -1, -1, -1, -1}, aux;

    for(int line = 0; line < height_band; line++) {

        // Create the respective line of the binary map of eligibles pixels
        int mask_line[width_band];
        read_line_tiff(maskLC, mask_line, line);

        for(int column = 0; column < width_band; column++) {

            if(mask_line[column] == true) { //Verify if the pixel passed the land cover test

                vector<double> ndvi_neighbors;
                vector<double> ts_neighbors;
                vector<double> albedo_neighbors;

                double pixel_value;

                aux = line % 7;

                if(relation[aux] != line) {

                    read_line_tiff(surface_temperature, bufferTS[aux], line);
                    read_line_tiff(ndvi, bufferNDVI[aux], line);
                    read_line_tiff(albedo, bufferAlb[aux], line);
                    relation[aux] = line;

                }

                if(!isnan(bufferNDVI[aux][column])){

                    for(int i = -3; i <= 3; i++){

                        for(int j = -3; j <= 3; j++){

                            // Add for the NDVI, TS and Albedo the value of neighbors pixels into the respective vector

                            if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band) {

                                aux = (line + j) % 7;

                                if(relation[aux] != (line + j)) {

                                    read_line_tiff(surface_temperature, bufferTS[aux], line + j);
                                    read_line_tiff(ndvi, bufferNDVI[aux], line + j);
                                    read_line_tiff(albedo, bufferAlb[aux], line + j);
                                    relation[aux] = line + j;

                                }

                                pixel_value = bufferNDVI[aux][column+i]; 
                                if(!isnan(pixel_value))
                                    ndvi_neighbors.push_back(pixel_value);
                                
                                pixel_value = bufferTS[aux][column + i]; 
                                if(!isnan(pixel_value))
                                    ts_neighbors.push_back(pixel_value);

                                pixel_value = bufferAlb[aux][column + i];
                                if(!isnan(pixel_value))
                                    albedo_neighbors.push_back(pixel_value);

                            }

                        }

                    }

                    // Do the calculation of the dispersion measures from the NDVI, TS and Albedo

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
                    
                        sumNDVI += (ndvi_neighbors[i] - meanNDVI) * (ndvi_neighbors[i] - meanNDVI);
                        sumTS += (ts_neighbors[i] - meanTS) * (ts_neighbors[i] - meanTS);
                        sumAlb += (albedo_neighbors[i] - meanAlb) * (albedo_neighbors[i] - meanAlb);

                    }

                    sdNDVI = sqrt(sumNDVI / ndvi_neighbors.size());
                    sdTS = sqrt(sumTS / ts_neighbors.size());
                    sdAlb = sqrt(sumAlb / albedo_neighbors.size());

                    cvNDVI = sdNDVI / meanNDVI;
                    cvAlb = sdAlb / meanAlb;


                    // Check if the pixel is eligible
                    mask_line[column] = (cvNDVI < 0.25) && (cvAlb < 0.25) && (sdTS < 1.5);

                } else {

                    mask_line[column] = false;

                }

            }

        }

        write_line_tiff(output, mask_line, line);

    }

    for(int i = 0; i < 7; i++){
        free(bufferTS[i]);
        free(bufferNDVI[i]);
        free(bufferAlb[i]);
    }
    free(bufferAlb);
    free(bufferNDVI);
    free(bufferTS);

}

/**
 * @brief   Removes from the binary TIFF gave as input, groups of pixel with value 1, that have less pixels than a specified value.
 * @param   input: A binary TIFF to be processed.
 * @param   output: A binary TIFF.
 * @param   groupSize: The inferior limit of the group of pixels size.
 */
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

/**
 * @brief  Computes the HO.
 * @param  net_radiation_line[]: Array containing the specified line from the Rn computation.
 * @param  soil_heat_flux[]: Array containing the specified line from the G computation.
 * @param  width_band: Band width.
 * @param  ho_line[]: Auxiliary array for save the calculated value of HO for the line.
 */
void hoCalc(double net_radiation_line[], double soil_heat_flux[], int width_band, double ho_line[]){

    for(int col = 0; col < width_band; col++)
        ho_line[col] = net_radiation_line[col] - soil_heat_flux[col];

};

/**
 * @brief   Select a pair of pixels, which will be used as hot and cold pixels of the SEBAL Algorithm.
 * @note    For further information, check https://www.sciencedirect.com/science/article/pii/S0034425717302018.
 * @param   ndvi: NDVI TIFF.
 * @param   surface_temperature: TS TIFF.
 * @param   albedo: Albedo TIFF.
 * @param   net_radiation: Rn TIFF.
 * @param   soil_heat: G TIFF.
 * @param   landCover: Land Cover TIFF.
 * @param   height_band: Band height.
 * @param   width_band: Band width.
 * @param   output_path: A path where will be written auxiliary TIFFs generated in the process.
 * @retval  A pair struct of pixel, where the first one is the hot pixel selected, and the second is the cold one.
 */
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
    vector<Candidate> listTS;

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

            if(all_condition[col] && !isnan(ndvi_line[col])){

                listTS.push_back(Candidate(ndvi_line[col],
                                           surface_temperature_line[col],
                                           net_radiation_line[col],
                                           soil_heat_line[col],
                                           ho_line[col],
                                           line, col));

            }

        }

    }

    vector<Candidate> listNDVI (listTS);

    sort(listTS.begin(), listTS.end(), compare_candidate_temperature);
    sort(listNDVI.begin(), listNDVI.end(), compare_candidate_ndvi);

    double ts_min = listTS[0].temperature, ts_max = listTS[listTS.size() - 1].temperature;
    double ndvi_min = listNDVI[0].ndvi, ndvi_max = listNDVI[listNDVI.size() - 1].ndvi;
    int binTS = int(ceil((ts_max - ts_min)/0.25)); //0.25 is TS bin size
    int binNDVI = int(ceil((ndvi_max - ndvi_min)/0.01)); //0.01 is ndvi bin size

    vector<Candidate> histTS[binTS], final_histTS;

    for(Candidate c : listTS) {
        
        int pos = int(ceil((c.temperature - ts_min)/0.25));
        histTS[pos > 0 ? pos-1 : 0].push_back(c);

    }

    for(int i = 0; i < binTS; i++) {

        if(histTS[i].size() > 50) {

            for(Candidate c : histTS[i])
                final_histTS.push_back(c);

        }

    }

    vector<Candidate> histNDVI[binNDVI], final_histNDVI;
    for(Candidate c : listNDVI) {
        
        int pos = int(ceil((c.ndvi - ndvi_min)/0.01));
        histNDVI[pos > 0 ? pos-1 : 0].push_back(c);

    }

    for(int i = 0; i < binNDVI; i++) {

        if(histNDVI[i].size() > 50) {

            for(Candidate c : histNDVI[i])
                final_histNDVI.push_back(c);

        }

    }

    // Select cold pixel
    int pixel_count = 0, n1 = 1, n2 = 1, ts_pos, ndvi_pos, beginTs = 0, beginNDVI = final_histNDVI.size() - 1;
    vector<Candidate> coldPixels;
    while (pixel_count < 10 && !(n2 == 10 && n1 == 10)) {

        ts_pos = int(floor(n1/100.0 * final_histTS.size()));
        ndvi_pos = int(floor((100 - n2)/100.0 * final_histNDVI.size()));

        for(int i = beginTs; i <= ts_pos && pixel_count < 10; i++) {

            for(int j = beginNDVI; j >= ndvi_pos && pixel_count < 10; j--) {

                if(equals(final_histTS[i], final_histNDVI[j])){
                    
                    coldPixels.push_back(final_histTS[i]);
                    pixel_count++;

                }

            }

        }
	
	    beginTs = ts_pos;
	    beginNDVI = ndvi_pos;

        if(n2 < 10) n2++;
        else if(n1 < 10){
            n1++;
            beginNDVI = final_histNDVI.size() - 1;
        }

    }

    //Select hot pixel
    pixel_count = 0, n1 = 1, n2 = 1;
    vector<Candidate> hotPixels;
    beginTs = final_histTS.size() - 1, beginNDVI = 0;
    while (pixel_count < 10 && !(n2 == 10 && n1 == 10)) {
        
        ts_pos = int(floor((100 - n1)/100.0 * final_histTS.size()));
        ndvi_pos = int(floor(n2/100.0 * final_histNDVI.size()));

        for(int i = beginNDVI; i <= ndvi_pos && pixel_count < 10; i++) {

            for(int j = beginTs; j >= ts_pos && pixel_count < 10; j--) {

                if(equals(final_histTS[j], final_histNDVI[i])){

                    hotPixels.push_back(final_histTS[j]);
                    pixel_count++;

                }

            }

        }

	    beginTs = ts_pos;
	    beginNDVI = ndvi_pos;

        if(n2 < 10) n2++;
        else if(n1 < 10){
            n1++;
            beginTs = final_histTS.size() - 1;
        }

    }

    sort(coldPixels.begin(), coldPixels.end(), compare_candidate_ndvi);
    sort(hotPixels.begin(), hotPixels.end(), compare_candidate_temperature);

    return {hotPixels[hotPixels.size() - 1], coldPixels[coldPixels.size() - 1]};
   
}
