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

            if(pixel_value == AGR) { //Verify if the pixel passed the land cover test

                mask_line[column] = true;

                for(int i = -3; i <= 3 && mask_line[column]; i++){

                    for(int j = -3; j <= 3 && mask_line[column]; j++){

                        // Add for the NDVI, TS and Albedo the value of neighbors pixels into the respective vector

                        if (column + i >= 0 && column + i < width_band && line + j >= 0 && line + j < height_band) {

                            pixel_value = read_position_tiff(ndvi, column + i, line + j);
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

void testHomogeneity(TIFF* ndvi, TIFF* surface_temperature, TIFF* albedo, TIFF* mask){

    uint32 height_band, width_band;
    TIFFGetField(ndvi, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(ndvi, TIFFTAG_IMAGEWIDTH, &width_band);

    for(int line = 0; line < height_band; line++) {

        // Create the respective line of the binary map of eligibles pixels
        double mask_line[width_band];
        read_line_tiff(mask, mask_line, line);

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
                double meanNDVI, meanTS, menAlb;
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

        write_line_tiff(mask, mask_line, line);

    }

}