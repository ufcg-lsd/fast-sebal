#include "candidate.h"

/**
 * @brief  Empty constructor, all attributes are initialized with 0.
 */
Candidate::Candidate(){
    this->ndvi = 0;
    this->temperature = 0;
    this->net_radiation = 0;
    this->soil_heat_flux = 0;
    this->ho = 0;
    this->line = 0;
    this->col = 0;
    this->negative_neighbour = 0;
    this->coefficient_variation = 0;
    this->zom = 0;
    this->ustar = 0;
}

Candidate::Candidate(const Candidate& c){
    this->ndvi = c.ndvi;
    this->temperature = c.temperature;
    this->net_radiation = c.net_radiation;
    this->soil_heat_flux = c.soil_heat_flux;
    this->ho = c.ho;
    this->line = c.line;
    this->col = c.col;
    this->negative_neighbour = c.negative_neighbour;
    this->coefficient_variation = c.coefficient_variation;
    this->zom = c.zom;
    this->ustar = c.ustar;
};

/**
 * @brief  Constructor with initialization values to attributes.
 * @param  ndvi: Pixel's NDVI.
 * @param  temperature: Pixel's surface temperature.
 * @param  net_radiation: Pixel's net radiation.
 * @param  soil_heat_flux: Pixel's soil heat flux.
 * @param  ho: Pixel's ho.
 * @param  line: Pixel's line on TIFF.
 * @param  col: Pixel's column on TIFF.
 */
Candidate::Candidate(double ndvi, double temperature, double net_radiation, double soil_heat_flux, double ho, int line, int col){
    this->ndvi = ndvi;
    this->temperature = temperature;
    this->net_radiation = net_radiation;
    this->soil_heat_flux = soil_heat_flux;
    this->ho = ho;
    this->line = line;
    this->col = col;
    this->negative_neighbour = 0;
    this->coefficient_variation = 0;
    this->zom = 0;
    this->ustar = 0;
}

/**
 * @brief  Calculates a initial value for Pixel's aerodynamic resistance. Adding this value to attribute aerodynamic resistance.
 * @param  u200: Wind speed at 200 m.
 * @param  A_ZOM: Coefficient A.
 * @param  B_ZOM: Coefficient B.
 * @param  VON_KARMAN: Karman's constant.
 */
void Candidate::setAerodynamicResistance(double u200, double A_ZOM, double B_ZOM, double VON_KARMAN){
    this->zom = exp(A_ZOM + B_ZOM * this->ndvi);
    this->ustar = (VON_KARMAN * u200)/log(200/this->zom);
    this->aerodynamic_resistance.push_back(log(20)/(this->ustar * VON_KARMAN));
}

/**
 * @brief  Counts how many neighbors pixels of the Candidate pixel have a negative value of NDVI.
 * @note   A pixel is a neighbour pixel if it is in a radius of 105 meters from the Candidate pixel.
 * @param  *ndvi: NDVI TIFF.
 */
void Candidate::extract_negative_neighbour(TIFF *ndvi){

    uint32 height_band, width_band;
    TIFFGetField(ndvi, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(ndvi, TIFFTAG_IMAGEWIDTH, &width_band);

    double pixel_value;
    int cont = 1;
    for(int i = -3; i <= 2; i++){
        for(int j = -3; j <= 2; j++){

            if (this->col + i >= 0 && this->col + i < width_band && this->line + j >= 0 && this->line + j < height_band) {

                pixel_value = read_position_tiff(ndvi, this->col + i, this->line + j);
                cont++;
                if(!isnan(pixel_value) && pixel_value < 0)
                    this->negative_neighbour++;

            }

        }
    }

}

/**
 * @brief  Calculates the coefficient of variation (CV) of NDVI for neighbors pixels of the Candidate pixel.
 * @note   A pixel is a neighbour pixel if it is in a radius of 105 meters from the Candidate pixel.
 * @param  *ndvi: NDVI TIFF.
 */
void Candidate::extract_coefficient_variation(TIFF *ndvi){

    uint32 height_band, width_band;
    TIFFGetField(ndvi, TIFFTAG_IMAGELENGTH, &height_band);
    TIFFGetField(ndvi, TIFFTAG_IMAGEWIDTH, &width_band);

    vector<double> values_pixels_neighbors;
    double pixel_value;
    int cont = 1;
    for(int i = -3; i <= 2; i++){
        for(int j = -3; j <= 2; j++){

            if (this->col + i >= 0 && this->col + i < width_band && this->line + j >= 0 && this->line + j < height_band) {

                pixel_value = read_position_tiff(ndvi, this->col + i, this->line + j);
                cont++;
                if(!isnan(pixel_value))
                    values_pixels_neighbors.push_back(pixel_value);

            }
            
        }
    }
    
    double mean, sd;
    double sum = 0;

    for(int i = 0; i < values_pixels_neighbors.size(); i++)
        sum += values_pixels_neighbors[i];
    
    mean = sum / values_pixels_neighbors.size();

    sum = 0;
    for(int i = 0; i < values_pixels_neighbors.size(); i++)
        sum += (values_pixels_neighbors[i] - mean) * (values_pixels_neighbors[i] - mean);
    
    sd = sqrt(sum / values_pixels_neighbors.size());
    this->coefficient_variation = sd / mean;
}

/**
 * @brief  Prints the data contained at the struct.
 */
void Candidate::toString(){
    string toString;
    printf("NDVI: %.10lf\n", this->ndvi);
    printf("TS: %.10lf\n", this->temperature);
    printf("Rn: %.10lf\n", this->net_radiation);
    printf("Ustar: %.10lf\n", this->ustar);
    printf("Z0m: %.10lf\n", this->zom);
    printf("G: %.10lf\n", this->soil_heat_flux);
    printf("HO: %.10lf\n", this->ho);
    printf("Negative neighbour: %d\n", this->negative_neighbour);
    printf("Coefficient variation: %.10lf\n", this->coefficient_variation);
    printf("Line: %d\n", this->line);
    printf("Col: %d\n", this->col);
    
    if(this->aerodynamic_resistance.size() > 0){
        printf("Aerodynamic resistance:\n");
        for(unsigned i = 0; i < this->aerodynamic_resistance.size(); i++){
            printf("Aerodynamic resistance[%i]: %.10f\n", i, this->aerodynamic_resistance[i]);
        }
    }

    printf("\n");
}

bool equals(Candidate a, Candidate b){

    return (a.col == b.col) && (a.line == b.line);

}

/**
 * @brief  Compares two Candidates based upon their surface temperature.
 * @param  a: First candidate.
 * @param  b: Second candidate.
 * @retval TRUE if second candidate is greater than first one, and FALSE otherwise.
 */
bool compare_candidate_temperature(Candidate a, Candidate b){
    return a.temperature < b.temperature;
}

/**
 * @brief  Compares two Candidates based upon their surface temperature.
 * @param  a: First candidate.
 * @param  b: Second candidate.
 * @retval TRUE if second candidate is greater than first one, and FALSE otherwise.
 */
bool compare_candidate_ndvi(Candidate a, Candidate b){
    return a.ndvi < b.ndvi;
}

/**
 * @brief  Compares two Candidates based upon their HO.
 * @param  a: First candidate.
 * @param  b: Second candidate.
 * @retval TRUE if second candidate is greater than first one, and FALSE otherwise.
 */
bool compare_candidate_ho(Candidate a, Candidate b){
    return a.ho < b.ho;
}
