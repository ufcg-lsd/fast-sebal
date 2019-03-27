#include "candidate.h"

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

void Candidate::setAerodynamicResistance(double u200, double A_ZOM, double B_ZOM, double VON_KARMAN){
    this->zom = exp(A_ZOM + B_ZOM * this->ndvi);
    this->ustar = (VON_KARMAN * u200)/log(200/this->zom);
    this->aerodynamic_resistance.push_back(log(20)/(this->ustar * VON_KARMAN));
}

void Candidate::extract_negative_neighbour(TIFF *ndvi){
    double pixel_value;
    int cont = 0;
    for(int i = -3; i <= 3; i++){
        for(int j = -3; j <= 3; j++){
            pixel_value = read_position_tiff(ndvi, this->col + i, this->line + j);
            printf("Pixel value %d: %.8lf\n", cont, pixel_value);
            cont++;
            if(!isnan(pixel_value) && pixel_value < 0)
                this->negative_neighbour++;
        }
    }
    //printf("Negative neighbour: %d\n", this->negative_neighbour);
}

void Candidate::extract_coefficient_variation(TIFF *ndvi){
    vector<double> values_pixels_neighbours;
    double pixel_value;
    int cont = 0;
    for(int i = -3; i <= 3; i++){
        for(int j = -3; j <= 3; j++){
            pixel_value = read_position_tiff(ndvi, this->col + i, this->line + j);
            printf("Pixel value %d: %.8lf\n", cont, pixel_value); //DEBUG
            cont++;
            if(!isnan(pixel_value))
                values_pixels_neighbours.push_back(pixel_value);
        }
    }

    double mean, sd;
    double sum = 0;

    for(int i = 0; i < values_pixels_neighbours.size(); i++)
        sum += values_pixels_neighbours[i];
    
    mean = sum / values_pixels_neighbours.size();

    sum = 0;
    for(int i = 0; i < values_pixels_neighbours.size(); i++)
        sum += (values_pixels_neighbours[i] - mean) * (values_pixels_neighbours[i] - mean);

    sd = sqrt(sum / values_pixels_neighbours.size());

    this->coefficient_variation = sd / mean;
}

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
    
    if(this->aerodynamic_resistance.size() > 0){
        printf("Aerodynamic resistance:\n");
        for(unsigned i = 0; i < this->aerodynamic_resistance.size(); i++){
            printf("Aerodynamic resistance[%i]: %.10f\n", i, this->aerodynamic_resistance[i]);
        }
    }

    printf("\n");
}

bool compare_candidate_temperature(Candidate a, Candidate b){
    return a.temperature < b.temperature;
}

bool compare_candidate_ho(Candidate a, Candidate b){
    return a.ho < b.ho;
}