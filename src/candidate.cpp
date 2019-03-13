#include "candidate.h"

Candidate::Candidate(){
    this->ndvi = 0;
    this->temperature = 0;
    this->net_radiation = 0;
    this->soil_heat_flux = 0;
    this->ho = 0;
}

Candidate::Candidate(double ndvi, double temperature, double net_radiation, double soil_heat_flux, double ho){
    this->ndvi = ndvi;
    this->temperature = temperature;
    this->net_radiation = net_radiation;
    this->soil_heat_flux = soil_heat_flux;
    this->ho = ho;
}

void Candidate::setAerodynamicResistance(double u200, double A_ZOM, double B_ZOM, double VON_KARMAN){
    this->zom = exp(A_ZOM + B_ZOM * this->ndvi);
    this->ustar = (VON_KARMAN * u200)/log(200/this->zom);
    this->aerodynamic_resistance.push_back(log(2/0.1)/(this->ustar * VON_KARMAN));
}

bool compare_candidate_temperature(Candidate a, Candidate b){
    return a.temperature < b.temperature;
}

bool compare_candidate_ho(Candidate a, Candidate b){
    return a.ho < b.ho;
}