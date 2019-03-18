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

/*
    double ndvi, temperature, ustar;
	double net_radiation, soil_heat_flux, ho, zom;

	vector< double > aerodynamic_resistance;
*/

void Candidate::toString(){
    string toString;
    printf("Candidate\n");
    printf("NDVI: %.2f\n", this->ndvi);
    printf("TS: %.2f\n", this->temperature);
    printf("RS: %.2f\n", this->net_radiation);
    printf("Ustar: %.2f\n", this->soil_heat_flux);
    printf("Z0m: %.2f\n", this->zom);
    printf("G: %.2f\n", this->soil_heat_flux);
    printf("HO: %.2f\n", this->ho);
    
    if(this->aerodynamic_resistance.size() > 0){
        printf("Aerodynamic resistance:\n");
        for(unsigned i = 0; i < this->aerodynamic_resistance.size(); i++){
            printf("Aerodynamic resistance[%i]: %.2f\n", i, this->aerodynamic_resistance[i]);
        }
    }else{
        cout << "Else do toString" << endl;
    }

    printf("\n");   
}

bool compare_candidate_temperature(Candidate a, Candidate b){
    return a.temperature < b.temperature;
}

bool compare_candidate_ho(Candidate a, Candidate b){
    return a.ho < b.ho;
}