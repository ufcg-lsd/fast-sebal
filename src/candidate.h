#pragma once

#include "types.h"

struct Candidate{
	double ndvi, temperature, ustar;
	double net_radiation, soil_heat_flux, ho, zom;

	vector< double > aerodynamic_resistance;

	Candidate();
	Candidate(double ndvi, double temperature, double net_radiation, double soil_heat_flux, double ho);
	void setAerodynamicResistance(double u200, double A_ZOM, double B_ZOM, double VON_KARMAN);
	string toString();
};

bool compare_candidate_temperature(Candidate a, Candidate b);

bool compare_candidate_ho(Candidate a, Candidate b);