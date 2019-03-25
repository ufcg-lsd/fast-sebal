#pragma once

#include "utils.h"

struct Candidate{
	double ndvi, temperature, ustar;
	double net_radiation, soil_heat_flux, ho, zom;
	int line, col;
	int negative_neighbour;
	double coefficient_variation;

	vector< double > aerodynamic_resistance;

	Candidate();
	Candidate(double ndvi, double temperature, double net_radiation, double soil_heat_flux, double ho, int line, int col);
	void setAerodynamicResistance(double u200, double A_ZOM, double B_ZOM, double VON_KARMAN);
	void extract_negative_neighbour(TIFF *ndvi);
	void extract_coefficient_variation(TIFF *ndvi);
	void toString();
};

bool compare_candidate_temperature(Candidate a, Candidate b);

bool compare_candidate_ho(Candidate a, Candidate b);