#pragma once

#include <tiffio.h>
#include <math.h>
#include <iostream>
#include <string.h>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include <time.h>
#include <chrono>
#include <stdlib.h>

using namespace std;

// CONSTANTS DECLARATION

// Epsilon
const double EPS = 1e-7;

// Not a number
const double NaN = -sqrt(-1.0);

// Pi
const double PI = acos(-1);

// Karman's constant
const double VON_KARMAN = 0.41;

// Earth's gravity
const double GRAVITY = 9.81;

// Atmospheric density
const double RHO = 1.15;

// Specific heat of air
const double SPECIFIC_HEAT_AIR = 1004;

// Solar constant
const double GSC = 0.082;
