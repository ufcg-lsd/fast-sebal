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

using namespace std;

const double EPS = 1e-7;
const double NaN = -sqrt(-1.0);
const double PI = acos(-1);
const double VON_KARMAN = 0.41;
const double GRAVITY = 9.81;
const double RHO = 1.15; //Air density
const double SPECIFIC_HEAT_AIR = 1004;
const double GSC = 0.082; //Solar constant