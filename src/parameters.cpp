#include "parameters.h"

/**
 * @brief  Empty constructor. Setting all attributes to 0.
 */
MTL::MTL(){
    this->number_sensor = 0;
    this->julian_day = 0;
    this->year = 0;
    this->sun_elevation = 0;
    this->rad_add_10 = 0;
    this->rad_mult_10 = 0;
};

/**
 * @brief  Constructor receiving the metadata path.
 * @param  metadata_path: Metadata file path.
 */
MTL::MTL(string metadata_path){
    map<string, string> mtl;

    ifstream in(metadata_path);
    if(!in.is_open() || !in){
        cerr << "Open metadata problem!" << endl;
        exit(2);
    }

    string line;
    while(getline(in, line)){
        stringstream lineReader(line);
        string token;
        vector<string> nline;
        while(lineReader >> token)
            nline.push_back(token);
        
        if(nline.size() >= 3) {
            mtl[nline[0]] = nline[2];
        }
    }

    in.close();

    char julian_day[3];
    julian_day[0] = mtl["LANDSAT_SCENE_ID"][14];
    julian_day[1] = mtl["LANDSAT_SCENE_ID"][15];
    julian_day[2] = mtl["LANDSAT_SCENE_ID"][16];

    char year[4];
    year[0] = mtl["LANDSAT_SCENE_ID"][10];
    year[1] = mtl["LANDSAT_SCENE_ID"][11];
    year[2] = mtl["LANDSAT_SCENE_ID"][12];
    year[3] = mtl["LANDSAT_SCENE_ID"][13];

    int hours = atoi(mtl["SCENE_CENTER_TIME"].substr(1, 2).c_str());
    int minutes = atoi(mtl["SCENE_CENTER_TIME"].substr(4, 2).c_str());

    this->number_sensor = atoi(new char(mtl["LANDSAT_SCENE_ID"][3]));
    this->julian_day = atoi(julian_day);
    this->year = atoi(year);
    this->sun_elevation = atof(mtl["SUN_ELEVATION"].c_str());
    this->distance_earth_sun = atof(mtl["EARTH_SUN_DISTANCE"].c_str());
    this->image_hour = (hours + minutes / 60.0) * 100;

    if(this->number_sensor == 8){
        rad_mult_10 = atof(mtl["RADIANCE_MULT_BAND_10"].c_str());
        rad_add_10 = atof(mtl["RADIANCE_ADD_BAND_10"].c_str());
    }
};

/**
 * @brief  Constructor.
 * @param  number_sensor: Number of the satellite sensor.sensor_path
 * @param  year: Year of image.
 */
Sensor::Sensor(int number_sensor, int year){
    string sensor_path = capture_parameter_path(number_sensor, year);
    load_parameter_values(sensor_path);
}

/**
 * @brief  Get the path of the correct sensor parameters based on his number.
 * @param  number_sensor: Number of the satellite sensor.
 * @param  year: Year of image. 
 * @retval Path to sensor parameter file.
 */
string Sensor::capture_parameter_path(int number_sensor, int year){
    switch(number_sensor){
        case 8:
            return "./parameters/LC.data";
            break;
        case 7:
            return "./parameters/ETM.data";
            break;
        case 5:
            if(year < 1992)
                return "./parameters/TM1.data";
            else
                return "./parameters/TM2.data";
            break;
        default:
            cerr << "Sensor problem" << endl;
            exit(6);
    }
}

/**
 * @brief  Loads the sensor parameters.
 * @param  sensor_path: Path to sensor parameter file.
 */
void Sensor::load_parameter_values(string sensor_path){
    ifstream in(sensor_path);
    if(!in.is_open() || !in) {
        cerr << "Open sensor parameters problem" << endl;
        exit(1);
    }

    string line, token;
    for(int i = 1; i < 8; i++){
        getline(in, line);
        stringstream lineReader(line);
        vector<string> nline;
        while(lineReader >> token)
            nline.push_back(token);

        this->parameters[i][this->GRESCALE] = atof(nline[0].c_str());
        this->parameters[i][this->BRESCALE] = atof(nline[1].c_str());
        this->parameters[i][this->ESUN] = atof(nline[2].c_str());
        this->parameters[i][this->WB] = atof(nline[3].c_str());
    }
    
    in.close();
};

/**
 * @brief  Empty constructor. Set temperature_image to 0.
 */
Station::Station(){
    this->temperature_image = 0;
};

/**
 * @brief  Constructor.
 * @param  station_data_path: Weather station data file.
 * @param  image_hour: Image hour.
 */
Station::Station(string station_data_path, double image_hour){
    ifstream in(station_data_path);
    if(!in.is_open() || !in){
        cerr << "Open station data problem!" << endl;
        exit(2);
    }

    string line;
    while(getline(in, line)){
        istringstream lineReader(line);
        vector<string> nline;
        string token;
        while(getline(lineReader, token, ';'))
            nline.push_back(token);

        if(nline.size())
            this->info.push_back(nline);
    }

    in.close();

    if(this->info.size() < 1){
        cerr << "Station data empty!" << endl;
        exit(12);
    }

    double diff = fabs(atof(this->info[0][2].c_str()) - image_hour);
    this->temperature_image = atof(this->info[0][6].c_str());
    this->v6 = atof(this->info[0][5].c_str());
    this->v7_max = atof(this->info[0][6].c_str());
    this->v7_min = atof(this->info[0][6].c_str());
    this->latitude = atof(this->info[0][3].c_str());
    this->longitude = atof(this->info[0][4].c_str());

    for(int i = 1; i < this->info.size(); i++){

        v7_max = max(v7_max, atof(this->info[i][6].c_str()));
        v7_min = min(v7_min, atof(this->info[i][6].c_str()));

        if(fabs(atof(this->info[i][2].c_str()) - image_hour) < diff){
            diff = fabs(atof(this->info[i][2].c_str()) - image_hour);
            this->temperature_image = atof(this->info[i][6].c_str());
            this->v6 = atof(this->info[i][5].c_str());
        }
    }
};