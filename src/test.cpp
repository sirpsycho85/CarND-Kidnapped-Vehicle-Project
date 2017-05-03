#include <iostream>
#include <math.h>
#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

void TestConvertToParticleCoordinates() {
	struct Particle particle;
	particle.id = 0;
	particle.x = -1;
	particle.y = -1;
	particle.theta = M_PI*0.5;
	particle.weight = 1;

	Map map;
	if (!read_map_data("data/map_2.txt", map)) {
		cout << "Error: Could not open map file" << endl;
	}
	
	vector<Map::single_landmark_s> landmark_list = map.landmark_list;

	for(int i=0; i < landmark_list.size(); ++i) {
		cout<< "x=" << landmark_list[i].x_f
			<< " y=" << landmark_list[i].y_f
			<< endl;
	}

	vector<Map::single_landmark_s> converted_landmarks;

	ParticleFilter pf;

	converted_landmarks = pf.ConvertToParticleCoordinates(particle,landmark_list);
	for(int i=0; i < converted_landmarks.size(); ++i) {
		cout<< "x=" << converted_landmarks[i].x_f
			<< " y=" << converted_landmarks[i].y_f
			<< endl;
	}
}

void TestPrediction() {
	double x_init = 0;
	double y_init = 0;
	double theta_init = M_PI/2;
	double sigma_pos [3] = {0.1, 0.1, 0.0};

	ParticleFilter pf;
	pf.init(x_init, y_init, theta_init, sigma_pos);

	double delta_t = 1;
	double velocity = 1;
	double yawrate = 0;
	pf.prediction(delta_t, sigma_pos, velocity, yawrate);
}

int main() {

	//TestConvertToParticleCoordinates();
	//TestPrediction();
	//Tested association looking at one particle one timestep for a few LMs
	//Tested landmarks within range in excel sheet
	return 0;
}
