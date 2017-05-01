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

void TestSimpleCoordianteTransform() {
	double xp, yp, theta_m_to_p, xlm, ylm;
	
	xp = 1;
	yp = 2;
	theta_m_to_p = M_PI*0.5;
	xlm = 2;
	ylm = 0;

	ParticleFilter::SimpleCoordianteTransform(xp, yp, theta_m_to_p, xlm, ylm);
}

int main() {

	TestConvertToParticleCoordinates();
	//TestSimpleCoordianteTransform();
	return 0;
}
