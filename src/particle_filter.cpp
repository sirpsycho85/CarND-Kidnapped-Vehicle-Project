/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <cmath>
#include "particle_filter.h"

using namespace std;

default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	
	num_particles = 10;
	particles.resize(num_particles);
	weights.resize(num_particles);

	double std_x, std_y, std_theta;
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];

	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for (int i = 0; i < num_particles; ++i) {
		struct Particle temp_particle;
		temp_particle.id = i;
		temp_particle.x = dist_x(gen);
		temp_particle.y = dist_y(gen);
		temp_particle.theta = dist_theta(gen);
		temp_particle.weight = 1;

		particles[i] = temp_particle;
		weights[i] = temp_particle.weight;
	}

	is_initialized = true;

	//TODO: do I need a new struct each time? how assignment works...
	//TODO: do I need to create new default_random_engine each time?
	// having an separate array of weights helps you pass that as a parameter to the function discrete_distribution
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

	/*cout << "v = " << velocity << "\tdt = " << delta_t
		<< "\tx = " << particles[1].x
		<< "\ty = " << particles[1].y
		<< endl;*/

	for(int i = 0; i < num_particles; ++i) {
		double x, y, theta;
		x = particles[i].x;
		y = particles[i].y;
		theta = particles[i].theta;

		particles[i].x += velocity/yaw_rate
			* (sin(theta + yaw_rate*delta_t) - sin(theta));

		particles[i].y += velocity/yaw_rate
			* (cos(theta) - cos(theta + yaw_rate*delta_t));

		particles[i].theta += yaw_rate*delta_t;

		normal_distribution<double> dist_x(0, std_pos[0]);
		normal_distribution<double> dist_y(0, std_pos[1]);
		normal_distribution<double> dist_theta(0, std_pos[2]);

		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);

		//TODO why measurement uncertainties, not process uncertanties. Reconcile with Kalman.

		//TODO when this is called from main why does it say "noiseless"

		//TODO why is error always the same, even after clean and recompile.
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	/*
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	predicted is distances from one particle's perspective to all map landmarks within range
		so the landmark IDs are known
	observations is actual measurements from lidar
	perform nearest neighbor with those inputs
		so you can assign the landmark ID to the actual measurement

	for every actual measurement from lidar
		find the predicted that's closest
		assign the landmark id from the predicted to the actual measurement

	*/

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {

	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
