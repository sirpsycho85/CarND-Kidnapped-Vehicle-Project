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
	now you have for each particle, a mapping between actual measurements and landmarks
		how close these are corresponds to the weight you give each particle

	for every actual measurement from lidar
		find the predicted that's closest
		assign the landmark id from the predicted to the actual measurement

	*/
	double min_dist, current_dist;
	int nearest_landmark_id;

	for(int i = 0; i < observations.size(); ++i) {
		
		min_dist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
		nearest_landmark_id = predicted[0].id;
		
		for(int j = 1; j < predicted.size(); ++j) {
			current_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if(current_dist < min_dist) {
				nearest_landmark_id = predicted[j].id;
			}
		}

		observations[i].id = nearest_landmark_id;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	
	std::vector<LandmarkObs> predicted;

	vector<Map::single_landmark_s> landmark_list = map_landmarks.landmark_list;

	for(int p_num = 0; p_num < num_particles; ++p_num) {
		//translation
		double p_x = particles[p_num].x;
		double p_y = particles[p_num].y;
		double p_theta = particles[p_num].x;

		for(int o_num = 0; o_num < observations.size(); ++o_num) {
			double o_x = observations[o_num].x;
			double o_y = observations[o_num].y;
			
			double temp_x = o_x*cos(p_theta) + o_y*sin(p_theta) + p_x;
			double temp_y = o_x*sin(p_theta) + o_y*cos(p_theta) + p_y;
			
			observations[o_num].x = temp_x;
			observations[o_num].y = temp_y;
		}

		//for each landmark...
		for(int lm_num = 0; lm_num < landmark_list.size(); ++lm_num) {
			double lm_dist = dist(p_x, p_y, landmark_list[lm_num].x_f, landmark_list[lm_num].y_f);
			if(lm_dist < sensor_range) {
				struct LandmarkObs predicted_landmark;

				//TODO - need to add measurements not absolute positions right?
				predicted_landmark.x = landmark_list[lm_num].x_f;
				predicted_landmark.y = landmark_list[lm_num].y_f;
				predicted_landmark.id = landmark_list[lm_num].id_i;
				
				predicted.push_back(predicted_landmark);
			}
		}
	}

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

	/*
	For each particle
		Transform observations into map coordiante system
		For each landmark
			Calculate distance to each actual landmark (predict measurements)
			If it's within sensor_range, add as a LandmarkObs to 'predicted'
		Call data association to update measurements with landmark IDs
		Update weights with multi-variate Gaussian
			Note: The same landmark within range may be the nearest for a measurement
			Note: Might be convenient in dataAssociation to create a new vector with actual landmarks
			That way you don't have to search through them again
		Normalize weights
	*/

	// TODO: handle empty 'predicted' or 'observations'
	// TODO: how do you know the map y points downward?
	// TODO: rotation matrix in article disagrees with wiki
	// TODO: not clear if the transformation is correct. See what happens to x for a measurement after translating, given some theta for a particle
	// 			also distances seem large
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
