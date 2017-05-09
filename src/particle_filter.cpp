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
#include <map>
#include <math.h>
#include <list>
#include <algorithm>
#include "Eigen/Dense"
#include "particle_filter.h"

using namespace std;

//TODO: get original .h file and make appropriate changes

// declaring some functions up top bc don't want to touch .h file for autograder
vector<Map::single_landmark_s> GetLandmarksWithinRange(struct Particle particle, vector<Map::single_landmark_s> landmark_list, double sensor_range);
vector<LandmarkObs> CastLandmarksAsObservations(vector<Map::single_landmark_s> landmark_list);
void UpdateParticleWeights(struct Particle &particle, vector<LandmarkObs> converted_landmarks_observations, vector<LandmarkObs> associated_observations, double std_landmark[]);
LandmarkObs GetLandmarkObservationById(vector<LandmarkObs> landmark_observations, int id);
double MultivariatePDF(Eigen::VectorXd mean, Eigen::MatrixXd covar, Eigen::VectorXd measurement);
void NormalizeParticleWeights(vector<Particle> &particles, vector<double> &weights);
bool CompareByParticleWeights(struct Particle particleA, struct Particle particleB);

default_random_engine gen(std::random_device{}());

double largest_product_of_probabilities;

void ParticleFilter::init(double x, double y, double theta, double std[]) {

	num_particles = 5;
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

	//TODO: initialize to 1 or 1/num_particles?
	// having an separate array of weights helps you pass that as a parameter to the function discrete_distribution
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

	for(int i = 0; i < num_particles; ++i) {
		
		double theta = particles[i].theta;

		if(fabs(yaw_rate) < 0.0001) {
			yaw_rate = 0.0001;
		}

		particles[i].x += velocity/yaw_rate
			* (sin(theta + yaw_rate*delta_t) - sin(theta));

		particles[i].y += velocity/yaw_rate
			* (cos(theta) - cos(theta + yaw_rate*delta_t));

		particles[i].theta += yaw_rate*delta_t;

		//TODO: seems like the 

		normal_distribution<double> dist_x(0, std_pos[0]);
		normal_distribution<double> dist_y(0, std_pos[1]);
		normal_distribution<double> dist_theta(0, std_pos[2]);

		particles[i].x += dist_x(gen);
		particles[i].y += dist_y(gen);
		particles[i].theta += dist_theta(gen);

		//TODO why measurement uncertainties, not process uncertanties. Reconcile with Kalman.

		//TODO when this is called from main why does it say "noiseless"
	}
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], vector<LandmarkObs> observations, Map map_landmarks) {

	// TODO: are predicted observations to pass to multi-variate the distances or the absolute locations?
	// TODO: normalize weights
	// TODO: handle empty 'predicted' or 'observations'
	// TODO: how do you know the map y points downward?
	// TODO: rotation matrix in article disagrees with wiki...
	// Note: Might be convenient in dataAssociation to create a new vector with actual landmarks
	// 		That way you don't have to search through them again

	/*for each particle:
		landmarks_in_range = GetLandmarksWithingSensorRange(particle,landmarks,sensorRange)
		converted_landmarks = ConvertToParticleCoordinates(particle,remaining_landmarks)
		chosen_matched_observations = observations // just a copy for input into next function
		dataAssociation(converted_landmarks ,chosen_matched_observations )
		UpdateParticleWeights(particle,observations, chosen_matched_observations )
			http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c*/

	largest_product_of_probabilities = 0;

	for(int p_num = 0; p_num < num_particles; ++p_num) {
		
		vector<Map::single_landmark_s> landmarks_in_range, converted_landmarks;
		vector<LandmarkObs> converted_landmarks_observations, associated_observations;

		landmarks_in_range = GetLandmarksWithinRange(particles[p_num], map_landmarks.landmark_list, sensor_range);

		converted_landmarks = ParticleFilter::ConvertToParticleCoordinates(particles[p_num], landmarks_in_range);

		converted_landmarks_observations = CastLandmarksAsObservations(converted_landmarks);

		associated_observations = observations;

		dataAssociation(converted_landmarks_observations, associated_observations);
		UpdateParticleWeights(particles[p_num],converted_landmarks_observations, associated_observations, std_landmark);
		weights[p_num] = particles[p_num].weight;
	}

	NormalizeParticleWeights(particles, weights);
}


vector<Map::single_landmark_s> GetLandmarksWithinRange(struct Particle particle, vector<Map::single_landmark_s> landmark_list, double sensor_range) {

	vector<Map::single_landmark_s> landmarks_in_range;

	for(int lm_num = 0; lm_num < landmark_list.size(); ++lm_num) {
		
		double lm_dist = dist(particle.x, particle.y, landmark_list[lm_num].x_f, landmark_list[lm_num].y_f);

		if(lm_dist < sensor_range) {
			landmarks_in_range.push_back(landmark_list[lm_num]);
		}
	}

	return landmarks_in_range;
}


vector<Map::single_landmark_s> ParticleFilter::ConvertToParticleCoordinates(struct Particle particle, vector<Map::single_landmark_s> landmark_list) {
	
	//https://math.stackexchange.com/questions/65059/converting-between-two-frames-of-reference-given-two-points
	
	vector<Map::single_landmark_s> converted_landmarks;

	double p_x = particle.x;
	double p_y = particle.y;
	double p_theta = particle.theta;

	for(int lm_num = 0; lm_num < landmark_list.size(); ++lm_num) {
		double lm_x = landmark_list[lm_num].x_f;
		double lm_y = landmark_list[lm_num].y_f;
		
		double temp_x = (lm_x-p_x)*cos(p_theta) + (lm_y-p_y)*sin(p_theta);
		double temp_y = -1*(lm_x-p_x)*sin(p_theta) + (lm_y-p_y)*cos(p_theta);

		Map::single_landmark_s converted_landmark;
		converted_landmark.id_i = landmark_list[lm_num].id_i;
		converted_landmark.x_f = temp_x;
		converted_landmark.y_f = temp_y;

		converted_landmarks.push_back(converted_landmark);
	}

	return converted_landmarks;
}


vector<LandmarkObs> CastLandmarksAsObservations(vector<Map::single_landmark_s> landmark_list) {
	
	vector<LandmarkObs> observations;
	
	for(int lm_num = 0; lm_num < landmark_list.size(); ++lm_num) {
		LandmarkObs obs;
		obs.x = landmark_list[lm_num].x_f;
		obs.y = landmark_list[lm_num].y_f;
		obs.id = landmark_list[lm_num].id_i;
		observations.push_back(obs);
	}

	return observations;
}


void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
	
	//TODO: why are associations always the same

	double min_dist;
	double current_dist;
	int nearest_landmark_id;
	double nearest_x, nearest_y;

	for(int i = 0; i < observations.size(); ++i) {
		
		min_dist = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y);
		nearest_landmark_id = predicted[0].id;
		nearest_x = predicted[0].x;
		nearest_y = predicted[0].y;
		
		for(int j = 1; j < predicted.size(); ++j) {
			
			current_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			
			if(current_dist < min_dist) {
				nearest_landmark_id = predicted[j].id;
				nearest_x = predicted[j].x;
				nearest_y = predicted[j].y;
				min_dist = current_dist;
			}
		}

		observations[i].id = nearest_landmark_id;
	}
}


void UpdateParticleWeights(struct Particle &particle, 
	vector<LandmarkObs> converted_landmarks_observations, 
	vector<LandmarkObs> associated_observations,
	double std_landmark[]) {
	

	// for each in associated_observations
	// find the landmark from converted_landmarks_observations that matches ID2
	// then figure out likelihood and multiply it into the product
	// remember to update both partile.weight and weights[].

	//TODO: largest_product_of_probabilities which is product for "best" particle should get larger with more particles...

	double new_weight = 1;

	for(int i = 0; i < associated_observations.size(); ++i) {

		LandmarkObs myLandmarkObs = GetLandmarkObservationById(converted_landmarks_observations,
			associated_observations[i].id);

		Eigen::VectorXd mean = Eigen::VectorXd(2);
		Eigen::MatrixXd covar = Eigen::MatrixXd(2,2);
		Eigen::VectorXd x = Eigen::VectorXd(2);

		mean << myLandmarkObs.x,myLandmarkObs.y;
		covar << std_landmark[0],0,0,std_landmark[1];
		x << associated_observations[i].x,associated_observations[i].y;

		double prob =  MultivariatePDF(mean, covar, x);
		new_weight = new_weight * prob;
	}

	particle.weight = new_weight;

	if(new_weight > largest_product_of_probabilities) {
			largest_product_of_probabilities = new_weight;
	}
}


LandmarkObs GetLandmarkObservationById(vector<LandmarkObs> landmark_observations, int id) {

	for(int i = 0; i < landmark_observations.size(); ++i) {
		if(landmark_observations[i].id == id) {
			return landmark_observations[i];
		}
	}
}


double MultivariatePDF(Eigen::VectorXd mean, Eigen::MatrixXd covar, Eigen::VectorXd measurement) {
	double term1, term2;

	term1 = exp(-0.5 * (measurement - mean).transpose() * covar.inverse() * (measurement - mean));
	term2 = sqrt((2*M_PI*covar).determinant());
	
	return term1/term2;
}


void NormalizeParticleWeights(vector<Particle> &particles, vector<double> &weights) {
	double sum_of_weights = 0;
	
	for(int p_num = 0; p_num < particles.size(); ++p_num) {
		sum_of_weights = sum_of_weights + weights[p_num];
	}

	for(int p_num = 0; p_num < particles.size(); ++p_num) {
		particles[p_num].weight = particles[p_num].weight/sum_of_weights;
		weights[p_num] = weights[p_num]/sum_of_weights;
	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//BACKLOG: maps, iterators. Shouldn't need this intermediate new_weights vector.

    discrete_distribution<> d(weights.begin(), weights.end());

    map<int, int> m;
    for(int n=0; n<100; ++n) {
        ++m[d(gen)];
    }

    vector<double> new_weights;
    for(auto p : m) {
        new_weights.push_back(p.second);
    }

    vector<Particle> new_particles = {};
    for(int i = 0; i < particles.size(); ++i) {
    	struct Particle old_p = particles[i];
    	struct Particle new_p;
		new_p.id = old_p.id;
		new_p.x = old_p.x;
		new_p.y = old_p.y;
		new_p.theta = old_p.theta;
		new_p.weight = new_weights[i];
    }

    particles = new_particles;
    weights = new_weights;
}


bool CompareByParticleWeights(struct Particle particleA, struct Particle particleB) {
	return particleA.weight > particleB.weight;
}


void ParticleFilter::write(string filename) {
	// You don't need to modify this file.
	ofstream dataFile;
	dataFile.open(filename, ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
