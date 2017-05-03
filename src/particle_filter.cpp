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
#include "Eigen/Dense"
#include "particle_filter.h"

using namespace std;

// declaring some functions up top bc don't want to touch .h file for autograder
vector<Map::single_landmark_s> GetLandmarksWithinRange(struct Particle particle, vector<Map::single_landmark_s> landmark_list, double sensor_range);
vector<LandmarkObs> CastLandmarksAsObservations(vector<Map::single_landmark_s> landmark_list);
void UpdateParticleWeights(struct Particle particle, vector<LandmarkObs> converted_landmarks_observations, vector<LandmarkObs> associated_observations);
LandmarkObs GetLandmarkObservationById(vector<LandmarkObs> landmark_observations, int id);

default_random_engine gen(std::random_device{}());

struct normal_random_variable {
	// Credit: http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
    normal_random_variable(Eigen::MatrixXd const& covar)
        : normal_random_variable(Eigen::VectorXd::Zero(covar.rows()), covar)
    {}

    normal_random_variable(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
        : mean(mean)
    {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
        transform = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::VectorXd mean;
    Eigen::MatrixXd transform;

    Eigen::VectorXd operator()() const
    {
        static std::mt19937 gen{ std::random_device{}() };
        static std::normal_distribution<> dist;

        return mean + transform * Eigen::VectorXd{ mean.size() }.unaryExpr([&](double x) { return dist(gen); });
    }
};

void ParticleFilter::init(double x, double y, double theta, double std[]) {

	bool verbose = false;

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

		if(verbose) {
			cout << "Particle " << i 
			<< " x = " << temp_particle.x << " y = " << temp_particle.y 
			<< " theta = " << temp_particle.theta
			<< endl;}
	}

	is_initialized = true;

	//TODO: initialize to 1 or 1/num_particles?
	// having an separate array of weights helps you pass that as a parameter to the function discrete_distribution
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

	bool verbose = false;

	for(int i = 0; i < num_particles; ++i) {
		
		double theta = particles[i].theta;

		if(verbose){cout << "x = " << particles[i].x
			<< " y = " << particles[i].y
			<< " theta = " << particles[i].theta
			<< " v*t = " << delta_t *velocity
			<< endl;}

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

		if(verbose){cout << "x = " << particles[i].x
			<< " y = " << particles[i].y
			<< endl;}

		//TODO why measurement uncertainties, not process uncertanties. Reconcile with Kalman.

		//TODO when this is called from main why does it say "noiseless"
	}
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		vector<LandmarkObs> observations, Map map_landmarks) {

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
	
	//Other notes
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

		//make a copy of the observations
			//for each landmark
			//  if within range
			//    convert to particle coordianate system
			//    dataAssiciation: update copy of observations with landmark IDs
			//translation

	bool verbose = false;
	for(int p_num = 0; p_num < num_particles; ++p_num) {

		if(verbose) {cout<<"particle: " << p_num<<endl;}
		
		vector<Map::single_landmark_s> landmarks_in_range, converted_landmarks;
		vector<LandmarkObs> converted_landmarks_observations, associated_observations;

		landmarks_in_range = GetLandmarksWithinRange(particles[p_num], map_landmarks.landmark_list, sensor_range);

		converted_landmarks = ParticleFilter::ConvertToParticleCoordinates(particles[p_num], landmarks_in_range);

		converted_landmarks_observations = CastLandmarksAsObservations(converted_landmarks);

		associated_observations = observations;

		dataAssociation(converted_landmarks_observations, associated_observations);

		UpdateParticleWeights(particles[p_num],converted_landmarks_observations, associated_observations);

		if(verbose) {cout<<endl;}
	}
}


vector<Map::single_landmark_s> GetLandmarksWithinRange(struct Particle particle, vector<Map::single_landmark_s> landmark_list, double sensor_range) {
	bool verbose = false;

	vector<Map::single_landmark_s> landmarks_in_range;
	
	if(verbose) {cout << "landmarks in range: ";}

	for(int lm_num = 0; lm_num < landmark_list.size(); ++lm_num) {
		
		double lm_dist = dist(particle.x, particle.y, landmark_list[lm_num].x_f, landmark_list[lm_num].y_f);

		if(lm_dist < sensor_range) {
			landmarks_in_range.push_back(landmark_list[lm_num]);

			if(verbose) {cout << landmark_list[lm_num].id_i << " ";}
		}
	}

	if(verbose) {cout << endl;}
	return landmarks_in_range;
}


vector<Map::single_landmark_s> ParticleFilter::ConvertToParticleCoordinates(struct Particle particle, vector<Map::single_landmark_s> landmark_list) {
	
	//https://math.stackexchange.com/questions/65059/converting-between-two-frames-of-reference-given-two-points

	bool verbose = false;
	
	vector<Map::single_landmark_s> converted_landmarks;

	if(verbose) {cout << "landmarks in particle coordinates: ";}

	double p_x = particle.x;
	double p_y = particle.y;
	double p_theta = particle.theta;

	for(int lm_num = 0; lm_num < landmark_list.size(); ++lm_num) {
		double lm_x = landmark_list[lm_num].x_f;
		double lm_y = landmark_list[lm_num].y_f;

		if(verbose) {cout << "theta: " << p_theta << " old: " << lm_x << " " << lm_y;}
		
		double temp_x = (lm_x-p_x)*cos(p_theta) + (lm_y-p_y)*sin(p_theta);
		double temp_y = -1*(lm_x-p_x)*sin(p_theta) + (lm_y-p_y)*cos(p_theta);

		Map::single_landmark_s converted_landmark;
		converted_landmark.id_i = landmark_list[lm_num].id_i;
		converted_landmark.x_f = temp_x;
		converted_landmark.y_f = temp_y;

		converted_landmarks.push_back(converted_landmark);

		if(verbose) {cout << " new: " << temp_x << " " << temp_y;}
	}

	if(verbose) {cout<<endl;}

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

	bool verbose = false;
	if(verbose) {cout<<"data association:" << endl;}

	double min_dist;
	double current_dist;
	int nearest_landmark_id;
	double nearest_x, nearest_y;

	for(int i = 0; i < observations.size(); ++i) {

		if(verbose) {cout << "obs: " << i 
			<< " x: " << observations[i].x << " y: " << observations[i].y
			<< endl;}
		
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

		if(verbose) {cout<<" dist: " << min_dist
			<< " landmark: " << nearest_landmark_id
			<< " nearest x: " << nearest_x << " nearest y: " << nearest_y
			<< endl;}
	}
}


void UpdateParticleWeights(struct Particle particle, vector<LandmarkObs> converted_landmarks_observations, vector<LandmarkObs> associated_observations) {
	/*
	int size = 2;
	Eigen::MatrixXd covar(size,size);
	covar << 1, .5,
	        .5, 1;

	normal_random_variable sample { covar };

	std::cout << sample() << std::endl;
	std::cout << sample() << std::endl;
	*/

	// for every measurement, find the landmark from converted observations that matches ID.
	// then figure out likelihood and multiply it into the product
	// remember to update both partile.weight weights[]

	//double product_of_probabilities = 1;

	bool verbose = false;

	for(int i = 0; i < associated_observations.size(); ++i) {
		//get by id from converted_landmark_observations
		if(verbose) {
			LandmarkObs myLandmarkObs = GetLandmarkObservationById(converted_landmarks_observations,associated_observations[i].id);
			cout << myLandmarkObs.id << " ";}
	}
	if(verbose) {cout << endl;}
}


LandmarkObs GetLandmarkObservationById(vector<LandmarkObs> landmark_observations, int id) {
	return landmark_observations[id];
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
