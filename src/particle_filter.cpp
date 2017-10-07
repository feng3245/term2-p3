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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "map.h"
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
default_random_engine gen;
num_particles = 8;
particles = std::vector<Particle>(num_particles);
int i = 0;

normal_distribution<double> dist_x(x, std[0]);

normal_distribution<double> dist_y(y, std[1]);

normal_distribution<double> dist_theta(theta, std[2]);



while(i<num_particles)
{

Particle p = Particle();

p.id = i;
p.x = dist_x(gen);
p.y = dist_y(gen);
p.theta = dist_theta(gen);
p.weight = 1.0;
particles.at(i) = p;
i++;
}
is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
int i = 0;
default_random_engine gen;


while(i<num_particles)
{

double x, y, theta;

Particle p = particles.at(i);
theta = p.theta;
double thetanorm = p.theta;

if(fabs(yaw_rate)>0.0001)
{
x = p.x+velocity/yaw_rate*(sin(thetanorm+yaw_rate*delta_t)-sin(thetanorm));
y = p.y+velocity/yaw_rate*(cos(thetanorm)-cos(thetanorm+yaw_rate*delta_t));
}
else
{
x = p.x + velocity*cos(theta+yaw_rate)*delta_t;
y = p.y + velocity*sin(theta+yaw_rate)*delta_t;
}
theta = theta+yaw_rate*delta_t;
normal_distribution<double> dist_x(x, std_pos[0]);

normal_distribution<double> dist_y(y, std_pos[1]);

normal_distribution<double> dist_theta(theta, std_pos[2]);


p.x = dist_x(gen);
p.y = dist_y(gen);
p.theta = dist_theta(gen);
particles.at(i) = p;

i++;
}


}

void ParticleFilter::dataAssociation(std::vector<Map::single_landmark_s> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
int i =0;
for(i=0;i<observations.size();i++)
{
LandmarkObs observation = observations.at(i);
double distance = -1;
int j=0;
double predictedx, predictedy;
for(j=0; j<predicted.size();j++)
{
//std::cout << "observation x "<< observation.x << " observation y "<< observation.y << std::endl;
//std::cout << "pred x " << predicted.at(j).x_f << " pred y "<< predicted.at(j).y_f << endl;
//std::cout << "predicted id " << predicted.at(j).id_i << endl;
double tempdistance = pow(observations.at(i).x - predicted.at(j).x_f,2)+pow(observations.at(i).y - predicted.at(j).y_f,2);
if(distance <0 || tempdistance < distance)
{
predictedx = predicted.at(j).x_f;
predictedy = predicted.at(j).y_f;
distance = tempdistance;
observation.id = predicted.at(j).id_i;
//std::cout << "distance is " << distance << "id is "<< observation.id <<endl;

}
}
observations.at(i) = observation;
//std::cout << "obs id " << observations.at(i).id << endl;
//cout << " assocations Landmark " << observation.id << " obs x,y "<< observation.x << observation.y << " landmark x,y " <<predictedx << " " << predictedy << endl << endl;
 
}
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
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

// Transformation
int j = 0;
int i = 0;

for(j =0 ;j < num_particles; j++)
{
std::vector<LandmarkObs> translatedObservations = std::vector<LandmarkObs>(); 
std::vector<double> sense_x = vector<double>();
std::vector<double> sense_y = vector<double>();
for(i =0 ;i < observations.size(); i++ )
{
double theta = particles.at(j).theta;
double xpart = particles.at(j).x;
double ypart = particles.at(j).y;
double obsx = observations.at(i).x;
double obsy = observations.at(i).y;
double x_map = xpart+ (cos(theta)*obsx)-(sin(theta)*obsy);
double y_map = ypart +(sin(theta)*obsx)+(cos(theta)*obsy);
sense_x.push_back(x_map);
sense_y.push_back(y_map);
LandmarkObs translatedObservation = LandmarkObs();
translatedObservation.x = x_map;
translatedObservation.y = y_map;
translatedObservations.push_back(translatedObservation);
//cout << "Transform obs x,y "<< obsx << " " << obsy << " Tobs x,y "<< x_map << " " <<y_map << endl;
}
dataAssociation(map_landmarks.landmark_list,translatedObservations);
int k = 0;
std::vector<int> associations = std::vector<int>();
for(k= 0; k < translatedObservations.size();k++)
{
associations.push_back(translatedObservations.at(k).id);

}
particles.at(j)=SetAssociations(particles.at(j),associations, sense_x, sense_y);
}
//update weights
int l =0;
for(l =0; l< particles.size(); l++)
{
	Particle p = particles.at(l);
	int m = 0;
	double finalweight = 1.0;
	for(m = 0; m < p.sense_x.size(); m++)
	{
		double x_obs = p.sense_x.at(m);
		double y_obs = p.sense_y.at(m);
		double mu_x = 0;
		double mu_y = 0;
		int z = 0;
		for(z = 0; z<map_landmarks.landmark_list.size(); z++)
		{
			if(map_landmarks.landmark_list.at(z).id_i == p.associations.at(m))
			{
				mu_x = map_landmarks.landmark_list.at(z).x_f;
				mu_y = map_landmarks.landmark_list.at(z).y_f;
			}
		}
//		std::cout << "x obs "<< x_obs <<std::endl;
//		std::cout << "y obs "<< y_obs <<std::endl;
//		std::cout << "landmark x " << mu_x <<endl;
//		std::cout << "landmark y " << mu_y <<endl;
		double exponent = ((x_obs - mu_x)*(x_obs - mu_x))/(2 * std_landmark[0]*std_landmark[0]) + ((y_obs - mu_y)*(y_obs - mu_y))/(2 * std_landmark[1]*std_landmark[1]);
//		std::cout << "exp " << exponent << std::endl;
		finalweight *= (1.0/(2*M_PI*std_landmark[0]*std_landmark[1]))*exp(-exponent);
//std::cout << "associations "<< getAssociations(p) << endl;

	}
	p.weight = finalweight;
	particles[l] = p;
//std::cout << "final weight is" << finalweight <<endl;
}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
std::vector<Particle> resampledParticles = std::vector<Particle>(num_particles);
int i =0;
double totalprob = 0;

for(i =0; i< num_particles; i++)
{
	totalprob += particles.at(i).weight;
}
std::vector<double> wnorm = std::vector<double>(num_particles);
int j = 0;
for(j = 0; j< num_particles; j++)
	wnorm.at(j) = particles.at(j).weight/totalprob;
int k =0;
double wheelend =0;
std::vector<double> wheeledw = std::vector<double>(num_particles);
for (k = 0; k < num_particles; k++)
{
	wheelend = wheelend+wnorm.at(k);
	wheeledw.at(k) = wheelend;
	resampledParticles.at(k) = particles.at(k);
}
int l = 0;
for(l = 0; l < num_particles; l++)
{
double prob = ((double) rand()) / (double) RAND_MAX;
int m = 0;
while(wheeledw.at(m)<prob)
m++;

particles.at(l) = resampledParticles.at(m);
}
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
