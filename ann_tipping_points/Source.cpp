#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>

#include "individual.h"

const int popsize = 5000;
std::vector<Individual> pop(5000); // population size: 5000
//Varied Parameters
// A and B were varied between 0 and 1 (A+B = 1) to investigate deterministic vs stochastic env., not found to have a great effect 

float A = 1.f;   //Deterministic scaling constant, default value 
float B = 0.f;   //stochastic scaling constant, default value

//R between 1 and 100000, P between 0 and 1
float R = 10.f;  //Environmental variation
float P = 0.99f;   //predictability
///
const auto& env_dist = std::uniform_real_distribution<float>(-1.f, 1.f); // not explicitly stated in botero 2015
float E;
float Cue;

int gmax = 50000;
int tmax = 5;

// mutation prob
std::bernoulli_distribution mut_event(0.001); // mutation probability
// mutation size
std::cauchy_distribution<double> m_shift(0.0, 0.01); // how much of mutation

// reproduction
void reproduction()
{
	// make fitness vec
	std::vector<double> fitness_vec;
	float max = 0.f; float min = 0.f;
	for (int a = 0; a < popsize; a++)
	{
		fitness_vec.push_back(1.0 / static_cast<double> (pop[a].mismatch));
	}

	// make temp pop vector, position and energy vectors
	std::vector<Individual> pop2(popsize);

	// assign parents
	for (int a = 0; a < popsize; a++) {

		std::discrete_distribution<> weighted_lottery(fitness_vec.begin(), fitness_vec.end());
		int parent_id = weighted_lottery(rnd::reng);

		// replicate ann_dev
		pop2[a].ann_dev = pop[parent_id].ann_dev;
		pop2[a].ann_life = pop[parent_id].ann_life;

		// reset mismatch
		pop2[a].mismatch = 0.f;
		// process a baseline -- this is development
		pop2[a].update_I_g(Cue);

		// mutate ann_dev
		for (auto& w : pop2[a].ann_dev) {
			// probabilistic mutation of ANN
			if (mut_event(rnd::reng)) {
				w += static_cast<float> (m_shift(rnd::reng));
			}
		}

		// mutate ann_life
		for (auto& w : pop2[a].ann_life) {
			// probabilistic mutation of ANN
			if (mut_event(rnd::reng)) {
				w += static_cast<float> (m_shift(rnd::reng));
			}
		}

	}

	//overwrite old pop
	pop = pop2;

}

// output reaction norm

int main() {
	// generations
	for (int g = 0; g < gmax; g++) {
		std::cout << "gen = " << g << "\n";

		for (int t = 0; t < tmax; t++) {
			//update environment
			E = A * std::sin((2 * M_PI * (g * tmax + t)) / (tmax * R)) + B * env_dist(rnd::reng);
			//calculate cue
			Cue = std::normal_distribution<float>(P * E, (1 - P) / 3)(rnd::reng);
			/// Is Cue calculated once for the whole population, or per individual?

			//individual update during lifetime
			for (int i = 0; i < pop.size(); ++i) {

				pop[i].update_I_t(Cue);
				pop[i].update_mismatch(E);

			}

		}
		reproduction();
		
	}




	return 0;
}








