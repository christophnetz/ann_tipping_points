#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <fstream>
#include <string>
#include "individual.h"

const int popsize = 5000;
//Varied Parameters
// A and B were varied between 0 and 1 (A+B = 1) to investigate deterministic vs stochastic env., not found to have a great effect 

float A = 1.f;   //Deterministic scaling constant, default value 
float B = 0.f;   //stochastic scaling constant, default value

//R between 1 and 100000, P between 0 and 1
float R = 10.f;  //Environmental variation
float P = 0.99f;   //predictability
///
const auto& env_dist = std::uniform_real_distribution<float>(-1.f, 1.f); // not explicitly stated in botero 2015

int gmax = 50000;
const int gext = 1000;
int tmax = 5;

// magic number of offsprin
float q = 2.2f;

//Mutation
const float mrate = 0.001f;
const float mmean = 0.f;
const float mshape = 0.05f;

//Costs
const float kd = 0.02f;
const float ka = 0.01f;
const float tau = 0.25;
//float q = 2.2;

//Cue range for reaction norm
const float cue_max = 1.f;
const float cue_min = 0.f;
const float cue_inc = 0.05f;

//R between 1 and 100000, P between 0 and 1
std::vector<float> vecR = { 1.f, powf(10.f, 0.5f), 10.f, powf(10.f, 1.5f), 100.f, powf(10.f, 2.5f), 1000.f, powf(10.f, 3.5f), 10000.f, powf(10.f, 4.5f), 100000.f };
std::vector<float> vecP = { 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f };

// standard vector of cues
std::vector<float> vec_cues = { 0.f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.f };

// adjust pop size
void adjust_popsize(std::vector<Individual>& tmp_pop, const int targetsize) {

	while (static_cast<int>(tmp_pop.size()) < targetsize) {
		int duplicate = std::uniform_int_distribution<int>(0, tmp_pop.size() - 1)(rnd::reng);
		tmp_pop.push_back(tmp_pop[duplicate]);
	}

	while (static_cast<int>(tmp_pop.size()) > targetsize) {
		int remove = std::uniform_int_distribution<int>(0, tmp_pop.size() - 1)(rnd::reng);
		tmp_pop.erase(tmp_pop.begin() + remove);
	}

}


/// reproduction with fixed pop size
std::vector<Individual> reproduction(std::vector<Individual>& pop) {

  //Calculate fitness
  std::vector<float> fitness;
  std::vector<Individual> tmp_pop;


  for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

    fitness.push_back(pop[i].calculate_fitness_ann(ka, tau));

  }

  float mean_fitness = accumulate(fitness.begin(), fitness.end(), 0.f) / static_cast<float>(fitness.size());

  //Reproduction
  for (int i = 0; i < static_cast<int>(pop.size()); ++i) 
  {
    if (fitness[i] > 0.f) 
    {
      int off = std::poisson_distribution<int>(fitness[i] / mean_fitness)(rnd::reng);
      for (int j = 0; j < off; ++j) {
        tmp_pop.push_back(pop[i]);
      }
    }
  }

  // must have mutation on the temporary population
  for (int i = 0; i < static_cast<int>(tmp_pop.size()); ++i)
  {
  	tmp_pop[i].mutate();
  }

  // adjust pop size
  adjust_popsize(tmp_pop, pop.size());

  return tmp_pop;
}

/// function for free reproduction
std::vector<Individual> free_reproduction(std::vector<Individual>& pop) {

	//Calculate fitness
	std::vector<float> fitness;
	std::vector<Individual> tmp_pop;


	for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

		fitness.push_back(pop[i].calculate_fitness_ann(ka, tau));

	}

	//Reproduction according to botero
	for (int i = 0; i < static_cast<int>(pop.size()); ++i) {
		if (fitness[i] > 0.f) {
			int off = std::poisson_distribution<int>(q * fitness[i] / 1.0f)(rnd::reng); //Wmax is equal to 1, right?
			for (int j = 0; j < off; ++j) {
				tmp_pop.push_back(pop[i]);
			}

		}

	}

	//incorporated pop adjust
	while (static_cast<int>(tmp_pop.size()) > popsize) {
		int remove = std::uniform_int_distribution<int>(0, tmp_pop.size() - 1)(rnd::reng);
		tmp_pop.erase(tmp_pop.begin() + remove);
	}

	return tmp_pop;
}


/// plot output
void print_reaction_norm(const float R, const float P, std::vector<Individual>& pop,
	const std::vector<float> vec_cues)
{
	// get filename
	const std::string outfile = "data/data_ann_logR" + std::to_string(log10f(R)).substr(0, 3) + "_P" + std::to_string(P).substr(0, 3) + ".txt";
	std::ofstream ofs(outfile);

	ofs << "ind,baseline,cue,resp" << "\n";

	for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

		std::vector<float> vec_resp = pop[i].get_reaction(vec_cues);
		float base = pop[i].I_baseline;

		for (int j = 0; j < static_cast<int>(vec_cues.size()); j++) {

			ofs << i << "," << base << "," << vec_cues[j] << "," << vec_resp[j] << "\n";
		}
	}
	ofs.close();
}

/// print extinction data
void print_extinction_data(const float R, const float P, const float R_new, 
	const float P_new, const int which_gen)
{
	// filename for ofstream
	const std::string extinct_out = "data/extinction_data.csv";
	std::ofstream ext_ofs;

	// check if DOES NOT exist then write col names
	std::ifstream f(extinct_out.c_str());
	if (!f.good()) {
		ext_ofs.open(extinct_out, std::ofstream::out);
		ext_ofs << "R,P,R_new,P_new,gen_extinct" << "\n";
	}
	
	// then append data
	ext_ofs.open(extinct_out, std::ofstream::app);
	ext_ofs << log10f(R) << "," << P << "," 
		<< log10f(R_new) << "," << P_new << "," 
		<< which_gen << "\n";
	
}

/// function to evolve population
std::vector<Individual> evolve_pop(std::vector<Individual> pop, const float R, const float P)
{
	//Initialization
	std::vector<Individual> tmp_pop;
	float E = 0.f;
	float Cue;
	if (1.f - P == 0.f) {
		Cue = E;
	}
	else {
		Cue = std::normal_distribution<float>(P * E, ((1.f - P) / 3.f))(rnd::reng);
	}

	for (int i = 0; i < static_cast<int>(pop.size()); i++) {
		pop[i].update_I_g(Cue);
	}

	// print info
	std::cout << "R = " << R << " P = " << P << "\n";
	// generations
	for (int g = 0; g < gmax; g++)
	{
		for (int t = 0; t < tmax; t++) {
			//update environment
			E = A * std::sinf((2 * static_cast<float>(M_PI)* (static_cast<float>(g)* static_cast<float>(tmax) +
				static_cast<float>(t))) / (static_cast<float>(tmax)* R)) + B * env_dist(rnd::reng);
			//calculate cue
			if (1.f - P == 0.f) {
				Cue = E;
			}
			else {
				Cue = std::normal_distribution<float>(P * E, ((1.f - P) / 3.f))(rnd::reng);
			}

			//individual update during lifetime
			for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

				pop[i].update_I_t(Cue);
				pop[i].update_mismatch(E);

			}

		}
		pop = reproduction(pop);
	}
	std::cout << "pop evolved...\n";
	return pop;
}

/// function to shift population and assess extinction
void test_extinction(std::vector<Individual> pop, const float R, const float P, const float R_new, const float P_new)
{
	float E; float Cue;
	std::vector<Individual> tmp_pop;

	// implement phase maintenance
	int tnew = static_cast<int>(roundf(gmax * tmax * R_new / R));
	int g_init = tnew / tmax;
	int t_init = tnew % tmax;

	// continue evolution with new R and P
	{
		for (int g = g_init; g < g_init + gext; g++) {

			std::cout << "new R = " << log10f(R_new) << " new P = " << P_new << "\n";

			if ((g - g_init) % 100 == 0 || (g - g_init) == 0) {
				std::cout << "g = " << g - g_init << " popsize = " << pop.size() << "\n";
			}

			for (int t = t_init; t < t_init + tmax; t++) {
				//the generations are now shifted by t_init, to keep the phase displacement to a minimum.
				//All but pretty, but i don't see any other way and i think it works correctly

				//update environment

				E = std::sinf((2.f * static_cast<float>(M_PI)* static_cast<float>(static_cast<float>(g)* tmax + t)) / (static_cast<float>(tmax)* R));

				//calculate cue
				if (1.f - P == 0.f) {
					Cue = E;
				}
				else {
					Cue = std::normal_distribution<float>(P * E, ((1.f - P) / 3.f))(rnd::reng);
				}
				
				//individual update
				for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

					pop[i].update_I_t(Cue);         //Insulation
					pop[i].update_mismatch(E);      //phenotypic mismatch

				}

			}

			//Reproduction with flexible pop size
			tmp_pop = free_reproduction(pop);

			if (tmp_pop.size() == 0) {
				// print extinction data and break
				print_extinction_data(R, P, R_new, P_new, g - g_init);
				std::cout << "pop extinct! g = " << g - g_init << "\n\n";
				break;
			}

			//Mutation
			for (int i = 0; i < static_cast<int>(tmp_pop.size()); i++) {

				tmp_pop[i].mutate();
				tmp_pop[i].update_I_g(Cue);
			}

			std::swap(pop, tmp_pop);
			tmp_pop.clear();

			if (g == g_init + gext)
			{
				std::cout << "pop survived!\n";
				print_extinction_data(R, P, R_new, P_new, g);
			}
		}

	}
}

/// main function
int main() {

	for (int r = 0; r < static_cast<int>(vecR.size()); ++r) {
		float R = vecR[r];
		for (int p = 0; p < static_cast<int>(vecP.size()); ++p) {
			float P = vecP[p];

			// make agents and evolve and print evolved reaction norm
			std::vector<Individual> pop(popsize);
			pop = evolve_pop(pop, R, P);
			print_reaction_norm(R, P, pop, vec_cues);

			// test evolved agents for extinction
			for (int r_new = 0; r_new < static_cast<int>(vecR.size()); ++r_new) {
				float R_new = vecR[r_new];
				for (int p_new = 0; p_new < static_cast<int>(vecP.size()); ++p_new) {
					float P_new = vecP[p_new];
					test_extinction(pop, R, P, R_new, P_new);
				}
			}

		}
	}

	return 0;
}








