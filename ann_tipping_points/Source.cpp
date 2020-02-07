#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include<numeric>
#include <fstream>

#include "individual.h"



std::vector<Individual> reproduction(std::vector<Individual>& pop, float kd, float ka, float tau) {

  //Calculate fitness
  std::vector<float> fitness;
  std::vector<Individual> tmp_pop;


  for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

    fitness.push_back(pop[i].calculate_fitness(kd, ka, tau));

  }

  float mean_fitness = accumulate(fitness.begin(), fitness.end(), 0.f) / static_cast<float>(fitness.size());

  std::cout << mean_fitness << std::endl;

  //Reproduction
  for (int i = 0; i < static_cast<int>(pop.size()); ++i) {
    if (fitness[i] > 0.f) {
      int off = std::poisson_distribution<int>(fitness[i] / mean_fitness)(rnd::reng);
      for (int j = 0; j < off; ++j) {
        tmp_pop.push_back(pop[i]);
      }

    }

  }

  return tmp_pop;
}

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





int main() {



  //Mutation
  const float mrate = 0.001f;
  const float mmean = 0.f;
  const float mshape = 0.05f;

  //Costs
  const float kd = 0.02f;
  const float ka = 0.01f;
  const float tau = 0.25;
  //float q = 2.2;

  //Varied Parameters
  // A and B were varied between 0 and 1 (A+B = 1) to investigate deterministic vs stochastic env., not found to have a great effect 

  const float A = 1.f;   //Deterministic scaling constant, default value 
  const float B = 0.f;   //stochastic scaling constant, default value

  //R between 1 and 100000, P between 0 and 1
  float R = 10.f;  //Environmental variation
  float P = 1.f;   //predictability
  ///
  const auto& env_dist = std::uniform_real_distribution<float>(-1.f, 1.f); // not explicitly stated in botero 2015
  float E;
  float Cue;

  int gmax = 50000;
  int tmax = 5;


  std::vector<Individual> pop(5000); // population size: 5000
  std::vector<Individual> tmp_pop;

  //Initialization
  E = 0.f;
  if (1.f - P == 0.f) {
    Cue = E;
  }
  else {
    Cue = std::normal_distribution<float>(P * E, ((1.f - P) / 3.f))(rnd::reng);
  }

  for (int i = 0; i < static_cast<int>(pop.size()); i++) {
    pop[i].update_I_g(Cue);
  }


  for (int g = 0; g < gmax; g++) {

    std::cout << g <<"\t";


    for (int t = 0; t < tmax; t++) {

      //update environment
      E = A * std::sinf((2.f * static_cast<float>(M_PI)* static_cast<float>(static_cast<float>(g)* tmax + t)) / (static_cast<float>(tmax)* R)) + B * env_dist(rnd::reng);

      //calculate cue
      if (1.f - P == 0.f) {
        Cue = E;
      }
      else {
        Cue = std::normal_distribution<float>(P * E, ((1.f - P) / 3.f))(rnd::reng);
      }
      /// Is Cue calculated once for the whole population, or per individual?

      //individual update
      for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

        pop[i].update_I_t(Cue);         //Insulation
        pop[i].update_mismatch(E);      //phenotypic mismatch

      }

    }

    //Reproduction
    tmp_pop = reproduction(pop, kd, ka, tau);

    //Adjust pop size
    adjust_popsize(tmp_pop, pop.size());


    //Mutation
    for (int i = 0; i < static_cast<int>(tmp_pop.size()); i++) {

      tmp_pop[i].mutate(mrate, mmean, mshape);
      tmp_pop[i].update_I_g(Cue);
    }

    std::swap(pop, tmp_pop);
    tmp_pop.clear();
  }


  std::ofstream ofs("data.txt");

  ofs << "ind" << "\t" << "h" << "\t" << "I01" << "\t" << "I02" << "\t" << "b1" << "\t" << "b2" << "\t" << "s" << "\t" << "a" << "\n";

  for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

    ofs << i << "\t" << pop[i].h << "\t" << pop[i].I01 << "\t" << pop[i].I02 << "\t" << pop[i].b1 << "\t" << pop[i].b2 << "\t" << pop[i].s << "\t" << pop[i].a << "\n";

  }

  ofs.close();



  return 0;
}








