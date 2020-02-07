#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>

#include "individual.h"






int main() {

  std::vector<Individual> pop(5000); // population size: 5000
  //Varied Parameters
  // A and B were varied between 0 and 1 (A+B = 1) to investigate deterministic vs stochastic env., not found to have a great effect 
  
  float A = 1.f;   //Deterministic scaling constant, default value 
  float B = 0.f;   //stochastic scaling constant, default value

  //R between 1 and 100000, P between 0 and 1
  float R = 10.f;  //Environmental variation
  float P = 1.f;   //predictability
  ///
  const auto& env_dist = std::uniform_real_distribution<float>(-1.f, 1.f); // not explicitly stated in botero 2015
  float E;
  float Cue;

  int gmax = 50000;
  int tmax = 5;

  std::cout << "Hello world" << std::endl;

  for (int g = 0; g < gmax; g++) {

    for (int t = 0; t < tmax; t++) {
      //update environment
      E = A * std::sin( (2 * M_PI * (g * tmax + t))/(tmax * R)) + B * env_dist(rnd::reng);
      //calculate cue
      Cue = std::normal_distribution<float>(P * E, (1 - P) / 3)(rnd::reng);
      /// Is Cue calculated once for the whole population, or per individual?

      //individual update
      for (int i = 0; i < pop.size(); ++i) {

        pop[i].update_I_t(Cue);

      }


      //phenotypic mismatch


    }

    //Reproduction



  }




  return 0;
}








