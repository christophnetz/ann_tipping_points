#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include<numeric>
#include <fstream>
#include <string>

#include "individual.h"



//Mutation
const float mrate = 0.001f;
const float mmean = 0.f;
const float mshape = 0.05f;

//Costs
const float kd = 0.02f;
const float ka = 0.01f;
const float tau = 0.25;
const float q = 2.2f;

const int gmax = 50000;
const int gext = 1000;
const int tmax = 5;
const int popsize = 5000;

const float A = 1.f;   //Deterministic scaling constant, default value 
const float B = 0.f;   //stochastic scaling constant, default value
const auto& env_dist = std::uniform_real_distribution<float>(-1.f, 1.f); // not explicitly stated in botero 2015


std::vector<Individual> reproduction(std::vector<Individual>& pop) {

  //Calculate fitness
  std::vector<float> fitness;
  std::vector<Individual> tmp_pop;


  for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

    fitness.push_back(pop[i].calculate_fitness(kd, ka, tau));

  }

  float mean_fitness = accumulate(fitness.begin(), fitness.end(), 0.f) / static_cast<float>(fitness.size());

  //std::cout << mean_fitness << std::endl;

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


//for flexible pop size
std::vector<Individual> free_reproduction(std::vector<Individual>& pop) {

  //Calculate fitness
  std::vector<float> fitness;
  std::vector<Individual> tmp_pop;


  for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

    fitness.push_back(pop[i].calculate_fitness(kd, ka, tau));

  }


  //std::cout << mean_fitness << std::endl;

  //Reproduction
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
    int remove = std::uniform_int_distribution<int>(0, static_cast<int>(tmp_pop.size()) - 1)(rnd::reng);
    tmp_pop.erase(tmp_pop.begin() + remove);
  }

  return tmp_pop;
}

void adjust_popsize(std::vector<Individual>& tmp_pop, const int targetsize) {

  while (static_cast<int>(tmp_pop.size()) < targetsize) {
    int duplicate = std::uniform_int_distribution<int>(0, static_cast<int>(tmp_pop.size()) - 1)(rnd::reng);
    tmp_pop.push_back(tmp_pop[duplicate]);
  }

  while (static_cast<int>(tmp_pop.size()) > targetsize) {
    int remove = std::uniform_int_distribution<int>(0, static_cast<int>(tmp_pop.size()) - 1)(rnd::reng);
    tmp_pop.erase(tmp_pop.begin() + remove);
  }

}

void simulation1(std::vector<Individual>& pop, const float& P, const float& R) {

  std::vector<Individual> tmp_pop;
  float E;
  float Cue;
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

    //std::cout << g << "\t";

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
    tmp_pop = reproduction(pop);

    //Adjust pop size
    adjust_popsize(tmp_pop, static_cast<int>(pop.size()));


    //Mutation
    for (int i = 0; i < static_cast<int>(tmp_pop.size()); i++) {

      tmp_pop[i].mutate(mrate, mmean, mshape);
      tmp_pop[i].update_I_g(Cue);
    }

    std::swap(pop, tmp_pop);
    tmp_pop.clear();
  }


}

void simulation2(const float& P, const float& R, const float& Pold, const float& Rold, std::vector<Individual> pop, float& extinction, float& g_extinction) {
  float E = 0.f;
  float Cue = 0.f;
  std::vector<Individual> tmp_pop;

  //Phase maintenance
  int tnew = static_cast<int>(roundf(gmax * tmax * R / Rold));
  int g_init = tnew / tmax;
  int t_init = tnew % tmax;

  g_extinction += static_cast<float>(gext);

  for (int g = g_init; g < g_init + gext; g++) {

    //std::cout << g << "\t";


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
      /// Is Cue calculated once for the whole population, or per individual?

      //individual update
      for (int i = 0; i < static_cast<int>(pop.size()); ++i) {

        pop[i].update_I_t(Cue);         //Insulation
        pop[i].update_mismatch(E);      //phenotypic mismatch

      }

    }

    //Reproduction
    tmp_pop = free_reproduction(pop);

    if (tmp_pop.size() == 0) {
      extinction += 1.f;
      g_extinction += static_cast<float>(g - g_init);
      g_extinction -= static_cast<float>(gext);
      std::cout << "Extinction!!!    " << g_extinction << std::endl;
      break;
    }

    //Mutation
    for (int i = 0; i < static_cast<int>(tmp_pop.size()); i++) {

      tmp_pop[i].mutate(mrate, mmean, mshape);
      tmp_pop[i].update_I_g(Cue);
    }

    std::swap(pop, tmp_pop);
    tmp_pop.clear();
  }
}





int main() {





  //Varied Parameters
  // A and B were varied between 0 and 1 (A+B = 1) to investigate deterministic vs stochastic env., not found to have a great effect 



  //R between 1 and 100000, P between 0 and 1
  std::vector<float> vecR = { 1.f, powf(10.f, 0.5f), 10.f, powf(10.f, 1.5f), 100.f, powf(10.f, 2.5f), 1000.f, powf(10.f, 3.5f), 10000.f, powf(10.f, 4.5f), 100000.f };
  std::vector<float> vecP = { 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f };
  //float R = 10.f;  //Environmental variation
  //float P = 1.f;   //predictability
  ///



  const std::string outfile2 = "extinction_data.csv";
  std::ofstream ofs2(outfile2);
  ofs2 << "R,P,R_new,P_new,extinct,gen_extinct" << "\n";




  for (int r = 0; r < vecR.size(); ++r) {
    float R = vecR[r];
    for (int p = 0; p < vecP.size(); ++p) {
      float P = vecP[p];



      std::vector<Individual> pop(popsize); // population size: 5000

      simulation1(pop, P, R);

      //    std::string filetype = ".png";
      const std::string outfile = "data_logR" + std::to_string(log10f(R)).substr(0, 3) + "_P" + std::to_string(P).substr(0, 3) + ".txt";
      std::ofstream ofs(outfile);
      ofs << "ind" << "\t" << "h" << "\t" << "I01" << "\t" << "I02" << "\t" << "b1" << "\t" << "b2" << "\t" << "s" << "\t" << "a" << "\n";

      for (int i = 0; i < static_cast<int>(pop.size()); ++i) {
        ofs << i << "\t" << pop[i].h << "\t" << pop[i].I01 << "\t" << pop[i].I02 << "\t" << pop[i].b1 << "\t" << pop[i].b2 << "\t" << pop[i].s << "\t" << pop[i].a << "\n";
      }
      ofs.close();


      ///////////////////////////////////////////////////////
      ///Transition
      //Questions: shift population just once/ a hundred times?
      //What's this relative extinction rate?

      

      for (int r_new = 0; r_new < vecR.size(); ++r_new) {
        float R_new = vecR[r_new];
        for (int p_new = 0; p_new < vecP.size(); ++p_new) {
          float P_new = vecP[p_new];
          float extinction = 0.f;
          float g_extinction = 0.f;
          for (int sim2 = 0; sim2 < 10; sim2++) {
            simulation2(P_new, R_new, P, R, pop, extinction, g_extinction);
          }
          ofs2 << log10f(R) << "," << P << "," << log10f(R_new) << "," << P_new << "," << extinction / 10.f << "," << g_extinction / 10.f<< "\n";
          std::cout << g_extinction << std::endl;
        }

      }
      



    }

  }

  ofs2.close();


  return 0;
}








