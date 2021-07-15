#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <random>
#include "rnd.hpp"

struct Individual {

  Individual();


  void update_I_t(const float C);

  void update_I_g(float C); 

  inline void update_mismatch(const float E) { mismatch += abs(E - I_realized); }

  float calculate_fitness(float kd, float ka, float tau);


  void mutate(float mrate, float mmean, float mshape);


  //genes
  float h;
  float I01;
  float I02;
  float b1;
  float b2;
  float s;
  float a;

  //lifetime params
  bool h_coin;
  int n;
  float I_realized;
  float mismatch;

};



#endif
