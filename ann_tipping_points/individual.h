#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <random>
#include "rnd.hpp"

struct Individual {

  Individual() : h(1.f), mismatch(0.f) {
  
    const auto& dist1 = std::uniform_real_distribution<float>(0.f, 1.f);
    const auto& dist2 = std::uniform_real_distribution<float>(-1.f, 1.f);
    const auto& dist3 = std::uniform_real_distribution<float>(-2.f, 2.f);

    s = dist1(rnd::reng);
    a = dist1(rnd::reng);

    I01 = dist2(rnd::reng);
    I02 = dist2(rnd::reng);

    b1 = dist3(rnd::reng);
    b2 = dist3(rnd::reng);
    

  }


  void update_I_t(const float C) {
    if (s > 0.5) {
      if (a > 0) {
        if (h_coin) {
          I_realized = b1 * C + I01;
        }
        else {
          I_realized = b2 * C + I02;
        }
      }

    }
  
  }
  void update_I_g(float C) {
  
    h_coin = std::bernoulli_distribution(h)(rnd::reng);

    if (s > 0.5) {
      if (h_coin) {
        I_realized = b1 * C + I01;
      }
      else {
        I_realized = b2 * C + I02;
      }

    }
    else {
      if (h_coin) {
        I_realized = I01;
      }
      else {
        I_realized = I02;
      }

    }
  
  }



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
  float I_realized;
  float mismatch;

};



#endif