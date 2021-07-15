#include "individual.h"




Individual::Individual() {

  h = 1.f;
  mismatch = 0.f;
  n = 0;
  h_coin = true;
  I_realized = 0.f;

  auto dist1 = std::uniform_real_distribution<float>(0.f, 1.f);
  auto dist2 = std::uniform_real_distribution<float>(-1.f, 1.f);
  auto dist3 = std::uniform_real_distribution<float>(-2.f, 2.f);

  s = dist1(rnd::reng);

  a = dist1(rnd::reng);

  b1 = dist3(rnd::reng);
  b2 = dist3(rnd::reng);



  I01 = dist2(rnd::reng);
  I02 = dist2(rnd::reng);

}


void Individual::update_I_t(const float C) {
  if (s > 0.5f) {
    if (a > 0.f) {

      if (a >= 1 || std::bernoulli_distribution(a)(rnd::reng)) {

        n++;

        if (h_coin) {
          I_realized = b1 * C + I01;
        }
        else {
          I_realized = b2 * C + I02;
        }
      }
    }

  }
}


void Individual::update_I_g(float C) {

  if (h >= 1) {
    h_coin = true;
  }
  else {
    h_coin = std::bernoulli_distribution(h)(rnd::reng);
  }

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


float Individual::calculate_fitness(float kd, float ka, float tau) {

  float fit;
  if (s > 0.5) {
    fit = exp(-tau * mismatch) - kd - n * ka;
  }
  else {
    fit = exp(-tau * mismatch);
  }
  mismatch = 0;
  n = 0;

  if (fit < 0.f)
    fit = 0.f;

  return fit;

}


void Individual::mutate(float mrate, float mmean, float mshape) {

  std::bernoulli_distribution mdist(mrate);
  std::normal_distribution<float> sdist(mmean, mshape);

  if (mdist(rnd::reng)) {
    h += sdist(rnd::reng);
  }
  if (mdist(rnd::reng)) {
    I01 += sdist(rnd::reng);
  }
  if (mdist(rnd::reng)) {
    I02 += sdist(rnd::reng);
  }
  if (mdist(rnd::reng)) {
    s += sdist(rnd::reng);
  }
  //Unclear: mutation of s before or after check?

  if (mdist(rnd::reng)) {
    b1 += sdist(rnd::reng);
  }
  if (mdist(rnd::reng)) {
    b2 += sdist(rnd::reng);
  }
  if (mdist(rnd::reng)) {
    a += sdist(rnd::reng);
  }



}
