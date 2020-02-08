#include "individual.h"

// if ann_life weights are different from 0, there is a cost

float Individual::calculate_fitness(float kd, float ka, float tau) {

	float tot_mut_dist;
	for (auto& w : pop2[a].ann_life) {
		tot_mut_dist += abs(w);
	}
  
  float fit;
  if (tot_mut_dist > 0.f) {
    fit = exp(-tau * mismatch) - kd - n * ka;
  }
  else {
    fit = exp(-tau * mismatch);
  }
  mismatch = 0.f;

  if (fit < 0.f)
    fit = 0.f;

  return fit;

}


