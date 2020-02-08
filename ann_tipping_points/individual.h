#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <random>
#include "rnd.hpp"
#include "ann.h"

using namespace ann;

// specify ann structure
using Ann = Network<float,
Layer< Neuron<1, activation::rtlu>, 3>, // for now, 1 input for env cues
Layer< Neuron<3, activation::rtlu>, 3>,
Layer< Neuron<3, activation::rtlu>, 1>  // one output phenotype value
>;

// individuals
struct Individual {

    Individual() : mismatch(0.f), ann_dev(0.f), ann_life(0.f), n(0.f) {};

    Ann ann_dev, ann_life;

    float I_baseline, I_realized;
    float mismatch, n;
    void update_I_g(const float C);
    void update_I_t(const float C);
    inline void update_mismatch(const float E) { mismatch += abs(E - I_realized); }
    float calculate_fitness(float kd, float ka, float tau);
    std::vector<float> get_reaction(const std::vector<float> vec_cues);
};



#endif
