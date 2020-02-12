#ifndef INDIVIDUAL_H_INCLUDED
#define INDIVIDUAL_H_INCLUDED

#include <random>
#include "rnd.hpp"
#include "ann2.hpp"

using namespace ann;

// specify ann structure
using Ann = Network<float,
Layer< Neuron<1, activation::identity>, 3>, // for now, 1 input for env cues
Layer< Neuron<3, activation::identity>, 1>  // one output phenotype value
>;

// individuals
struct Individual {
    Individual();

    Ann ann_dev, ann_life;

    float I_baseline, I_realized;
    float mismatch, n;
    void update_I_g(const float C);
    void update_I_t(const float C);
    inline void update_mismatch(const float E) { mismatch += abs(E - I_realized); }
    float calculate_fitness();
    std::vector<float> get_reaction(const std::vector<float> vec_cues);
};



#endif
