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

	Individual() : mismatch(0.f), ann_dev(0.f), ann_life(0.f) {};

	Ann ann_dev, ann_life;

	float I_baseline, I_realized;
	float mismatch;

	// response to development cue
	void update_I_g(const float C) {
		Ann::input_t inputs_g;
		inputs_g[0] = C;
		auto output_g = ann_dev(inputs_g);
		I_baseline = output_g[0];
	}

	// lifetime modification of baseline value
	void update_I_t(const float C) {
		Ann::input_t inputs_t;
		inputs_t[0] = C;
		auto output_t = ann_dev(inputs_t);
		I_realized = I_baseline + output_t[0];
	}
};



#endif