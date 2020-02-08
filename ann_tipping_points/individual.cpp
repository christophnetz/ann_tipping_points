#include "individual.h"

/// function development cue
void Individual::update_I_g(const float C) {
    Ann::input_t inputs_g;
    inputs_g[0] = C;
    auto output_g = ann_dev(inputs_g);
    I_baseline = output_g[0];
}

/// function lifetime cue
// lifetime modification of baseline value
void Individual::update_I_t(const float C) {
    Ann::input_t inputs_t;
    inputs_t[0] = C;
    auto output_t = ann_dev(inputs_t);
    I_realized = I_baseline + output_t[0];
}

/// function ann fitness 
float Individual::calculate_fitness(float kd, float ka, float tau) {

    float tot_mut_dist = 0.f;
    for (auto& w : ann_life) {
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

/// function reaction norm
std::vector<float> Individual::get_reaction(const std::vector<float> vec_cues){
    std::vector<float> vec_reaction;

    for(size_t i_cue = 0; i_cue < static_cast<size_t>(vec_cues.size()); i_cue++){
        // get ann output
        Ann::input_t inputs_t;
        inputs_t[0] = vec_cues[i_cue];
        auto output_t = ann_dev(inputs_t);
        I_realized = I_baseline + output_t[0];

        vec_reaction.push_back(I_realized);
    }
    return vec_reaction;
}
