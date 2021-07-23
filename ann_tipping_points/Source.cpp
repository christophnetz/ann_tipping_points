#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include<numeric>
#include <fstream>
#include <string>
//Version 2.1 of https://github.com/jarro2783/cxxopts
#include "cxxopts.hpp"

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
auto env_dist = std::uniform_real_distribution<float>(-1.f, 1.f); // not explicitly stated in botero 2015

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

    if (tmp_pop.size() > popsize)
    {
        //incorporated pop adjust
        std::shuffle(tmp_pop.begin(), tmp_pop.end(), rnd::reng);
        tmp_pop.resize(popsize);
    }
    return tmp_pop;
}

void simulation1(std::vector<Individual>& pop, const float& P, const float& R) {

    std::vector<Individual> tmp_pop(popsize);
    rndutils::mutable_discrete_distribution<> mut_dist;

    float E = 0.0;
    float Cue;
    //Initialization
    if (1.f - P == 0.f)
    {
        Cue = E;
    }
    else
    {
        Cue = std::normal_distribution<float>(P * E, ((1.f - P) / 3.f))(rnd::reng);
    }

    for (size_t i = 0; i < pop.size(); i++)
    {
        pop[i].update_I_g(Cue);
    }


    for (int g = 0; g < gmax; g++)
    {
        if (g % 100 == 0) {
            std::cout << "Gen" << g << std::endl;
        }
        for (int t = 0; t < tmax; t++)
        {
            //update environment
            E = A * std::sin((2.f * M_PI * (g * tmax + t)) / (tmax * R)) + B * env_dist(rnd::reng);

            //calculate cue
            if (1.f - P == 0.f)
            {
                Cue = E;
            }
            else
            {
                Cue = std::normal_distribution<float>(P * E, ((1.f - P) / 3.f))(rnd::reng);
            }

            //individual update
            for (size_t i = 0; i < pop.size(); ++i) {

                pop[i].update_I_t(Cue);         //Insulation
                pop[i].update_mismatch(E);      //phenotypic mismatch

            }

        }

        //Reproduction
        mut_dist.mutate_transform(pop.begin(),
                                  pop.end(),
                                  [](Individual& i) {return i.calculate_fitness(kd, ka, tau); });

        for (size_t i = 0; i != pop.size(); i++)
        {
            tmp_pop[i] = pop[mut_dist(rnd::reng)];
        }



        //Mutation
        for (int i = 0; i < static_cast<int>(tmp_pop.size()); i++)
        {
            tmp_pop[i].mutate(mrate, mmean, mshape);
            tmp_pop[i].update_I_g(Cue);
        }

        std::swap(pop, tmp_pop);
    }


}

void simulation2(const float& P_new,
                 const float& R_new,
                 const float& Rold,
                 std::vector<Individual> pop,
                 float& extinction,
                 float& g_extinction) {
    float E = 0.f;
    float Cue = 0.f;
    std::vector<Individual> tmp_pop;

    //Phase maintenance
    int tnew = static_cast<int>(roundf(gmax * tmax * R_new / Rold));
    int g_init = tnew / tmax;
    int t_init = tnew % tmax;

    g_extinction += gext;

    for (int g = g_init; g < g_init + gext; g++)
    {
        for (int t = t_init; t < t_init + tmax; t++)
        {
            //the generations are now shifted by t_init, to keep the phase displacement to a minimum.
            //All but pretty, but i don't see any other way and i think it works correctly

            //update environment

            E = std::sin((2.f * M_PI * (g * tmax + t)) / (tmax * R_new));

            //calculate cue
            if (1.f - P_new == 0.f)
            {
                Cue = E;
            }
            else
            {
                Cue = std::normal_distribution<float>(P_new * E, ((1.f - P_new) / 3.f))(rnd::reng);
            }
            /// Is Cue calculated once for the whole population, or per individual?

            //individual update
            for (size_t i = 0; i < pop.size(); ++i)
            {
                pop[i].update_I_t(Cue);         //Insulation
                pop[i].update_mismatch(E);      //phenotypic mismatch
            }
        }

        //Reproduction
        tmp_pop = free_reproduction(pop);

        if (tmp_pop.size() == 0) {
            extinction += 1.f;
            g_extinction += g - g_init;
            g_extinction -= gext;
            //std::cout << "Extinction!!!    " << g_extinction << std::endl;
            break;
        }

        //Mutation
        for (int i = 0; i < static_cast<int>(tmp_pop.size()); i++)
        {
            tmp_pop[i].mutate(mrate, mmean, mshape);
            tmp_pop[i].update_I_g(Cue);
        }

        pop = std::move(tmp_pop);
    }
}





int main(int argc, char** argv) {


    cxxopts::Options options("",
                             "Insert the parameters for the simualtion and see if you can get a mutational switch to evolve");
    options.add_options()
            ("P,predictability", "the predictability of the environment", cxxopts::value<float>())
            ("R,repeatability", "the repeatability of the environment", cxxopts::value<float>())
            ("s, step_p", "the step at whihc to gnerate the vectpr of P values", cxxopts::value<float>())
            ("S, step_r", "the step at whihc to gnerate the vectpr of R values", cxxopts::value<float>())
            ("m, max_p", "the max value of P at which to test extinction", cxxopts::value<float>())
            ("M, max_r", "the max value of R at which to test extinction", cxxopts::value<float>())
            ;

    auto results = options.parse(argc, argv);

    //Varied Parameters
    // A and B were varied between 0 and 1 (A+B = 1) to investigate deterministic vs stochastic env., not found to have a great effect


    //Create 2 very big arrays of R and P so you can always find the value of R and P given by the command line
    int size_par_vecs = 100;
    std::vector<float> vecR(size_par_vecs);
    float start_R = 0.0;
    float step_R = results["step_r"].as<float>();
    std::generate(vecR.begin(),vecR.end(),[&](){return start_R = start_R + step_R;});
    std::for_each(vecR.begin(),vecR.end(),[](float i){return std::pow(10,i);});

    std::vector<float> vecP(size_par_vecs);
    float start_P = 0.0;
    float step_P = results["step_p"].as<float>();
    std::generate(vecP.begin(),vecP.end(),[&](){return start_P = start_P + step_P;});

    float original_R = results["R"].as<float>();
    float P = results["P"].as<float>();

    std::vector<Individual> pop(popsize); // population size: 5000

    simulation1(pop, P, original_R);

    const std::string outfile = "data_logR" +
            std::to_string(log10f(original_R)).substr(0, 6) +
            "_P" +
            std::to_string(P).substr(0, 6) +
            ".txt";

    std::ofstream ofs(outfile);
    ofs << "ind" << "\t" <<
           "h" << "\t" <<
           "I01" << "\t" <<
           "I02" << "\t" <<
           "b1" << "\t" <<
           "b2" << "\t" <<
           "s" << "\t" <<
           "a" << "\n";

    for (int i = 0; i < static_cast<int>(pop.size()); ++i)
    {
        ofs << i << "\t" <<
               pop[i].h << "\t" <<
               pop[i].I01 << "\t" <<
               pop[i].I02 << "\t" <<
               pop[i].b1 << "\t" <<
               pop[i].b2 << "\t" <<
               pop[i].s << "\t" <<
               pop[i].a << "\n";
    }
    ofs.close();

    ///////////////////////////////////////////////////////
    ///Transition
    //Questions: shift population just once/ a hundred times?
    //What's this relative extinction rate?


    const std::string outfile2 = "extinction_data_R" +
            std::to_string(log10f(original_R)).substr(0, 3) +
            "_P" +
            std::to_string(P).substr(0, 3) +
            ".csv";

    std::ofstream ofs2(outfile2);
    ofs2 << "R,P,R_new,P_new,extinct,gen_extinct" << "\n";

    float max_R = std::pow(10, results["max_r"].as<float>()) + 0.0001;
    float max_P = results["max_p"].as<float>() + 0.0001;


    for (int r_new = 0; vecR[r_new] < max_R; ++r_new)
    {
        float R_new = vecR[r_new];
        std::cout << "###iterating extinction R " << R_new << std::endl;

        for (int p_new = 0; vecP[p_new] < max_P; ++p_new)
        {
            float P_new = vecP[p_new];
            std::cout << "iterating extinction P " << P_new << std::endl;
            float extinction = 0.f;
            float g_extinction = 0.f;

            for (int sim2 = 0; sim2 < 10; sim2++)
            {
                simulation2(P_new, R_new, original_R, pop, extinction, g_extinction);
            }
            {
                ofs2 << log10f(original_R) << "," <<
                        P << "," << log10f(R_new) << "," <<
                        P_new << "," << extinction / 10.f << "," <<
                        g_extinction / 10.f << "\n";
            }
        }

    }
    ofs2.close();

    return 0;
}
