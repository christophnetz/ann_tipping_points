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
    for (int i = 0; i < static_cast<int>(pop.size()); ++i)
    {
        if (fitness[i] > 0.f)
        {
            int off = std::poisson_distribution<int>(fitness[i] / mean_fitness)(rnd::reng);

            for (int j = 0; j < off; ++j)
            {
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

    if(tmp_pop.size() > popsize)
    {
        //incorporated pop adjust
        std::shuffle(tmp_pop.begin(), tmp_pop.end(), rnd::reng);
        tmp_pop.resize(popsize);
    }
    return tmp_pop;
}

void adjust_popsize(std::vector<Individual>& tmp_pop, const size_t targetsize) {

    while (tmp_pop.size() < targetsize) {
        int duplicate = std::uniform_int_distribution<int>(0, tmp_pop.size() - 1)(rnd::reng);
        tmp_pop.push_back(tmp_pop[duplicate]);
    }

    if(tmp_pop.size() > targetsize)
    {
        //incorporated pop adjust
        std::shuffle(tmp_pop.begin(), tmp_pop.end(), rnd::reng);
        tmp_pop.resize(targetsize);
    }

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
        for (int t = 0; t < tmax; t++)
        {
            //update environment
            E = A * std::sin((2.f * M_PI * (g * tmax + t)) / (tmax* R)) + B * env_dist(rnd::reng);

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

        mut_dist.mutate_transform(pop.begin(),
                                  pop.end(),
                                  [](Individual& i) {return i.calculate_fitness(kd, ka, tau);});

        for( size_t i = 0; i != pop.size(); i++)
        {
            tmp_pop[i] = pop[mut_dist(rnd::reng)];
        }
        //    //Reproduction
        //    pop = reproduction(pop);

        //    //Adjust pop size
        //    adjust_popsize(tmp_pop, popsize);


        //Mutation
        for (int i = 0; i < static_cast<int>(tmp_pop.size()); i++)
        {
            tmp_pop[i].mutate(mrate, mmean, mshape);
            tmp_pop[i].update_I_g(Cue);
        }

        std::swap(pop, tmp_pop);
        if(g % 1000 == 0)
        {
            std::cout<< "gen: " << g << std::endl;
        }
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





int main(int argc, char ** argv) {


    cxxopts::Options options("",
                             "Insert the parameters for the simualtion and see if you can get a mutational switch to evolve");
    options.add_options()
            ("P,predictability", "the predictability of the environment", cxxopts::value<float>())
            ("R,repeatability", "the repeatability of the environment", cxxopts::value<float>())
            ;

    auto results = options.parse(argc,argv);

    //Varied Parameters
    // A and B were varied between 0 and 1 (A+B = 1) to investigate deterministic vs stochastic env., not found to have a great effect



    //R between 1 and 100000, P between 0 and 1
    std::vector<float> vecR = { 1.f,
                                3.16228f,
                                10.f,
                                31.6228f,
                                100.f,
                                316.228f,
                                1000.f,
                                3162.28f,
                                10000.f,
                                31622.8f,
                                100000.f
                              };

    std::vector<float> vecP = { 0.0f,
                                0.1f,
                                0.2f,
                                0.3f,
                                0.4f,
                                0.5f,
                                0.6f,
                                0.7f,
                                0.8f,
                                0.9f,
                                1.0f
                              };

    float original_R = results["R"].as<float>();
    float P = results["P"].as<float>();

    std::vector<Individual> pop(popsize); // population size: 5000

    simulation1(pop, P, original_R);

    const std::string outfile = "data_logR" +
            std::to_string(log10f(original_R)).substr(0, 3) +
            "_P" +
            std::to_string(P).substr(0, 3) +
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

    for (size_t r_new = 0; r_new < vecR.size(); ++r_new)
    {
        float R_new = vecR[r_new];
        for (size_t p_new = 0; p_new < vecP.size(); ++p_new)
        {
            float P_new = vecP[p_new];
            float extinction = 0.f;
            float g_extinction = 0.f;
            for (int sim2 = 0; sim2 < 10; sim2++)
            {
                simulation2(P_new, R_new, original_R, pop, extinction, g_extinction);
            }
            ofs2 << log10f(original_R) << "," <<
                    P << "," << log10f(R_new) << "," <<
                    P_new << "," << extinction / 10.f << "," <<
                    g_extinction / 10.f<< "\n";
            //std::cout << g_extinction << std::endl;
        }
    }
    ofs2.close();

    return 0;
}
