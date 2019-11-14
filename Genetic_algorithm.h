//
// Created by cheshirecat on 11/14/19.
//

#ifndef GENETIC_ALGORITHM_GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_GENETIC_ALGORITHM_H

#include <vector>

class Individual
{
public:
    std::vector<double> chromosome;
    double fitness;
    double chance; // is chance to be selected at crossover
};

class Genetic_algorithm {
public:
    Genetic_algorithm(double gMax, double gMin, double (*fitnessCalc)(std::vector<double>), unsigned int crossoverP,
                      double classicMutChan, double sigmaMutChan, unsigned int populationSize);

private:
    //g_max >= g_min
    double g_max; // gene upper bound
    double g_min; // gene lower bound
    double (*fitness_calc) (std::vector<double>);
    unsigned int crossover_p; //amount of cut points at crossover. >= 1;
    //classic_mut_chan + sigma_mut_chan <= 1
    double classic_mut_chan; //chance of classic mutation
    double sigma_mut_chan; //chance of sigma mutation;
    unsigned int population_size;
    std::vector<Individual> population;
    std::vector<Individual*> child; //child at last iteration
};


#endif //GENETIC_ALGORITHM_GENETIC_ALGORITHM_H