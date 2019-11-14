//
// Created by cheshirecat on 11/14/19.
//

#ifndef GENETIC_ALGORITHM_GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_GENETIC_ALGORITHM_H

#include <vector>
#include <random>

class Individual
{
public:
    std::vector<double> chromosome;
    double fitness;
    double chance; // is chance to be selected at crossover

    explicit Individual(unsigned int chrom_size) : fitness{0}, chance{0}, chromosome(chrom_size) { }
};

class Genetic_algorithm {
public:
    Genetic_algorithm(double gMin, double gMax, double (*fitnessCalc)(std::vector<double>&), unsigned int crossoverP,
                      double classicMutChan, double sigmaMutChan, unsigned int populationSize,
                      unsigned int chromosomeSize);

    [[nodiscard]] const std::vector<Individual> &getPopulation() const;
    [[nodiscard]] const std::vector<Individual *> &getChild() const;


private:
    //g_max >= g_min
    double g_max; // gene upper bound
    double g_min; // gene lower bound
    double (*fitness_calc) (std::vector<double>&);
    unsigned int crossover_p; //amount of cut points at crossover. >= 1;
    //classic_mut_chan + sigma_mut_chan <= 1
    double classic_mut_chan; //chance of classic mutation
    double sigma_mut_chan; //chance of sigma mutation;
    unsigned int population_size;
    std::vector<Individual> population;
    std::vector<Individual*> child; //child at last iteration
    unsigned int chromosome_size;

    std::random_device rd;
    std::default_random_engine random_engine;


    void initialization();
    void suitability_rating();
};


#endif //GENETIC_ALGORITHM_GENETIC_ALGORITHM_H
