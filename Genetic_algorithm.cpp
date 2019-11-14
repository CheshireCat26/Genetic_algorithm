//
// Created by cheshirecat on 11/14/19.
//

#include "Genetic_algorithm.h"
#include <stdexcept>
#include <ctime>

using namespace std;

Genetic_algorithm::Genetic_algorithm(double gMax, double gMin, double (*fitnessCalc)(vector<double>&),
                                     unsigned int crossoverP, double classicMutChan, double sigmaMutChan,
                                     unsigned int populationSize, unsigned int chromosomeSize)
                                     : g_max(gMax), g_min(gMin), fitness_calc(fitnessCalc), crossover_p(crossoverP),
                                       classic_mut_chan(classicMutChan), sigma_mut_chan(sigmaMutChan),
                                       population_size(populationSize), chromosome_size(chromosomeSize)
{
    random_engine.seed(rd());
    if (gMax < gMax)
        throw runtime_error("Genetic_algorithm::Genetic_algorithm(): gMax must be <= gMin");
    if (classicMutChan + sigmaMutChan > 1)
        throw runtime_error("Genetic_algorithm::Genetic_algorithm(): classicMutChan + sigmaMutChan must be <= 1");

    initialization();
    suitability_rating();
}

void Genetic_algorithm::initialization() {
    uniform_real_distribution distribution(g_min, g_max);
    for (unsigned int i = 0; i < population_size; i++)
    {
        Individual individual(chromosome_size);
        for (unsigned int j = 0; j < chromosome_size; j++)
        {
            individual.chromosome[j] = distribution(random_engine);
        }
        population.emplace_back(individual);
    }
}

const vector<Individual> &Genetic_algorithm::getPopulation() const {
    return population;
}

void Genetic_algorithm::suitability_rating() {
    for (Individual& individual : population)
        individual.fitness = fitness_calc(individual.chromosome);
}

const vector<Individual *> &Genetic_algorithm::getChild() const {
    return child;
}
