//
// Created by cheshirecat on 11/14/19.
//

#include "Genetic_algorithm.h"
#include <stdexcept>

using namespace std;

Genetic_algorithm::Genetic_algorithm(double gMax, double gMin, double (*fitnessCalc)(std::vector<double>),
                                     unsigned int crossoverP, double classicMutChan, double sigmaMutChan,
                                     unsigned int populationSize) : g_max(gMax), g_min(gMin), fitness_calc(fitnessCalc),
                                                                    crossover_p(crossoverP),
                                                                    classic_mut_chan(classicMutChan),
                                                                    sigma_mut_chan(sigmaMutChan),
                                                                    population_size(populationSize)
{
    if (gMax < gMax)
        throw runtime_error("Genetic_algorithm::Genetic_algorithm(): gMax must be <= gMin");
    if (classicMutChan + sigmaMutChan > 1)
        throw runtime_error("Genetic_algorithm::Genetic_algorithm(): classicMutChan + sigmaMutChan must be <= 1");
}

