//
// Created by cheshirecat on 11/14/19.
//

#pragma once

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
    Genetic_algorithm(double gMin, double gMax, double (*fitnessCalc)(std::vector<double>&),
                      double classicMutChan, double sigmaMutChan, unsigned int populationSize,
                      unsigned int chromosomeSize, double crossoverCoef);

    [[nodiscard]] const std::vector<Individual> &getPopulation() const;
    [[nodiscard]] const std::vector<Individual *> &getChild() const;

    std::vector<double> study(unsigned int max_iteration, double good_error, double (*get_error)(std::vector<double>&));

private:
    //g_max >= g_min
    double g_max; // gene upper bound
    double g_min; // gene lower bound
    double (*fitness_calc) (std::vector<double>&);
    //classic_mut_chan + sigma_mut_chan <= 1
    double classic_mut_chan; //chance of classic mutation
    double sigma_mut_chan; //chance of sigma mutation;
    unsigned int population_size;
    std::vector<Individual> population;
    std::vector<Individual*> child; //child at last iteration
    unsigned int chromosome_size;
    double crossover_coef; //coefficent of crossover

    std::random_device rd;
    std::default_random_engine random_engine;

    double sum_fitness{}; //sum of fitness of population

    void initialization();
    void suitability_rating();
    void crossover();
    void calculate_sum_fitness();
    void selection();
    std::vector<double> generate_list_chances();
    //generate child for two parent; left parent gives left part of chromosome, right parent gives right part;
    std::vector<double> generate_chromosome(const std::vector<double>& left_parent, const std::vector<double>& right_parent,
                                            int crossover_point);
    //clear childer vector, but doesn't destroy objects witch it has.
    void clear_children();
    //return size of population to population_size
    void formation_new_population();
};
