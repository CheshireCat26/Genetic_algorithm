//
// Created by cheshirecat on 11/14/19.
//

#include "Genetic_algorithm.h"
#include <stdexcept>
#include <algorithm>
#include <random>
#include "my_utility.h"

using namespace std;

Genetic_algorithm::Genetic_algorithm(double gMax, double gMin, double (*fitnessCalc)(vector<double>&) ,
                                     double classicMutChan, double sigmaMutChan,
                                     unsigned int populationSize, unsigned int chromosomeSize, double crossoverCoef)
                                     : g_max(gMax), g_min(gMin), fitness_calc(fitnessCalc),
                                       classic_mut_chan(classicMutChan), sigma_mut_chan(sigmaMutChan),
                                       population_size(populationSize), chromosome_size(chromosomeSize),
                                       crossover_coef{crossoverCoef}
{
    random_engine.seed(rd());
    if (gMax < gMax)
        throw runtime_error("Genetic_algorithm::Genetic_algorithm(): gMax must be <= gMin");
    if (classicMutChan + sigmaMutChan > 1)
        throw runtime_error("Genetic_algorithm::Genetic_algorithm(): classicMutChan + sigmaMutChan must be <= 1");
    if (crossover_coef > 1)
        throw runtime_error("Genetic_algorithm::Genetic_algorithm(): crossover_coef must be <= 1");
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

const vector<Individual> &Genetic_algorithm::getPopulation() const
{
    return population;
}

void Genetic_algorithm::suitability_rating()
{
    for (Individual& individual : population)
        individual.fitness = fitness_calc(individual.chromosome);
}

std::vector<double> Genetic_algorithm::study(unsigned int max_iteration, double good_error,
                                             double (*get_error)(std::vector<double> &))
{
    //1 step
    initialization();
    //2 step
    suitability_rating();

    //prepare for calculate min error
    sort(population.begin(), population.end(),
            [](const Individual& a, const Individual& b){ return a.fitness > b.fitness;});
    calculate_sum_fitness();
    unsigned  int cur_iteration{1};
    //calc min error
    double cur_min_error = get_error(population[0].chromosome);
    //loop
    while (cur_iteration <= max_iteration && cur_min_error > good_error)
    {
        //1 loop step
        selection();
        //2 loop step
        crossover();
        //3 loop step
        mutation();
        //4 loop step
        suitability_rating();
        //5 loop step
        formation_new_population();
        //prepare for calculate min error
        sort(population.begin(), population.end(),
                [](const Individual& a, const Individual& b){ return a.fitness > b.fitness;});
        //calc min error
        cur_min_error = get_error(population[0].chromosome);
        //end loop operations
        calculate_sum_fitness();
        cur_iteration++;
    }

    return population[0].chromosome;
}

void Genetic_algorithm::crossover()
{
    vector<double> chances = generate_list_chances();
    discrete_distribution parent_distribution(chances.begin(), chances.end());
    uniform_int_distribution<int> crossover_p_distribution(1, population[0].chromosome.size() - 1); //generate crossover points
    for (int i{0}; i < round(population_size * crossover_coef); i++)
    {
        //choice 2 parents
        Individual& parent1 = population[parent_distribution(random_engine)];
        Individual& parent2 = population[parent_distribution(random_engine)];
        //generate crossover_points
        int crossover_point = crossover_p_distribution(random_engine);
        //create first child
        Individual child1 = Individual(parent1.chromosome.size());
        child1.chromosome = generate_chromosome(parent1.chromosome, parent2.chromosome, crossover_point);
        child1.child = true;
        //create second child
        Individual child2 = Individual(parent1.chromosome.size());
        child2.chromosome = generate_chromosome(parent2.chromosome, parent1.chromosome, crossover_point);
        child2.child = true;

        population.emplace_back(child1);

        population.emplace_back(child2);
    }

}

void Genetic_algorithm::calculate_sum_fitness()
{
    sum_fitness = 0;
    for (Individual& individual : population)
        sum_fitness += individual.fitness;
}

//calc chances
void Genetic_algorithm::selection()
{
    //calc sigma
    double aver_fitness = sum_fitness/population.size();
    double a = 2.5;
    double sum{0};
    for (auto & i : population)
        sum += i.fitness - aver_fitness;
    double sigma = sqrt(sum/(population.size() - 1));

    for (Individual& individual : population)
    {
        //calc sigma fitness
        double sigma_fitness = individual.fitness + (aver_fitness - a * sigma);
        if (sigma_fitness < 0)
            sigma_fitness = 0;
        individual.chance = sigma_fitness / sum_fitness;
    }

}

std::vector<double> Genetic_algorithm::generate_list_chances()
{
    vector<double> chances;
    for (const Individual& individual : population)
        chances.emplace_back(individual.chance);
    return  chances;
}

std::vector<double> Genetic_algorithm::generate_chromosome(const std::vector<double> &left_parent,
                                                           const std::vector<double> &right_parent,int crossover_point)
{
    std::vector<double> chrom(left_parent.size());
    for (int j{0}; j < crossover_point; j++)
        chrom[j] = left_parent[j];
    for (int j{crossover_point}; j < left_parent.size(); j++)
        chrom[j] = right_parent[j];
    return chrom;
}

void Genetic_algorithm::formation_new_population()
{
    while (population.size() > population_size)
        population.erase(population.end() - 1);
}

void Genetic_algorithm::mutation()
{
    vector<Individual *> child;
    for (Individual& individual : population)
        if (individual.child)
            child.emplace_back(&individual);

    classic_mutation(child);
    sigma_mitation(child);

    for(Individual& individual1 : population)
        if (individual1.child)
            individual1.child = false;
}

void Genetic_algorithm::classic_mutation(std::vector<Individual*> child)
{
    const double small_value = 0.001;
    uniform_int_distribution distribution(1, static_cast<int>(child[0]->chromosome.size()));
    uniform_int_distribution bin_distribution(0, 1);
    for (int i{0}; i < round(child.size() * classic_mut_chan); i++)
    {
        Individual& mut_child = get_random_child(child);
        //mutation of random number of genes
        int number_genes = distribution(random_engine);
        vector<int> id_genes;
        for(int j{0}; j < number_genes; j++)
        {
            id_genes.emplace_back(get_random_gene(id_genes));
            if (bin_distribution(random_engine))
                mut_child.chromosome[id_genes[id_genes.size() - 1]] += small_value;
            else
                mut_child.chromosome[id_genes[id_genes.size() - 1]] -= small_value;
        }
    }
}

Individual &Genetic_algorithm::get_random_child(std::vector<Individual*>& child)
{
    uniform_int_distribution distribution(0, static_cast<int>(child.size() - 1));
    while (true)
    {
        int id_child = distribution(random_engine);
        if (child[id_child] != nullptr)
        {
            Individual* ret = child[id_child];
            child[id_child] = nullptr;
            return *ret;
        }
    }
}

int Genetic_algorithm::get_random_gene(std::vector<int> taken_gene) {
    uniform_int_distribution distribution(0, static_cast<int>(population[0].chromosome.size() - 1));
    while (true)
    {
        int id_gene = distribution(random_engine);
        bool taken = false;
        for (int i = 0; i < taken_gene.size() && !taken; i++)
            if (taken_gene[i] == id_gene)
                taken = true;

        if (!taken)
            return  id_gene;
    }
}

void Genetic_algorithm::sigma_mitation(std::vector<Individual*> child)
{
    //TODO implement
}

