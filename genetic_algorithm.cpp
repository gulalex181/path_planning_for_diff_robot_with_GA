#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>

#include "genetic_algorithm.h"

// std::default_random_engine e(std::time(NULL));
// std::uniform_int_distribution<int> random_gen(0, RAND_MAX);

// Random seed
std::random_device rd;

// Initialize Mersenne Twister pseudo-random number generator
std::mt19937 gen(rd());

// Generate pseudo-random numbers
// uniformly distributed in range
std::uniform_int_distribution<> dis(0, 1000);
// std::normal_distribution<> dis{0, 1000};
// std::student_t_distribution<> dis{10};

int REPRESENTATION_SIZE; // bits (integer + fractional part size)
int FRACTIONAL_SIZE;     // bits (fractional part size)
int POPULATION_SIZE;     // chromosomes
int MAX_GENERATIONS;     // number of generations to evolve
double XOVER_PROB;       // crossover probability
double MUTATION_PROB;    // mutation probability
int FUNCTION_ORDER;      // number of the function's params
int ELITIST_COUNT;       // number of the worst chromosomes to replace with the best one

std::vector<double> GA(
    std::map<std::string, std::string> settings,
    double f(std::vector<double>),
    void on_best_chromosome_found(std::vector<double>)
) {

    GA_set_settings(settings);
    
    struct GA_population* population = (struct GA_population*)calloc(1, sizeof(struct GA_population));
    struct GA_population* population_tmp = (struct GA_population*)calloc(1, sizeof(struct GA_population));
    
    GA_population_init(population, POPULATION_SIZE);
    GA_population_init(population_tmp, POPULATION_SIZE);

    GA_population_evaluate_fitness(population, f, -1.0);
    GA_population_select_best_chromosome(population, on_best_chromosome_found);

    // GA_population_print(population);
    GA_population_print_state(population);
    
    while (population->gen < MAX_GENERATIONS) {
    
        population->gen++;

        GA_population_apply_selection(population, population_tmp);
        GA_population_apply_crossover(population);
        GA_population_apply_mutation(population);
       
        GA_population_evaluate_fitness(population, f, -1.0);
        GA_population_select_best_chromosome(population, on_best_chromosome_found);

        // GA_population_print(population);
        GA_population_print_state(population);
        
        GA_population_apply_elitism(population, ELITIST_COUNT);

        if (population->chromosomes[POPULATION_SIZE].fitness == 0.0) {

            break;

        }
    
    }
    
    // GA_population_print(population);

    return GA_chromosome_binary_to_decimals(population->chromosomes[POPULATION_SIZE].binary);

}

void GA_set_settings(
    std::map<std::string, std::string> settings
) {

    auto search = settings.find("representation_size");
    
    if (search != settings.end()) {

        REPRESENTATION_SIZE = std::stoi(search->second);

    }

    search = settings.find("fractional_size");
    
    if (search != settings.end()) {

        FRACTIONAL_SIZE = std::stoi(search->second);

    }

    search = settings.find("population_size");
    
    if (search != settings.end()) {

        POPULATION_SIZE = std::stoi(search->second);

    }

    search = settings.find("max_generations");
    
    if (search != settings.end()) {

        MAX_GENERATIONS = std::stoi(search->second);

    }

    search = settings.find("xover_prob");
    
    if (search != settings.end()) {

        XOVER_PROB = std::stod(search->second);

    }

    search = settings.find("mutation_prob");
    
    if (search != settings.end()) {

        MUTATION_PROB = std::stod(search->second);

    }

    search = settings.find("function_order");
    
    if (search != settings.end()) {

        FUNCTION_ORDER = std::stoi(search->second);

    }

    search = settings.find("elitism_count");
    
    if (search != settings.end()) {

        ELITIST_COUNT = std::stoi(search->second);

    }

}

double GA_randomize(double min, double max) {

    // int random_int = std::rand();
    // int random_int = random_gen(e);
    int random_int = std::abs(static_cast<int>(dis(gen)));

    // std::cout << random_int << std::endl;

    // return (((double)(random_int % 100) / 100.0) * (max - min) + min);
    return (((double)(random_int % 100) / 100.0) * (max - min) + min);
    // return (((double)(random_int % 100) / 100.0) * (max - min) + min);

}

void GA_population_init(struct GA_population* population, int psize) {

    population->size = psize;
    population->gen = 0;
    population->chromosomes =
        (struct GA_chromosome*)calloc(population->size + 1, sizeof(struct GA_chromosome));

    for (int i = 0; i < population->size + 1; ++i) {

        GA_chromosome_init(&population->chromosomes[i]);

    }

    /**
     * Set the first chromosome in the population as the best one.
     * The last entry in the population is the best chromosome
     */

    population->best_chromosome_id = 0;
    GA_chromosome_copy(
        &population->chromosomes[POPULATION_SIZE],
        &population->chromosomes[0]
    );

}

void GA_population_print(struct GA_population* population) {

    for (int i = 0; i < population->size; ++i) {

        std::cout << "Chromosome #" << i << ": Binary = " << GA_chromosome_print(&population->chromosomes[i]) << std::endl
            << "Fitness = " << population->chromosomes[i].fitness << std::endl << std::endl;

    }

    std::cout << "Best chromosome: Binary = " << GA_chromosome_print(&population->chromosomes[population->size]) << std::endl
        << "Fitness = " << population->chromosomes[population->size].fitness << std::endl << std::endl;

}

void GA_chromosome_init(struct GA_chromosome* chromosome) {

    /* Create a binary chromosome representation with random genes */

    chromosome->binary = (int*)calloc(REPRESENTATION_SIZE * FUNCTION_ORDER, sizeof(int));

    for (int i = 0; i < REPRESENTATION_SIZE * FUNCTION_ORDER; ++i) {
    
        chromosome->binary[i] = ((GA_randomize(0, 1)) < 0.5) ? 0 : 1;

    }

    /* Fitness values */

    chromosome->fitness  = -INFINITY;

}

void GA_chromosome_copy(struct GA_chromosome* chromosome_to, struct GA_chromosome* chromosome_from) {

    chromosome_to->fitness = chromosome_from->fitness;

    for (int i = 0; i < REPRESENTATION_SIZE * FUNCTION_ORDER; ++i) {
    
        chromosome_to->binary[i] = chromosome_from->binary[i];

    }

}

std::string GA_chromosome_print(struct GA_chromosome* chromosome) {

    std::string binary_str = "";

    for (int i = 0; i < REPRESENTATION_SIZE * FUNCTION_ORDER; ++i) {
    
        binary_str += std::to_string(chromosome->binary[i]);

    }

    return binary_str;

}

std::vector<double> GA_chromosome_binary_to_decimals(int* binary) {

    std::vector<double> decimals;
    
    double decimal;

    int integer_part;
    int fractional_part;

    int integer_index;
    int fractional_index;

    int sign_index, start, end;

    int exponent;

    // std::cout << "--------------" << std::endl;

    for (int arg = 0; arg < FUNCTION_ORDER; ++arg) {

        sign_index = arg * REPRESENTATION_SIZE;
        start = arg * REPRESENTATION_SIZE + 1; // one bit for sign
        end = arg * REPRESENTATION_SIZE + REPRESENTATION_SIZE - FRACTIONAL_SIZE;

        // std::cout << "sign_index = " << sign_index << std::endl;
        // std::cout << "integer start = " << start << std::endl;
        // std::cout << "integer end = " << end << std::endl;

        integer_part = 0;
        fractional_part = 0;

        integer_index = 0;
        fractional_index = 0;

        /* Integer part */

        for (int i = start; i < end; ++i) {
        
            // example:
            // Let function arguments count = 2, fraction size = 4
            // Chromosome = 111110000001101101101110
            // Then first argument is 111110000001 where integer part = 11111000, fractional = 0001
            // Then second argument is 101101101110 where integer part = 10110110, fractional = 1110
            // 
            // For the first argument:
            // Absolute value of 11111000 is 1111000
            // Start index = 1
            // End index = 8
            // 
            // 6543210
            // 1111000 = 1 * 2^6 + 1 * 2^5 + 1 * 2^4 + 1 * 2^3 + 0 * 2^2 + 0 * 2^1 + 0 * 2^0
            exponent = (end - start - 1) - integer_index;
            integer_part += binary[i] * std::pow(2, exponent);
            // std::cout << "i = " << i << ", binary[i] = " << binary[i] << ", exponent = " << exponent << ", integer_part = " << integer_part << std::endl;

            ++integer_index;
        
        }

        // std::cout << "integer part = " << integer_part << std::endl;

        /* Fractional part */

        start = arg * REPRESENTATION_SIZE + REPRESENTATION_SIZE - FRACTIONAL_SIZE;
        end = arg * REPRESENTATION_SIZE + REPRESENTATION_SIZE;

        // std::cout << "fractional start = " << start << std::endl;
        // std::cout << "fractional end = " << end << std::endl;

        for (int i = start; i < end; ++i) {
        
            exponent = (end - start - 1) - fractional_index;
            fractional_part += binary[i] * std::pow(2, exponent);

            ++fractional_index;
        
        }

        // std::cout << "fractional_part = " << fractional_part << std::endl;

        decimal = (double)integer_part + std::stod("0." + std::to_string(fractional_part));

        // std::cout << "decimal = " << decimal << std::endl;

        if (binary[sign_index] == 1) {

            decimal *= -1;

        }

        // std::cout << "decimal (+ sign) = " << decimal << std::endl;

        decimals.push_back(decimal);

    }

    // std::cout << "--------------" << std::endl;

    return decimals;

}

void GA_population_evaluate_fitness(
    struct GA_population* population,
    double f(std::vector<double>),
    double weight
) {

    for (int i = 0; i < population->size; ++i) {

        population->chromosomes[i].fitness =
            weight * f(GA_chromosome_binary_to_decimals(population->chromosomes[i].binary));

    }

}

void GA_population_select_best_chromosome(
    struct GA_population* population,
    void on_best_chromosome_found(std::vector<double>)
) {

    /* Search the best chromosome */

    for (int i = 1; i < population->size; ++i) {
        
        if (population->chromosomes[i].fitness >
            population->chromosomes[POPULATION_SIZE].fitness) {
        
            population->best_chromosome_id = i;
            GA_chromosome_copy(
                &population->chromosomes[POPULATION_SIZE],
                &population->chromosomes[population->best_chromosome_id]
            );
        
        }

    }

    on_best_chromosome_found(GA_chromosome_binary_to_decimals(population->chromosomes[POPULATION_SIZE].binary));

}

void GA_population_print_state(struct GA_population* population) {

    std::cout << "Generation #" << population->gen
        << ": Best fitness = " << population->chromosomes[POPULATION_SIZE].fitness << std::endl;

}

void GA_population_apply_selection(struct GA_population* population, struct GA_population* population_tmp) {

    int index_1, index_2, index_3;

    /* Tournament selection */

    // std::cout << "--------------" << std::endl;
    // std::cout << "SELECTION [BEFORE]" << std::endl;
    // GA_population_print(population);
    // std::cout << "--------------" << std::endl;

    for (int i = 0; i < population->size; ++i) {

        index_1 = 0;
        index_2 = 0;
        index_3 = 0;

        while (index_1 == index_2 || index_1 == index_3 || index_2 == index_3) {

            index_1 = GA_randomize(0, POPULATION_SIZE - 1);
            index_2 = GA_randomize(0, POPULATION_SIZE - 1);
            index_3 = GA_randomize(0, POPULATION_SIZE - 1);

        }

        // std::cout << "Index_1 = " << index_1 << ", Index_2 = " << index_2 << ", Index_3 = " << index_3 << std::endl;

        if (population->chromosomes[index_1].fitness > population->chromosomes[index_2].fitness) {

            if (population->chromosomes[index_1].fitness > population->chromosomes[index_3].fitness) {

                GA_chromosome_copy(
                    &population_tmp->chromosomes[i],
                    &population->chromosomes[index_1]
                );

            } else {

                GA_chromosome_copy(
                    &population_tmp->chromosomes[i],
                    &population->chromosomes[index_3]
                );

            }

        } else {

            if (population->chromosomes[index_2].fitness > population->chromosomes[index_3].fitness) {
            
                GA_chromosome_copy(
                    &population_tmp->chromosomes[i],
                    &population->chromosomes[index_2]
                );

            } else {

                GA_chromosome_copy(
                    &population_tmp->chromosomes[i],
                    &population->chromosomes[index_3]
                );

            }

        }
    
    }

    // Once the new population is created copy it back in the working var
    for (int i = 0 ; i < population->size; ++i) {

        population->chromosomes[i] = population_tmp->chromosomes[i];
    
    }

    // std::cout << "--------------" << std::endl;
    // std::cout << "SELECTION [AFTER]" << std::endl;
    // GA_population_print(population);
    // std::cout << "--------------" << std::endl;

}

void GA_population_apply_crossover(struct GA_population* population) {

    // Counter of members chosen
    int cnt = 0;

    // Crossover point index
    int crossover_point;

    // Crossover probability
    double crossover_probability = 0.0f;

    // Parent container
    struct GA_chromosome* p1;

    // std::cout << "--------------" << std::endl;
    // std::cout << "CROSSOVER [BEFORE]" << std::endl;
    // GA_population_print(population);
    // std::cout << "--------------" << std::endl;

    // Crossover loop
    for (int i = 0; i < population->size; ++i) {
    
        crossover_probability = GA_randomize(0, 1);
    
        // std::cout << "crossover_probability = " << crossover_probability << std::endl;

        if (crossover_probability < XOVER_PROB) {

            cnt++;
            
            if (cnt % 2 == 0) {

                crossover_point = floor(GA_randomize(2, REPRESENTATION_SIZE * FUNCTION_ORDER - 3));

                // std::cout << "CROSS POINT INDEX = " << crossover_point << std::endl;

                int tmp[crossover_point];

                for (int j = 0; j < REPRESENTATION_SIZE * FUNCTION_ORDER; ++j) {

                    if (j < crossover_point) {
                    
                        tmp[j] = p1->binary[j];
                        p1->binary[j] = population->chromosomes[i].binary[j];
                        population->chromosomes[i].binary[j] = tmp[j];

                    }


                }

            } else {

                p1 = &population->chromosomes[i];
            
            }
        
        }
    
    }

    // std::cout << "--------------" << std::endl;
    // std::cout << "CROSSOVER [AFTER]" << std::endl;
    // GA_population_print(population);
    // std::cout << "--------------" << std::endl;

    // std::cout << "cnt = " << cnt << std::endl;

    // for (int i = 0; i < population->size; ++i) {

    //     population->chromosomes[i].decimal = GA_chromosome_encode_binary(population->chromosomes[i].binary);

    // }
}

void GA_population_apply_mutation(struct GA_population* population) {

    double mutation_probability = 0.0f;

    // std::cout << "--------------" << std::endl;
    // std::cout << "MUTATION [BEFORE]" << std::endl;
    // GA_population_print(population);
    // std::cout << "--------------" << std::endl;

    for (int i = 0; i < population->size; ++i) {
    
        for (int j = 0; j < REPRESENTATION_SIZE * FUNCTION_ORDER; ++j) {

            mutation_probability = GA_randomize(0, 1);

            // std::cout << "mutation_probability = " << mutation_probability << std::endl;

            if (mutation_probability < MUTATION_PROB) {

                population->chromosomes[i].binary[j] = 1 - population->chromosomes[i].binary[j];

            }

        }

    }

    // std::cout << "--------------" << std::endl;
    // std::cout << "MUTATION [AFTER]" << std::endl;
    // GA_population_print(population);
    // std::cout << "--------------" << std::endl;

}

void GA_population_apply_elitism(struct GA_population* population, int count) {

    struct GA_chromosome* best  = (struct GA_chromosome*)calloc(1, sizeof(struct GA_chromosome));
    struct GA_chromosome* worst = (struct GA_chromosome*)calloc(1, sizeof(struct GA_chromosome));

    int best_id = 0, worst_id = 0;

    GA_chromosome_init(best);
    GA_chromosome_init(worst);

    best->fitness  = -INFINITY;
    worst->fitness = INFINITY;

    /* Search the best and the worst chromosomes */

    // std::cout << "--------------" << std::endl;
    // std::cout << "ELITIST [BEFORE]" << std::endl;
    // GA_population_print(population);
    // std::cout << "--------------" << std::endl;

    for (int i = 0; i < population->size - 1; ++i) {

        if (population->chromosomes[i].fitness > population->chromosomes[i + 1].fitness) {

            if (population->chromosomes[i].fitness > best->fitness) {

                best->fitness = population->chromosomes[i].fitness;
                best_id = i;
            
            }
            
            if (population->chromosomes[i + 1].fitness < worst->fitness) {
            
                worst->fitness = population->chromosomes[i + 1].fitness;
                worst_id = i + 1;
            
            }
        
        } else {

            if (population->chromosomes[i].fitness < worst->fitness) {

                worst->fitness = population->chromosomes[i].fitness;
                worst_id = i;
            
            }
            
            if (population->chromosomes[i + 1].fitness > best->fitness) {

                best->fitness = population->chromosomes[i + 1].fitness;
                best_id = i + 1;
            
            }
        
        }
    
    }

    /* if best chromosome from the new population is better than */
    /* the best chromosome from the previous population, then    */
    /* copy the best from the new population; else replace the   */
    /* worst chromosome from the current population with the     */
    /* best one from the previous generation                     */

    if (best->fitness > population->chromosomes[POPULATION_SIZE].fitness) {
    
        GA_chromosome_copy(
            &population->chromosomes[POPULATION_SIZE],
            &population->chromosomes[best_id]
        );
    
    } else {

        GA_chromosome_copy(
            &population->chromosomes[worst_id],
            &population->chromosomes[POPULATION_SIZE]
        );

    }

    /**
     * Replace `count` the worst chromosomes with the best one
     */
    
    worst->fitness = INFINITY;

    for (int j = 0; j < count; ++j) {

        for (int i = 0; i < population->size - 1; ++i) {

            if (population->chromosomes[i].fitness < worst->fitness) {
            
                worst->fitness = population->chromosomes[i].fitness;
                worst_id = i;
            
            }
        
        }

        GA_chromosome_copy(
            &population->chromosomes[worst_id],
            &population->chromosomes[POPULATION_SIZE]
        );

    }


    // std::cout << "--------------" << std::endl;
    // std::cout << "ELITIST [AFTER]" << std::endl;
    // GA_population_print(population);
    // std::cout << "--------------" << std::endl;

}