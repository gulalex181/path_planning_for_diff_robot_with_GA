#include <vector>

/* Structures */

// Chromosome abstraction
struct GA_chromosome {

    // Binary chromosome representation
    int* binary;

    // // Converted binary to decimal value
    // double decimal;

    // The fitness value of the chromosome
    double fitness;

    // // Relative fitness
    // double rfitness;

    // // Cumulative fitness
    // double cfitness;

};

// Chromosomes population
struct GA_population {

    // Size in chromosomes
    int size;

    // Chromosomes in the population
    struct GA_chromosome* chromosomes;

    // Current evolution generation
    int gen;

    // Index of the fittest chromosome
    int best_chromosome_id;

};

/* Function declarations */

/**
 * \brief Numerical genetic optimization algorithms
 */
std::vector<double> GA(
    std::map<std::string, std::string> settings,
    double f(std::vector<double>),
    void on_best_chromosome_found(std::vector<double>)
);

/**
 * \brief Set settings of the algorithm
 */
void GA_set_settings(
    std::map<std::string, std::string> settings
);

/**
 * \brief Initialize a chromosome
 */
void GA_chromosome_init(struct GA_chromosome* chromosome);

/**
 * \brief Copy a chromosome
 */
void GA_chromosome_copy(struct GA_chromosome* chromosome_to, struct GA_chromosome* chromosome_from);

/**
 * \brief Print a chromosome
 */
std::string GA_chromosome_print(struct GA_chromosome* chromosome);

/**
 * \brief Initialize a chromosomes population
 * with given parameters
 */
void GA_population_init(struct GA_population* population, int psize);

/**
 * \brief Print a population
 */
void GA_population_print(struct GA_population* population);

/**
 * \brief Encoding function from binary
 * to decimal given a binary string
 */
// double GA_chromosome_encode_binary(int *b);
std::vector<double> GA_chromosome_binary_to_decimals(int* binary);

/**
 * \brief Fitness evaluation function, takes a user defined
 * function and computes fitness of every chromosome as product of
 * that function and wight. By default genetic algorithm searches
 * the chromosomes with the biggest fitness. And if we need to minimize
 * a function we need to use a weight (e.g. -1).
 * 
 * \param population Population of chromosomes
 * \param f Function to maximize
 * \param weight Weight the function needs to be multiplied by
 */
void GA_population_evaluate_fitness(
    struct GA_population* population,
    double f(std::vector<double>),
    double weight
);

/**
 * \brief Select the best (fittest) chromosome
 * in the population
 */
void GA_population_select_best_chromosome(
    struct GA_population* population,
    void on_best_chromosome_found(std::vector<double>)
);

/**
 * \brief Print the state of the current evolution generation
 */
void GA_population_print_state(struct GA_population* population);

/**
 * \brief selection function using the
 * elitist model in which only the best
 * chromosome survives - Winner-Take-All
 */
void GA_population_apply_selection(struct GA_population* p, struct GA_population* newp);

/**
 * \brief Apply the single point crossover
 * operator which takes 2 parents
 */
void GA_population_apply_crossover(struct GA_population* population);

/**
 * \brief Apply mutation - random uniform mutation of the genes
 */
void GA_population_apply_mutation(struct GA_population* population);

/**
 * \brief Random value generator within the bounds
 */
double GA_randomize(double min, double max);

/**
 * \brief Apply elitism so that if the previous
 * best chromosome is better than the current
 * generation best the first will replace
 * the worst chromosome in the current generation.
 */
void GA_population_apply_elitism(struct GA_population* population, int count);
