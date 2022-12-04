// To build, compile and run project:
// 1. remove `build` folder if it exists
// 2. $ mkdir build
// 3. $ cd build
// 4. $ cmake ..
// 5. $ make
// 6. $ ./object_control 

// To compile and run (in `build` folder):
// $ make && ./object_control

#include <vector>

/* Initial conditions */

std::vector<double> x_start = {10, 10, 0};

double t_start = 0.0;
double t_stop = 40.0;
double t_step = 0.05; // integration step

/* Terminal conditions */

std::vector<double> x_stop = {0, 0, 0};
double tolerance = 0.5;

/* Control settings */

double u_1_limit_min = -10;
double u_1_limit_max =  10;

double u_2_limit_min = -10;
double u_2_limit_max =  10;

int u_steps_count = 30;

std::vector<double> u_1_k;
std::vector<double> u_2_k;

/* Constraint settings */

double zone_fill_step = 0.2;
double zone_x_1_max = 100;
double zone_x_2_max = 100;

std::vector<std::vector<double>> restricted_zone_1;
std::vector<std::vector<double>> restricted_zone_2;

std::vector<std::vector<double>> result;

/* Function declarations */

std::vector<std::vector<double>> model_calculate(
    std::vector<double> x_start,
    std::vector<double> x_stop,
    double tolerance,
    double t_start,
    double t_stop
);

// Right parts of the differential equations
double f_1(double x_1, double x_2, double x_3, double u_1, double u_2, double t);
double f_2(double x_1, double x_2, double x_3, double u_1, double u_2, double t);
double f_3(double x_1, double x_2, double x_3, double u_1, double u_2, double t);

// Constraints on the state
bool h_1(double x_1, double x_2, double x_3);
bool h_2(double x_1, double x_2, double x_3);

std::vector<std::vector<double>> get_restricted_zone_points(
    bool h(double, double, double),
    double x_1_max,
    double x_2_max
);

// Control signals
double u_1(double k, double t_step, double u_last, double lower_limit, double upper_limit);
double u_2(double k, double t_step, double u_last, double lower_limit, double upper_limit);
// double u_1(double t);
// double u_2(double t);

// Functions for genetic algorithm
void on_best_chromosome_found(std::vector<double> best_chromosome_in_generation);
double calculate_distance(std::vector<double> from, std::vector<double> to);
double f_GA(std::vector<double> q);

// Functions for gradient descent
double f_GD(std::vector<double> x);
double df_GD(int arg_num, std::vector<double> x);