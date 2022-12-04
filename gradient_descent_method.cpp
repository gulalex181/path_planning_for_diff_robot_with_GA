#include <iostream>
#include <cmath>
#include <vector>
#include <map>

#include "gradient_descent_method.h"

/* Optimization params */

namespace gradient_descent_method {

double alpha;
double tolerance;
int max_iterations;

std::vector<double> GD(
    std::map<std::string, std::string> settings,
    double df(int, std::vector<double>),
    std::vector<double> x_start
) {

    GD_set_settings(settings);

    int x_size = x_start.size();

    std::vector<double> x(x_size);
    std::vector<double> x_last(x_size);

    x_last = x_start;

    double distance;
    double distanceSquared;

    for (int i = 0; i < max_iterations; ++i) {

        distanceSquared = 0.0;

        std::cout << "iteration = " << i << ": ";

        for (int j = 0; j < x_size; ++j) {

            x.at(j) = x_last.at(j) - alpha * df(j, x_last);
            
            distanceSquared += std::pow(x.at(j) - x_last.at(j), 2);
            
            std::cout << "x" << j << " = " << x_last.at(j) << "\tdf = " << df(j, x_last) << "\t";

            x_last.at(j) = x.at(j);

        }

        distance = std::sqrt(distanceSquared);

        std::cout << "dist = " << distance << std::endl;

        if (distance < tolerance) break;

    }

    return x_last;

}


void GD_set_settings(
    std::map<std::string, std::string> settings
) {

    auto search = settings.find("alpha");
    
    if (search != settings.end()) {

        alpha = std::stod(search->second);

    }

    search = settings.find("tolerance");
    
    if (search != settings.end()) {

        tolerance = std::stod(search->second);

    }

    search = settings.find("max_iterations");
    
    if (search != settings.end()) {

        max_iterations = std::stoi(search->second);

    }

}

} // namespace gradient_descent_method