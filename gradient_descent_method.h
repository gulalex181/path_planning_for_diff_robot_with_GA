#include <vector>

namespace gradient_descent_method {

/* Function declarations */

// Numerical optimization methods
std::vector<double> GD(
    std::map<std::string, std::string> settings,
    double df(int, std::vector<double>),
    std::vector<double> x_start
);

void GD_set_settings(
    std::map<std::string, std::string> settings
);

} // namespace gradient_descent_method