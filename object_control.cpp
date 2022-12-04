#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <unistd.h>

#include "genetic_algorithm.h"
#include "gradient_descent_method.h"
#include "graphics.h"
#include "object_control.h"

int main() {

    u_1_k.resize(2 * u_steps_count);
    u_2_k.resize(2 * u_steps_count);

    restricted_zone_1 = get_restricted_zone_points(h_1, zone_x_1_max, zone_x_2_max);
    restricted_zone_2 = get_restricted_zone_points(h_2, zone_x_1_max, zone_x_2_max);

    /* Find the shortest path */

    GA({
        {"representation_size", "7"},
        {"fractional_size", "3"},
        {"population_size", "400"},
        {"max_generations", "200"},
        {"xover_prob", "0.7"},
        {"mutation_prob", "0.6"},
        {"function_order", std::to_string(2 * u_steps_count)},
        {"elitism_count", "50"},
    }, f_GA, on_best_chromosome_found);

    // // std::vector<double> a = {4.5, 0.5, 5.6, 10, 17};
    // // std::vector<double> a = {4.5, 0.5};
    // std::vector<double> a = {4.5};

    // /* Start point */

    // // std::vector<double> q_start = {20, 10, 5, 1, 4};
    // // std::vector<double> q_start = {20, 10};
    // std::vector<double> q_start = {20};

    // GD({
    //     {"alpha", "0.1"},
    //     {"tolerance", "0.00000001"},
    //     {"max_iterations", "100"},
    // }, df_GD, x_start);

    /* Plot results */

    // Control signals

    plot_2d<double>(result.at(0), result.at(4), {
        {"legend", "u_1(t)"},
        {"color", "green"},
        {"linestyle", "-"},
        {"markerstyle", ""},
        {"hold", "on"},
    });

    plot_2d<double>(result.at(0), result.at(5), {
        {"legend", "u_2(t)"},
        {"color", "blue"},
        {"linestyle", "-"},
        {"markerstyle", ""},
        {"hold", "off"},
    });

    // State variables

    plot_2d<double>(result.at(0), result.at(1), {
        {"legend", "x_1(t)"},
        {"color", "red"},
        {"linestyle", "-"},
        {"markerstyle", ""},
        {"hold", "on"},
    });

    plot_2d<double>(result.at(0), result.at(2), {
        {"legend", "x_2(t)"},
        {"color", "green"},
        {"linestyle", "-"},
        {"markerstyle", ""},
        {"hold", "on"},
    });

    plot_2d<double>(result.at(0), result.at(3), {
        {"legend", "x_3(t)"},
        {"color", "blue"},
        {"linestyle", "-"},
        {"markerstyle", ""},
        {"hold", "off"},
    });

    // Path with restricted zones
    
    plot_2d<double>(result.at(1), result.at(2), {
        {"legend", "Path"},
        {"color", "green"},
        {"linestyle", "-"},
        {"markerstyle", ""},
        {"hold", "on"},
    });

    plot_2d<double>(restricted_zone_1.at(0), restricted_zone_1.at(1), {
        {"legend", "Restricted zone 1"},
        {"color", "red"},
        {"linestyle", ""},
        {"markerstyle", "."},
        {"hold", "on"},
    });

    plot_2d<double>(restricted_zone_2.at(0), restricted_zone_2.at(1), {
        {"legend", "Restricted zone 2"},
        {"color", "red"},
        {"linestyle", ""},
        {"markerstyle", "."},
        {"hold", "off"},
    });

    return 0;

}

void on_best_chromosome_found(std::vector<double> best_chromosome_in_generation) {

    int middle_index = best_chromosome_in_generation.size() / 2 - 1;

    for (int i = 0; i <= middle_index; ++i) {

        u_1_k.at(i) = best_chromosome_in_generation.at(i);
        u_2_k.at(i) = best_chromosome_in_generation.at(middle_index + 1 + i);

    }

    result = model_calculate(
        x_start,
        x_stop,
        tolerance,
        t_start,
        t_stop
    );

    std::vector<double> t   = result.at(0);
    std::vector<double> x_1 = result.at(1);
    std::vector<double> x_2 = result.at(2);
    std::vector<double> x_3 = result.at(3);

    int last_index = t.size() - 1;

    double elapsed_time = t.at(last_index);
    double position_error = calculate_distance(
        {x_1.at(last_index), x_2.at(last_index), x_3.at(last_index)},
        x_stop
    );

    int points_in_restricted_zone = 0;

    for (int i = 0; i <= last_index; ++i) {

        points_in_restricted_zone += (h_1(x_1.at(i), x_2.at(i), x_3.at(i)) + h_2(x_1.at(i), x_2.at(i), x_3.at(i)));

    }

    std::cout << "elapsed_time = " << (elapsed_time - t_step) << ", ";
    std::cout << "position_error = " << position_error << ", ";
    std::cout << "points_in_restricted_zone = " << points_in_restricted_zone << ": ";

    plot_2d<double>(x_1, x_2, {
        {"legend", "Path"},
        {"color", "green"},
        {"linestyle", "-"},
        {"markerstyle", ""},
        // {"hold", "on"},
    });

    // plot_2d<double>(restricted_zone_1.at(0), restricted_zone_1.at(1), {
    //     {"legend", "Restricted zone 1"},
    //     {"color", "red"},
    //     {"linestyle", ""},
    //     {"markerstyle", "."},
    //     {"hold", "on"},
    // });

    // plot_2d<double>(restricted_zone_2.at(0), restricted_zone_2.at(1), {
    //     {"legend", "Restricted zone 2"},
    //     {"color", "red"},
    //     {"linestyle", ""},
    //     {"markerstyle", "."},
    // });

    // sleep(1);

}

double calculate_distance(std::vector<double> from, std::vector<double> to) {

    double distance = 0.0;

    for (int i = 0; i < to.size(); ++i) {

        distance += std::pow(to.at(i) - from.at(i), 2.0);

    }

    return std::sqrt(distance);

}

double f_GA(std::vector<double> k) {

    // u_i = k_i * delta_t + u_{i - 1}
    // GA gives us K coeffecient for u1 and u2

    int middle_index = k.size() / 2 - 1;

    for (int i = 0; i <= middle_index; ++i) {

        u_1_k.at(i) = k.at(i);
        u_2_k.at(i) = k.at(middle_index + 1 + i);

    }

    result = model_calculate(
        x_start,
        x_stop,
        tolerance,
        t_start,
        t_stop
    );

    std::vector<double> t   = result.at(0);
    std::vector<double> x_1 = result.at(1);
    std::vector<double> x_2 = result.at(2);
    std::vector<double> x_3 = result.at(3);

    int last_index = t.size() - 1;

    double elapsed_time = t.at(last_index);
    
    double position_error = calculate_distance(
        {x_1.at(last_index), x_2.at(last_index), x_3.at(last_index)},
        x_stop
    );

    int points_in_restricted_zone = 0;

    for (int i = 0; i <= last_index; ++i) {

        points_in_restricted_zone += (h_1(x_1.at(i), x_2.at(i), x_3.at(i)) + h_2(x_1.at(i), x_2.at(i), x_3.at(i)));

    }

    // std::cout << "elapsed_time" << elapsed_time << std::endl;
    // std::cout << "position_error" << 20 * position_error << std::endl;
    // std::cout << "points_in_restricted_zone" << 10 * points_in_restricted_zone << std::endl;

    if (elapsed_time < t_stop * 0.5) {

        return elapsed_time + position_error + 100 * points_in_restricted_zone;

    } else {

        return 10 * elapsed_time + position_error + 100 * points_in_restricted_zone;

    }

}

double df_GD(int arg_num, std::vector<double> k) {

    return 0.0;

}

double u_1(double k, double t_step, double u_last, double lower_limit, double upper_limit) {

    double u = k * t_step + u_last;

    // std::cout << "u_1 = " << u << std::endl;

    /* Control limits */

    if (u <= lower_limit) {

        u = lower_limit;

    }

    if (u >= upper_limit) {

        u = upper_limit;

    }

    return u;

}

double u_2(double k, double t_step, double u_last, double lower_limit, double upper_limit) {

    double u = k * t_step + u_last;

    // std::cout << "u_2 = " << u << std::endl;

    /* Control limits */

    if (u <= lower_limit) {

        u = lower_limit;

    }

    if (u >= upper_limit) {

        u = upper_limit;

    }

    return u;

}

std::vector<std::vector<double>> model_calculate(
    std::vector<double> x_start,
    std::vector<double> x_stop,
    double tolerance,
    double t_start,
    double t_stop
) {

    // State variables
    std::vector<double> x_1, x_2, x_3;

    // Control variables
    std::vector<double> u_1_, u_2_;

    // Time variable
    std::vector<double> t;

    /* Add initial conditionals to the vector */

    x_1.push_back(x_start.at(0));
    x_2.push_back(x_start.at(1));
    x_3.push_back(x_start.at(2));
    u_1_.push_back(0.0);
    u_2_.push_back(0.0);

    t.push_back(t_start);

    /* Compute variables */

    // Set iteration counter to zero
    int i = 0;
    double t_step_for_u = (t_stop - t_start) / u_steps_count;
    int u_step = 0;

    // Start loop
    while (t[i] <= t_stop
      && (std::abs(x_1[i] - x_stop.at(0)) > tolerance
       || std::abs(x_2[i] - x_stop.at(1)) > tolerance
       || std::abs(x_3[i] - x_stop.at(2)) > tolerance)) {
    // while (t[i] <= t_stop) {
        
        /* Increase size of the vector by one */

        x_1.resize(x_1.size() + 1);
        x_2.resize(x_2.size() + 1);
        x_3.resize(x_3.size() + 1);

        u_1_.resize(u_1_.size() + 1);
        u_2_.resize(u_2_.size() + 1);

        t.resize(t.size() + 1);

        /* Control signals */

        u_step = (int)(t[i] / t_step_for_u);

        u_1_.at(i + 1) = u_1(u_1_k[u_step], t_step, u_1_.at(i), u_1_limit_min, u_1_limit_max);
        u_2_.at(i + 1) = u_2(u_2_k[u_step], t_step, u_2_.at(i), u_2_limit_min, u_2_limit_max);

        /* Numerical method */

        double f_1_current = f_1(
            x_1[i],
            x_2[i],
            x_3[i],
            u_1_.at(i),
            u_2_.at(i),
            t[i]
        );

        x_1[i + 1] = x_1[i] + (t_step / 2) * (
            f_1_current
            + f_1(
                x_1[i]    + t_step * f_1_current,
                x_2[i]    + t_step * f_1_current,
                x_3[i]    + t_step * f_1_current,
                u_1_.at(i) + t_step * f_1_current,
                u_2_.at(i) + t_step * f_1_current,
                t[i] + t_step
            )
        );

        double f_2_current = f_2(
            x_1[i],
            x_2[i],
            x_3[i],
            u_1_.at(i),
            u_2_.at(i),
            t[i]
        );

        x_2[i + 1] = x_2[i] + (t_step / 2) * (
            f_2_current
            + f_2(
                x_1[i]    + t_step * f_2_current,
                x_2[i]    + t_step * f_2_current,
                x_3[i]    + t_step * f_2_current,
                u_1_.at(i) + t_step * f_2_current,
                u_2_.at(i) + t_step * f_2_current,
                t[i] + t_step
            )
        );

        double f_3_current = f_3(
            x_1[i],
            x_2[i],
            x_3[i],
            u_1_.at(i),
            u_2_.at(i),
            t[i]
        );

        x_3[i + 1] = x_3[i] + (t_step / 2) * (
            f_3_current
            + f_3(
                x_1[i]    + t_step * f_3_current,
                x_2[i]    + t_step * f_3_current,
                x_3[i]    + t_step * f_3_current,
                u_1_.at(i) + t_step * f_3_current,
                u_2_.at(i) + t_step * f_3_current,
                t[i] + t_step
            )
        );

        // x_1[i + 1] = modified_euler_method(x_1[i], x_2[i], x_3[i], u_1[i], u_2[i], t[i], t_step, &f_1);
        // x_2[i + 1] = modified_euler_method(x_1[i], x_2[i], x_3[i], u_1[i], u_2[i], t[i], t_step, &f_2);
        // x_3[i + 1] = modified_euler_method(x_1[i], x_2[i], x_3[i], u_1[i], u_2[i], t[i], t_step, &f_3);
        
        t[i + 1] = t[i] + t_step;

        /* Print results */

        // std::cout << "x_1 = " << x_1[i] << ",\tx_2 = " << x_2[i] << ",\tx_3 = " << x_3[i] << ",\tt = " << t[i] << std::endl;
        
        // Increase the iteration counter
        ++i;

    }

    /* Print results */

    // std::cout << "x_1 = " << x_1[i] << ",\tx_2 = " << x_2[i] << ",\tx_3 = " << x_3[i] << ",\tt = " << t[i] << std::endl;

    return {
        t,
        x_1,
        x_2,
        x_3,
        u_1_,
        u_2_
    };

}

double f_1(double x_1, double x_2, double x_3, double u_1, double u_2, double t) {

    return (0.5 * (u_1 + u_2) * std::cos(x_3));

}

double f_2(double x_1, double x_2, double x_3, double u_1, double u_2, double t) {

    return (0.5 * (u_1 + u_2) * std::sin(x_3));

}

double f_3(double x_1, double x_2, double x_3, double u_1, double u_2, double t) {

    return (0.5 * (u_1 - u_2));

}

bool h_1(double x_1, double x_2, double x_3) {

    return !((2.5 - std::sqrt(std::pow(x_1 - 2.5, 2) + std::pow(x_2 - 2.5, 2))) <= 0);

}

bool h_2(double x_1, double x_2, double x_3) {

    return !((2.5 - std::sqrt(std::pow(x_1 - 7.5, 2) + std::pow(x_2 - 7.5, 2))) <= 0);

}

std::vector<std::vector<double>> get_restricted_zone_points(
    bool h(double, double, double),
    double x_1_max,
    double x_2_max
) {

    std::vector<std::vector<double>> restricted_zone;
    std::vector<double> x_1_in_zone;
    std::vector<double> x_2_in_zone;

    for (double i = -x_1_max; i < x_1_max; i += zone_fill_step) {

        for (double j = -x_2_max; j < x_2_max; j += zone_fill_step) {
        
            // if a point is in one of restricted zones
            if (h(i, j, 0)) {

                x_1_in_zone.push_back(i);
                x_2_in_zone.push_back(j);

            }

        }

    }

    restricted_zone.push_back(x_1_in_zone);
    restricted_zone.push_back(x_2_in_zone);

    return restricted_zone;

}
