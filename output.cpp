#include "output.h"
#include "analytical.h"
#include "params.h"
#include <fstream>
#include <stdexcept>
#include <cmath>

void save_profile(const std::string& filename,
                 const std::vector<double>& u,
                 double dx,
                 double time) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << "# x\tu_analytical\tu_numerical\terror\n";
    double x_start = -params::a;
    int N = static_cast<int>(u.size());

    for (int i = 0; i < N; ++i) {
        double x = x_start + i * dx;
        double u_analytical = analytical_solution(x, time);
        double error = std::abs(u[i] - u_analytical);
        file << x << "\t" << u_analytical << "\t" << u[i] << "\t" << error << "\n";
    }
}

void save_max_error(const std::string& filename,
                   const std::vector<double>& max_errors,
                   double dt) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << "# time\tmax_error\n";
    for (size_t i = 0; i < max_errors.size(); ++i) {
        file << i * dt << "\t" << max_errors[i] << "\n";
    }
}

double calculate_max_error(const std::vector<double>& u,
                         double dx,
                         double time) {
    double max_error = 0.0;
    double x_start = -params::a;

    for (size_t i = 0; i < u.size(); ++i) {
        double x = x_start + i * dx;
        double error = std::abs(u[i] - analytical_solution(x, time));
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}