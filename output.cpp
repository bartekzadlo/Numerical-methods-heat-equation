#include "output.h"
#include "analytical.h"
#include <fstream>
#include <cmath>

void save_log_error_data(const std::string& filename,
                        const std::vector<double>& h_values,
                        const std::vector<double>& errors) {
    std::ofstream file(filename);
    if (!file) throw std::runtime_error("Cannot open file: " + filename);

    file << "# log10(h) log10(error)\n";
    for (size_t i = 0; i < h_values.size(); ++i) {
        file << std::log10(h_values[i]) << " " << std::log10(errors[i]) << "\n";
    }
}

double calculate_max_error(const std::vector<double>& u, double dx, double t) {
    double max_error = 0.0;
    double a = 6.0 * std::sqrt(params::D * t);
    double x_start = -a;

    for (size_t i = 0; i < u.size(); ++i) {
        double x = x_start + i * dx;
        double error = std::abs(u[i] - analytical_solution(x, t));
        if (error > max_error) {
            max_error = error;
        }
    }
    return max_error;
}