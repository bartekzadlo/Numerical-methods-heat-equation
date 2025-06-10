#include "output.h"
#include "analytical.h"
#include "params.h"
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>

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

void save_solutions_for_time(const std::string& filename,
                             const std::vector<double>& x,
                             const std::vector<double>& u_explicit,
                             const std::vector<double>& u_thomas,
                             const std::vector<double>& u_lu,
                             double t) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Nie można otworzyć pliku do zapisu: " << filename << "\n";
        return;
    }

    file << "# x u_explicit u_thomas u_lu analytical\n";

    int skip = 5; // zapisuj co 5-ty punkt, zmień według potrzeby

    for (size_t i = 0; i < x.size(); i += skip) {
        double u_analytical = analytical_solution(x[i], t);
        file << std::setprecision(10) << x[i] << " "
             << u_explicit[i] << " "
             << u_thomas[i] << " "
             << u_lu[i] << " "
             << u_analytical << "\n";
    }
}
