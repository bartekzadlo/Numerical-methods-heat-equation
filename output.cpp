#include "output.h"
#include "analytical.h"
#include "params.h"

#include <fstream>
#include <cmath>

void save_profile(const std::string& filename, const std::vector<double>& u,
                  double dx, double time) {
    std::ofstream file(filename);
    int N = static_cast<int>(u.size());
    double x_start = -params::a;

    for (int i = 0; i < N; ++i) {
        double x = x_start + i * dx;
        double u_analytical = analytical_solution(x, time);
        double error = std::abs(u[i] - u_analytical);
        file << x << "\t" << u_analytical << "\t" << u[i] << "\t" << error << "\n";
    }
}

void save_error(const std::string& filename,
                const std::vector<double>& errors,
                double dt) {
    std::ofstream file(filename);
    for (size_t i = 0; i < errors.size(); ++i) {
        file << i * dt << "\t" << errors[i] << "\n";
    }
}
