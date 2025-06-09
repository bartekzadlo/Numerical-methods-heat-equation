#include "params.h"
#include "analytical.h"
#include "explicit.h"
#include "implicit.h"
#include "output.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

void run_convergence_test() {
    std::vector<double> h_values = {0.1, 0.05, 0.025, 0.0125, 0.00625};
    std::vector<double> explicit_errors, thomas_errors, lu_errors;

    for (double h : h_values) {
        const double a = params::calculate_a();
        const int Nx = static_cast<int>(2 * a / h) + 1;
        const double dx = h;

        // Initialize solution vectors
        std::vector<double> u_explicit(Nx), u_thomas(Nx), u_lu(Nx);
        double x_start = -a;
        for (int i = 0; i < Nx; ++i) {
            double x = x_start + i * dx;
            double val = (x < 0.0) ? 1.0 : 0.0;
            u_explicit[i] = u_thomas[i] = u_lu[i] = val;
        }

        // Calculate time steps
        const double dt_explicit = params::lambda_explicit * dx * dx / params::D;
        const double dt_implicit = params::lambda_implicit * dx * dx / params::D;
        const int Nt_explicit = static_cast<int>(params::t_max / dt_explicit) + 1;
        const int Nt_implicit = static_cast<int>(params::t_max / dt_implicit) + 1;

        // Solve with different methods
        solve_explicit(u_explicit, dx, dt_explicit, Nt_explicit, params::D);
        solve_implicit_thomas(u_thomas, dx, dt_implicit, Nt_implicit, params::D);
        solve_implicit_lu(u_lu, dx, dt_implicit, Nt_implicit, params::D);

        // Calculate and store errors
        explicit_errors.push_back(calculate_max_error(u_explicit, dx, params::t_max));
        thomas_errors.push_back(calculate_max_error(u_thomas, dx, params::t_max));
        lu_errors.push_back(calculate_max_error(u_lu, dx, params::t_max));
    }

    // Save results for plotting
    save_log_error_data("error_explicit.dat", h_values, explicit_errors);
    save_log_error_data("error_thomas.dat", h_values, thomas_errors);
    save_log_error_data("error_lu.dat", h_values, lu_errors);

    std::cout << "Convergence test completed. Results saved to:\n"
              << "- error_explicit.dat\n"
              << "- error_thomas.dat\n"
              << "- error_lu.dat\n";
}

int main() {
    try {
        run_convergence_test();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}