#include "analytical.h"
#include "params.h"
#include "explicit.h"
#include "implicit.h"
#include "output.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>

void print_progress_header() {
    std::cout << "=============================================\n";
    std::cout << "       CONVERGENCE TEST - DIFFUSION EQ\n";
    std::cout << "=============================================\n";
    std::cout << "Parameters:\n";
    std::cout << "  D = " << params::D << " (diffusion coefficient)\n";
    std::cout << "  t_max = " << params::t_max << " (final time)\n";
    std::cout << "  lambda_explicit = " << params::lambda_explicit << "\n";
    std::cout << "  lambda_implicit = " << params::lambda_implicit << "\n";
    std::cout << "  a = " << params::calculate_a() << " (domain size)\n";
    std::cout << "=============================================\n\n";
}

void run_convergence_test() {
    std::vector<double> h_values = {0.2, 0.1, 0.05, 0.025, 0.0125};
    std::vector<double> explicit_errors, thomas_errors, lu_errors;

    print_progress_header();
    std::cout << "Starting convergence test for " << h_values.size() << " grid sizes...\n";

    for (size_t i = 0; i < h_values.size(); ++i) {
        double h = h_values[i];
        const double a = params::calculate_a();

        // --- ZMIANA: Oblicz Nx tak, aby x=0 zawsze było węzłem siatki ---
        const int Nx_half = static_cast<int>(std::ceil(a / h)); // Liczba punktów od x=0 do x=a
        const int Nx = 2 * Nx_half + 1;                        // Całkowita liczba punktów (nieparzysta)
        const double dx = (2.0 * a) / (Nx - 1);                 // Rzeczywisty krok siatki (może się różnić od h!)

        std::cout << "\n=== TEST CASE " << i+1 << "/" << h_values.size() << " ===\n";
        std::cout << "Requested h = " << h << ", Actual dx = " << dx << ", Nx = " << Nx << "\n";

        // Initialize solution vectors
        std::vector<double> u_explicit(Nx), u_thomas(Nx), u_lu(Nx);
        double x_start = -a;
        for (int j = 0; j < Nx; ++j) {
            double x = x_start + j * dx;
            if      (x <  0.0) u_explicit[j] = u_thomas[j] = u_lu[j] = 1.0;
            else if (x >  0.0) u_explicit[j] = u_thomas[j] = u_lu[j] = 0.0;
            else if (x = 0.0) u_explicit[j] = u_thomas[j] = u_lu[j] = 0.5;
        }

        // Calculate time steps (używamy rzeczywistego dx zamiast h!)
        const double dt_explicit = params::lambda_explicit * dx * dx / params::D;
        const double dt_implicit = params::lambda_implicit * dx * dx / params::D;
        const int Nt_explicit = static_cast<int>(params::t_max / dt_explicit) + 1;
        const int Nt_implicit = static_cast<int>(params::t_max / dt_implicit) + 1;

        std::cout << "Time steps - Explicit: " << Nt_explicit << " (dt = " << dt_explicit << ")\n";
        std::cout << "Time steps - Implicit: " << Nt_implicit << " (dt = " << dt_implicit << ")\n";

        // Solve with different methods
        auto start_time = std::chrono::high_resolution_clock::now();

        std::cout << "\nRunning EXPLICIT SCHEME...\n";
        solve_explicit(u_explicit, dx, dt_explicit, Nt_explicit, params::D);
        auto explicit_time = std::chrono::high_resolution_clock::now();
        std::cout << "Completed in "
                 << std::chrono::duration_cast<std::chrono::milliseconds>(explicit_time - start_time).count()
                 << " ms\n";

        std::cout << "\nRunning THOMAS ALGORITHM...\n";
        solve_implicit_thomas(u_thomas, dx, dt_implicit, Nt_implicit, params::D);
        auto thomas_time = std::chrono::high_resolution_clock::now();
        std::cout << "Completed in "
                 << std::chrono::duration_cast<std::chrono::milliseconds>(thomas_time - explicit_time).count()
                 << " ms\n";

        std::cout << "\nRunning LU DECOMPOSITION...\n";
        solve_implicit_lu(u_lu, dx, dt_implicit, Nt_implicit, params::D);
        auto lu_time = std::chrono::high_resolution_clock::now();
        std::cout << "Completed in "
                 << std::chrono::duration_cast<std::chrono::milliseconds>(lu_time - thomas_time).count()
                 << " ms\n";

        // Calculate and store errors
        explicit_errors.push_back(calculate_max_error(u_explicit, dx, params::t_max));
        thomas_errors.push_back(calculate_max_error(u_thomas, dx, params::t_max));
        lu_errors.push_back(calculate_max_error(u_lu, dx, params::t_max));

        std::cout << "\nCurrent errors:\n";
        std::cout << "  Explicit: " << explicit_errors.back() << "\n";
        std::cout << "  Thomas:   " << thomas_errors.back() << "\n";
        std::cout << "  LU:       " << lu_errors.back() << "\n";
    }

    // Save results for plotting
    save_log_error_data("error_explicit.dat", h_values, explicit_errors);
    save_log_error_data("error_thomas.dat", h_values, thomas_errors);
    save_log_error_data("error_lu.dat", h_values, lu_errors);

    std::cout << "\n=============================================\n";
    std::cout << " CONVERGENCE TEST COMPLETED SUCCESSFULLY\n";
    std::cout << "Results saved to:\n";
    std::cout << "- error_explicit.dat\n";
    std::cout << "- error_thomas.dat\n";
    std::cout << "- error_lu.dat\n";
    std::cout << "=============================================\n";
}

int main() {
    try {
        auto total_start = std::chrono::high_resolution_clock::now();
        run_convergence_test();
        auto total_end = std::chrono::high_resolution_clock::now();

        std::cout << "\nTotal execution time: "
                 << std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start).count()
                 << " seconds\n";
    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}