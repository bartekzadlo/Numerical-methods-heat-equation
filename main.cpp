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
#include <fstream>

void save_time_error_data(const std::string& filename,
                         const std::vector<double>& time_points,
                         const std::vector<double>& errors) {
    std::ofstream file(filename);
    if (!file) throw std::runtime_error("Cannot open file: " + filename);

    file << "# t max_error\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        file << time_points[i] << " " << errors[i] << "\n";
    }
}

void run_time_error_analysis() {
    const double D = params::D;
    const double t_max = params::t_max;
    const double dx = 0.0249567;
    const int Nx = 681;
    const double dt_explicit = 0.000249135;
    const double dt_implicit = 0.000622837;
    const double a = params::calculate_a();

    // Punkty czasowe co 0.1s
    const int num_time_points = static_cast<int>(t_max / 0.1) + 1;
    std::vector<double> time_points(num_time_points);
    for (int i = 0; i < num_time_points; ++i) {
        time_points[i] = i * 0.1;
    }

    std::vector<double> explicit_errors(num_time_points);
    std::vector<double> thomas_errors(num_time_points);
    std::vector<double> lu_errors(num_time_points);

    // Rozwiązanie początkowe
    std::vector<double> u_explicit(Nx), u_thomas(Nx), u_lu(Nx);
    double x_start = -a;
    for (int j = 0; j < Nx; ++j) {
        double x = x_start + j * dx;
        double u0 = (x < 0.0) ? 1.0 : (x > 0.0 ? 0.0 : 0.5);
        u_explicit[j] = u_thomas[j] = u_lu[j] = u0;
    }

    // Kopie robocze
    std::vector<double> u_explicit_copy = u_explicit;
    std::vector<double> u_thomas_copy = u_thomas;
    std::vector<double> u_lu_copy = u_lu;

    std::cout << "Rozpoczynam analizę błędów w czasie...\n";

    // === Metoda jawna ===
    std::cout << "\n[1] Metoda jawna (klasyczna)...\n";
    for (int i = 0; i < num_time_points; ++i) {
        double t_target = time_points[i];
        int steps = static_cast<int>(t_target / dt_explicit + 1e-9);
        double actual_t = steps * dt_explicit;

        solve_explicit(u_explicit_copy, dx, dt_explicit, steps, D);
        explicit_errors[i] = calculate_max_error(u_explicit_copy, dx, actual_t);

        std::cout << "t = " << std::fixed << std::setprecision(2) << actual_t
                  << ", max error = " << explicit_errors[i] << "\n";
    }

    // === Crank-Nicolson: Thomas ===
    std::cout << "\n[2] Crank-Nicolson (Thomas)...\n";
    for (int i = 0; i < num_time_points; ++i) {
        double t_target = time_points[i];
        int steps = static_cast<int>(t_target / dt_implicit + 1e-9);
        double actual_t = steps * dt_implicit;

        solve_implicit_thomas(u_thomas_copy, dx, dt_implicit, steps, D);
        thomas_errors[i] = calculate_max_error(u_thomas_copy, dx, actual_t);

        std::cout << "t = " << std::fixed << std::setprecision(2) << actual_t
                  << ", max error = " << thomas_errors[i] << "\n";
    }

    // === Crank-Nicolson: LU ===
    std::cout << "\n[3] Crank-Nicolson (LU)...\n";
    for (int i = 0; i < num_time_points; ++i) {
        double t_target = time_points[i];
        int steps = static_cast<int>(t_target / dt_implicit + 1e-9);
        double actual_t = steps * dt_implicit;

        solve_implicit_lu(u_lu_copy, dx, dt_implicit, steps, D);
        lu_errors[i] = calculate_max_error(u_lu_copy, dx, actual_t);

        std::cout << "t = " << std::fixed << std::setprecision(2) << actual_t
                  << ", max error = " << lu_errors[i] << "\n";
    }

    // Zapis danych do plików
    save_time_error_data("time_error_explicit.dat", time_points, explicit_errors);
    save_time_error_data("time_error_thomas.dat", time_points, thomas_errors);
    save_time_error_data("time_error_lu.dat", time_points, lu_errors);

    std::cout << "\nAnaliza zakończona. Dane zapisane do plików:\n";
    std::cout << " - time_error_explicit.dat\n";
    std::cout << " - time_error_thomas.dat\n";
    std::cout << " - time_error_lu.dat\n";
}


int main() {
    try {
        auto total_start = std::chrono::high_resolution_clock::now();
        run_time_error_analysis();
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