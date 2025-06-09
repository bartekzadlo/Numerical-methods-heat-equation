#include "params.h"
#include "analytical.h"
#include "explicit.h"
#include "implicit.h"
#include "output.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>

int main() {
    using namespace params;

    // Inicjalizacja wektorów
    std::vector<double> u_explicit(Nx);
    std::vector<double> u_implicit(Nx);
    double x_start = -a;

    // Warunek początkowy: U(x,0) = 1 dla x < 0, 0 dla x >= 0
    for (int i = 0; i < Nx; ++i) {
        double x = x_start + i * dx;
        u_explicit[i] = (x < 0.0) ? 1.0 : 0.0;
        u_implicit[i] = u_explicit[i];
    }

    try {
        // Rozwiązanie metodą jawną
        std::cout << "Solving explicit method (lambda = "
                  << params::D * dt_explicit / (dx * dx) << ")...\n";
        solve_explicit(u_explicit, dx, dt_explicit, Nt_explicit);

        // Rozwiązanie metodą Crank-Nicolson
        std::cout << "Solving Crank-Nicolson method (lambda = "
                  << params::D * dt_implicit / (dx * dx) << ")...\n";
        solve_implicit(u_implicit, dx, dt_implicit, Nt_implicit);

        // Zapis wyników
        save_profile("profile_explicit.dat", u_explicit, dx, t_max);
        save_profile("profile_implicit.dat", u_implicit, dx, t_max);

        // Obliczenie i zapis błędów
        std::vector<double> max_errors = {
            calculate_max_error(u_explicit, dx, t_max),
            calculate_max_error(u_implicit, dx, t_max)
        };
        save_max_error("max_errors.dat", max_errors, t_max);

        std::cout << "\nResults saved to:\n"
                  << "- profile_explicit.dat (explicit solution profile)\n"
                  << "- profile_implicit.dat (implicit solution profile)\n"
                  << "- max_errors.dat (max errors at t_max)\n\n"
                  << "Max errors at t = " << t_max << ":\n"
                  << "- Explicit: " << max_errors[0] << "\n"
                  << "- Implicit: " << max_errors[1] << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}