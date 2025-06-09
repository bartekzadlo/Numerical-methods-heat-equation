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

    std::vector<double> dx_values = {0.2, 0.1, 0.05, 0.025, 0.0125};
    std::vector<double> log_dx, log_err_explicit, log_err_thomas, log_err_lu;

    for (double dx_local : dx_values) {
        int Nx_local = static_cast<int>(2.0 * a / dx_local) + 1;
        double dt_local = 0.4 * dx_local * dx_local / D; // lambda < 0.5
        int Nt_local = static_cast<int>(t_max / dt_local);

        std::vector<double> u_init(Nx_local), u_exp(Nx_local), u_thomas(Nx_local), u_lu(Nx_local);

        // Warunek początkowy
        for (int i = 0; i < Nx_local; ++i) {
            double x = -a + i * dx_local;
            u_init[i] = (x < 0.0) ? 1.0 : 0.0;
        }

        u_exp = u_thomas = u_lu = u_init;

        std::cout << "dx = " << dx_local << ", dt = " << dt_local << ", Nx = " << Nx_local << ", Nt = " << Nt_local << "\n";

        // Rozwiązania
        solve_explicit(u_exp, dx_local, dt_local, Nt_local);
        solve_implicit(u_thomas, dx_local, dt_local, Nt_local);        // Thomas
        solve_implicit_lu(u_lu, dx_local, dt_local, Nt_local);         // LU

        // Rozwiązanie analityczne
        std::vector<double> u_exact = analytical_solution(Nx_local, dx_local, t_max);

        // Błędy
        double err_exp = calculate_max_error(u_exp, u_exact);
        double err_thomas = calculate_max_error(u_thomas, u_exact);
        double err_lu = calculate_max_error(u_lu, u_exact);

        log_dx.push_back(std::log10(dx_local));
        log_err_explicit.push_back(std::log10(err_exp));
        log_err_thomas.push_back(std::log10(err_thomas));
        log_err_lu.push_back(std::log10(err_lu));
    }

    // Zapis
    save_log_error_plot("error_explicit.dat", log_dx, log_err_explicit);
    save_log_error_plot("error_thomas.dat", log_dx, log_err_thomas);
    save_log_error_plot("error_lu.dat", log_dx, log_err_lu);

    std::cout << "Saved:\n- error_explicit.dat\n- error_thomas.dat\n- error_lu.dat\n";
    return 0;
}