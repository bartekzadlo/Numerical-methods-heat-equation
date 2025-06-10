#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "params.h"
#include "analytical.h"
#include "explicit.h"
#include "implicit.h"
#include "output.h"

void generate_solution_plots();

int main() {
    using namespace params;

    // Parametry siatki
    double a = calculate_a();
    double dx = 0.0249567;
    int Nx_half = static_cast<int>(std::ceil(a / dx));
    int Nx = 2 * Nx_half + 1;
    std::vector<double> x(Nx);
    double x_start = -a;

    for (int i = 0; i < Nx; ++i) {
        x[i] = x_start + i * dx;
    }

    // Parametry czasowe
    double dt_explicit = 0.000249135;
    double dt_implicit = 0.000622837;

    // Wygenerowanie wykresów rozwiązań dla różnych czasów
    generate_solution_plots();

    return 0;
}

void generate_solution_plots() {
    using namespace params;

    double a = calculate_a();
    double dx = 0.0249567;
    int Nx_half = static_cast<int>(std::ceil(a / dx));
    int Nx = 2 * Nx_half + 1;
    double x_start = -a;

    std::vector<double> x(Nx);
    for (int i = 0; i < Nx; ++i) {
        x[i] = x_start + i * dx;
    }

    const double dt_explicit = 0.000249135;
    const double dt_implicit = 0.000622837;

    std::vector<double> u_explicit(Nx);
    std::vector<double> u_thomas(Nx);
    std::vector<double> u_lu(Nx);

    // Warunki początkowe
    for (int i = 0; i < Nx; ++i) {
        if (x[i] < 0.0) u_explicit[i] = u_thomas[i] = u_lu[i] = 1.0;
        else if (x[i] > 0.0) u_explicit[i] = u_thomas[i] = u_lu[i] = 0.0;
        else u_explicit[i] = u_thomas[i] = u_lu[i] = 0.5;
    }

    std::vector<double> times = {0.0, 0.25 * t_max, 0.5 * t_max, 0.75 * t_max, t_max};

    for (double t : times) {
        int Nt_explicit = static_cast<int>(t / dt_explicit);
        int Nt_implicit = static_cast<int>(t / dt_implicit);

        std::vector<double> u_exp = u_explicit;
        std::vector<double> u_th = u_thomas;
        std::vector<double> u_l = u_lu;

        if (Nt_explicit > 0) solve_explicit(u_exp, dx, dt_explicit, Nt_explicit, D);
        if (Nt_implicit > 0) {
            solve_implicit_thomas(u_th, dx, dt_implicit, Nt_implicit, D);
            solve_implicit_lu(u_l, dx, dt_implicit, Nt_implicit, D);
        }

        std::string filename = "solution_t_" + std::to_string(static_cast<int>(t * 1000)) + ".dat";
        save_solutions_for_time(filename, x, u_exp, u_th, u_l, t);
        std::cout << "Saved solutions for t = " << t << " to " << filename << "\n";
    }
}
