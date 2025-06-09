#include "params.h"
#include "analytical.h"
#include "explicit.h"
#include "implicit.h"
#include "output.h"

#include <iostream>
#include <vector>

int main() {
    using namespace params;

    std::vector<double> u_explicit(Nx);
    std::vector<double> u_implicit(Nx);

    double x_start = -a;

    // Warunek początkowy: 1 dla x<0, 0 dla x>=0
    for (int i = 0; i < Nx; ++i) {
        double x = x_start + i * dx;
        u_explicit[i] = (x < 0.0) ? 1.0 : 0.0;
        u_implicit[i] = u_explicit[i];
    }

    std::cout << "Solving explicit method...\n";
    solve_explicit(u_explicit, dx, dt_explicit, Nt_explicit);

    std::cout << "Solving Crank-Nicolson method...\n";
    solve_implicit(u_implicit, dx, dt_implicit, Nt_implicit);

    save_profile("profile_explicit.dat", u_explicit, dx, t_max);
    save_profile("profile_implicit.dat", u_implicit, dx, t_max);

    std::cout << "Results saved to profile_explicit.dat and profile_implicit.dat\n";

    return 0;
}
