#include "explicit.h"
#include "params.h"

#include <vector>

// Metoda eksplicytna

void solve_explicit(std::vector<double>& u, double dx, double dt, int Nt) {
    const double lambda = params::D * dt / (dx * dx);
    int N = static_cast<int>(u.size());
    std::vector<double> u_new(N);

    for (int n = 0; n < Nt; ++n) {
        // Warunki brzegowe
        u_new[0] = 1.0;
        u_new[N - 1] = 0.0;

        for (int i = 1; i < N - 1; ++i) {
            u_new[i] = u[i] + lambda * (u[i + 1] - 2 * u[i] + u[i - 1]);
        }
        u = u_new;
    }
}
