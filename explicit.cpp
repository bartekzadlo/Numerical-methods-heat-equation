#include "explicit.h"
#include <stdexcept>

void solve_explicit(std::vector<double>& u, double dx, double dt, int Nt, double D) {
    const double lambda = D * dt / (dx * dx);
    const int N = u.size();

    if (lambda > 0.5) {
        throw std::invalid_argument("Warunek stabilności nie spełniony (lambda > 0.5)");
    }

    std::vector<double> u_new(N);

    for (int n = 0; n < Nt; ++n) {
        u_new[0] = 1.0;       // Warunek brzegowy lewy
        u_new[N-1] = 0.0;     // Warunek brzegowy prawy

        for (int i = 1; i < N-1; ++i) {
            u_new[i] = u[i] + lambda * (u[i+1] - 2*u[i] + u[i-1]);
        }

        std::swap(u, u_new);
    }
}