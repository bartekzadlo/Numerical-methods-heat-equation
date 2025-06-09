#include "implicit.h"
#include "params.h"
#include "linalg.h"

#include <vector>

void solve_implicit(std::vector<double>& u, double dx, double dt, int Nt) {
    int N = static_cast<int>(u.size());
    double lambda = params::D * dt / (dx * dx);

    // macierz trójdiagonalna A (Crank-Nicolson)
    std::vector<double> a(N - 2, -lambda / 2.0); // podprzekątna
    std::vector<double> b(N - 2, 1.0 + lambda);  // przekątna główna
    std::vector<double> c(N - 2, -lambda / 2.0); // nadprzekątna

    std::vector<double> d(N - 2, 0.0); // wektor prawej strony

    for (int n = 0; n < Nt; ++n) {
        for (int i = 1; i < N - 1; ++i) {
            double left = (i == 1) ? u[0] : u[i - 1];
            double right = (i == N - 2) ? u[N - 1] : u[i + 1];
            d[i - 1] = (1.0 - lambda) * u[i] + (lambda / 2.0) * (left + right);
        }

        // Warunki brzegowe
        d[0] += (lambda / 2.0) * 1.0;  // u[0] = 1.0
        d[N - 3] += (lambda / 2.0) * 0.0;  // u[N-1] = 0.0

        // Rozwiązanie układu
        thomas_solver(a, b, c, d);

        for (int i = 1; i < N - 1; ++i) {
            u[i] = d[i - 1];
        }

        u[0] = 1.0;
        u[N - 1] = 0.0;
    }
}
