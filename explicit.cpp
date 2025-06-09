#include "explicit.h"
#include "params.h"
#include <stdexcept>  // Do obsługi wyjątków

void solve_explicit(std::vector<double>& u, double dx, double dt, int Nt) {
    const double lambda = params::D * dt / (dx * dx);

    // Zabezpieczenie przed niestabilnością
    if (lambda > 0.5) {
        throw std::invalid_argument(
            "Explicit scheme requires lambda <= 0.5 (got lambda = " +
            std::to_string(lambda) + "). Adjust dt or dx."
        );
    }

    const int N = static_cast<int>(u.size());
    if (N < 2) {
        throw std::invalid_argument("Vector u must have at least 2 elements.");
    }

    // Podwójne buforowanie (u_prev -> aktualny stan, u_next -> nowy stan)
    std::vector<double> u_prev = u;  // Kopia danych wejściowych
    std::vector<double> u_next(N);   // Bufor na nowe wartości

    for (int n = 0; n < Nt; ++n) {
        // Warunki brzegowe (stałe Dirichleta)
        u_next[0] = 1.0;       // U(-a, t) = 1
        u_next[N - 1] = 0.0;   // U(+a, t) = 0

        // Główna pętla obliczeniowa
        for (int i = 1; i < N - 1; ++i) {
            u_next[i] = u_prev[i] + lambda * (u_prev[i + 1] - 2 * u_prev[i] + u_prev[i - 1]);
        }

        std::swap(u_prev, u_next);  // Zamiana buforów (bez kopiowania)
    }

    u = u_prev;  // Przekazanie wyniku do wejściowego wektora
}