#include "implicit.h"
#include "params.h"
#include "linalg.h"
#include <stdexcept>

void solve_implicit(std::vector<double>& u, double dx, double dt, int Nt) {
    int N = static_cast<int>(u.size());
    if (N < 3) throw std::invalid_argument("Vector size must be >= 3 for implicit scheme.");

    double lambda = params::D * dt / (dx * dx);

    // Inicjalizacja stałych macierzy (przeniesione przed pętlę)
    std::vector<double> a(N - 2, -lambda / 2.0); // podprzekątna
    std::vector<double> b(N - 2, 1.0 + lambda);  // przekątna
    std::vector<double> c(N - 2, -lambda / 2.0); // nadprzekątna
    std::vector<double> d(N - 2);                // prawa strona

    for (int n = 0; n < Nt; ++n) {
        // Wypełnij wektor prawej strony
        for (int i = 1; i < N - 1; ++i) {
            double left = (i == 1) ? u[0] : u[i - 1];
            double right = (i == N - 2) ? u[N - 1] : u[i + 1];
            d[i - 1] = (1.0 - lambda) * u[i] + (lambda / 2.0) * (left + right);
        }

        // Uwzględnij warunki brzegowe w prawiej stronie
        d[0] += (lambda / 2.0) * 1.0;    // Lewy brzeg: U[0] = 1.0
        d[N - 3] += (lambda / 2.0) * 0.0; // Prawy brzeg: U[N-1] = 0.0

        // Rozwiąż układ równań
        thomas_solver(a, b, c, d);

        // Zaktualizuj wektor u (pomiń brzegi)
        for (int i = 1; i < N - 1; ++i) {
            u[i] = d[i - 1];
        }

        // Ustaw warunki brzegowe explicite
        u[0] = 1.0;
        u[N - 1] = 0.0;
    }
}