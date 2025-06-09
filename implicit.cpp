#include "implicit.h"
#include "linalg.h"
#include <stdexcept>

void solve_implicit_thomas(std::vector<double>& u, double dx, double dt, int Nt, double D) {
    const int N = u.size();
    const double lambda = D * dt / (dx * dx);

    std::vector<double> a(N-2, -lambda/2.0); // Podprzekątna
    std::vector<double> b(N-2, 1.0 + lambda); // Przekątna
    std::vector<double> c(N-2, -lambda/2.0); // Nadprzekątna
    std::vector<double> d(N-2);              // Prawa strona

    for (int n = 0; n < Nt; ++n) {
        // Ustaw warunki brzegowe
        u[0] = 1.0;    // U(-a, t) = 1
        u[N-1] = 0.0;  // U(a, t) = 0

        for (int i = 1; i < N-1; ++i) {
            double left = u[i-1];
            double right = u[i+1];
            d[i-1] = (1.0 - lambda)*u[i] + (lambda/2.0)*(left + right);
        }

        thomas_solver(a, b, c, d);

        for (int i = 1; i < N-1; ++i) {
            u[i] = d[i-1];
        }
    }
}

void solve_implicit_lu(std::vector<double>& u, double dx, double dt, int Nt, double D) {
    const int N = u.size();
    const double lambda = D * dt / (dx * dx);
    const int system_size = N - 2;

    // Przygotowanie macierzy
    std::vector<std::vector<double>> A(system_size, std::vector<double>(system_size, 0.0));
    for (int i = 0; i < system_size; ++i) {
        A[i][i] = 1.0 + lambda;
        if (i > 0) A[i][i-1] = -lambda/2.0;
        if (i < system_size-1) A[i][i+1] = -lambda/2.0;
    }

    // Dekompozycja LU
    std::vector<std::vector<double>> L, U;
    lu_decompose(A, L, U);

    for (int n = 0; n < Nt; ++n) {
        // Ustaw warunki brzegowe
        u[0] = 1.0;    // U(-a, t) = 1
        u[N-1] = 0.0;  // U(a, t) = 0

        std::vector<double> b(system_size);
        for (int i = 0; i < system_size; ++i) {
            double left = u[i];     // i+1-1 = i
            double right = u[i+2];  // i+1+1 = i+2
            b[i] = (1.0 - lambda)*u[i+1] + (lambda/2.0)*(left + right);
        }

        std::vector<double> solution = lu_solve(L, U, b);

        for (int i = 0; i < system_size; ++i) {
            u[i+1] = solution[i];
        }
    }
}