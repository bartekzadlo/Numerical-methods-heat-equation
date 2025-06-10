#include "implicit.h"
#include "linalg.h"
#include "analytical.h"
#include "params.h"
#include <stdexcept>

void solve_implicit_thomas(std::vector<double>& u, double dx, double dt, int Nt, double D) {
    const int N = u.size();
    const double lambda = D * dt / (dx * dx);

    std::vector<double> a(N-2, -lambda/2.0);
    std::vector<double> b(N-2, 1.0 + lambda);
    std::vector<double> c(N-2, -lambda/2.0);
    std::vector<double> d(N-2);

    for (int n = 0; n < Nt; ++n) {
        // ZAPISZ WARUNKI BRZEGOWE Z POPRZEDNIEGO KROKU
        const double left_bc = u[0];     // Zachowaj u(0, t_n)
        const double right_bc = u[N-1];  // Zachowaj u(a, t_n)

        // OBLICZ PRAWĄ STRONĘ (użyj starych wartości)
        for (int i = 1; i < N-1; ++i) {
            double left = u[i-1];
            double right = u[i+1];
            d[i-1] = (1.0 - lambda)*u[i] + (lambda/2.0)*(left + right);
        }

        // DODAJ WPŁYW WARUNKÓW BRZEGOWYCH NA RHS
        d[0] += (lambda/2.0) * left_bc;       // Dla pierwszego równania
        d[N-3] += (lambda/2.0) * right_bc;    // Dla ostatniego równania

        // ROZWIĄŻ UKŁAD RÓWNAŃ
        thomas_solver(a, b, c, d);

        // AKTUALIZUJ ROZWIĄZANIE
        for (int i = 1; i < N-1; ++i) {
            u[i] = d[i-1];
        }

        // USTAW NOWE WARUNKI BRZEGOWE (są znane dla t_{n+1})
        u[0] = 1.0;    // U(-a, t_{n+1}) = 1
        u[N-1] = 0.0;  // U(a, t_{n+1}) = 0
    }
}

void solve_implicit_lu(std::vector<double>& u, double dx, double dt, int Nt, double D) {
    const int N = u.size();
    const double lambda = D * dt / (dx * dx);
    const int system_size = N - 2;
    const double a = params::calculate_a();
    const double x_left = -a;
    const double x_right = a;

    // Budowa macierzy A (dla punktów wewnętrznych)
    std::vector<std::vector<double>> A(system_size, std::vector<double>(system_size, 0.0));
    for (int i = 0; i < system_size; ++i) {
        A[i][i] = 1.0 + lambda;
        if (i > 0) A[i][i-1] = -lambda/2.0;
        if (i < system_size-1) A[i][i+1] = -lambda/2.0;
    }

    // Dekompozycja LU (tylko raz)
    std::vector<std::vector<double>> L, U;
    lu_decompose(A, L, U);

    // Prawa strona układu równań
    std::vector<double> b(system_size);

    for (int n = 0; n < Nt; ++n) {
        const double t_old = n * dt;
        const double t_new = (n + 1) * dt;

        // Oblicz nowe warunki brzegowe (dla t_{n+1})
        const double left_bc_new = analytical_solution(x_left, t_new);
        const double right_bc_new = analytical_solution(x_right, t_new);

        // Ustaw bieżące warunki brzegowe (dla t_n)
        u[0] = analytical_solution(x_left, t_old);
        u[N-1] = analytical_solution(x_right, t_old);

        // Oblicz prawą stronę układu
        for (int i = 0; i < system_size; ++i) {
            double left = u[i];     // u_{i}^{n}
            double center = u[i+1]; // u_{i+1}^{n}
            double right = u[i+2];  // u_{i+2}^{n}

            b[i] = (1.0 - lambda) * center + (lambda/2.0) * (left + right);
        }

        // Dodaj wkład od nowych warunków brzegowych
        b[0] += (lambda/2.0) * left_bc_new;
        if (system_size > 1) {
            b[system_size-1] += (lambda/2.0) * right_bc_new;
        }

        // Rozwiąż układ równań
        std::vector<double> solution = lu_solve(L, U, b);

        // Aktualizuj rozwiązanie dla punktów wewnętrznych
        for (int i = 0; i < system_size; ++i) {
            u[i+1] = solution[i];
        }

        // Ustaw nowe warunki brzegowe
        u[0] = left_bc_new;
        u[N-1] = right_bc_new;
    }
}
