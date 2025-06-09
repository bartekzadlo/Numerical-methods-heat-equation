#include "linalg.h"
#include <stdexcept>

// Metoda Thomasa - rozwiązuje układ równań z macierzą trójprzekątną
void thomas_solver(const std::vector<double>& a,
                   const std::vector<double>& b,
                   const std::vector<double>& c,
                   std::vector<double>& d) {
    int n = static_cast<int>(b.size());
    std::vector<double> c_prime(n, 0.0);

    // Etap wstępny (forward sweep) - modyfikacja współczynników
    c_prime[0] = c[0] / b[0];
    d[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double m = b[i] - a[i - 1] * c_prime[i - 1];
        if (i < n - 1)
            c_prime[i] = c[i] / m;
        d[i] = (d[i] - a[i - 1] * d[i - 1]) / m;
    }

    // Etap wsteczny (back substitution) - obliczenie rozwiązania
    for (int i = n - 2; i >= 0; --i) {
        d[i] = d[i] - c_prime[i] * d[i + 1];
    }
}

// Rozkład LU macierzy kwadratowej A na L (dolną) i U (górną)
void lu_decompose(const std::vector<std::vector<double>>& A,
                  std::vector<std::vector<double>>& L,
                  std::vector<std::vector<double>>& U) {
    int n = static_cast<int>(A.size());
    L.assign(n, std::vector<double>(n, 0.0));
    U.assign(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        // Obliczanie macierzy górnej U
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
                sum += L[i][j] * U[j][k];
            U[i][k] = A[i][k] - sum;
        }

        // Obliczanie macierzy dolnej L
        L[i][i] = 1.0; // Jedynki na diagonali L
        for (int k = i + 1; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
                sum += L[k][j] * U[j][i];
            if (U[i][i] == 0.0)
                throw std::runtime_error("LU decomposition: dzielenie przez zero (zerowy pivot)");
            L[k][i] = (A[k][i] - sum) / U[i][i];
        }
    }
}

// Rozwiązanie układu LUx = b, gdzie L i U to macierze z rozkładu LU
std::vector<double> lu_solve(const std::vector<std::vector<double>>& L,
                             const std::vector<std::vector<double>>& U,
                             const std::vector<double>& b) {
    int n = static_cast<int>(b.size());
    std::vector<double> y(n, 0.0);
    std::vector<double> x(n, 0.0);

    // Podstawianie do przodu: Ly = b
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
            sum += L[i][j] * y[j];
        y[i] = b[i] - sum;
    }

    // Podstawianie do tyłu: Ux = y
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j)
            sum += U[i][j] * x[j];
        if (U[i][i] == 0.0)
            throw std::runtime_error("LU solve: dzielenie przez zero (zerowy pivot)");
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}
