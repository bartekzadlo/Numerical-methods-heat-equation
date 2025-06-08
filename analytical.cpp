#include "analytical.h"
#include <cmath>
#include <stdexcept>

double analytical_solution(double x, double t) {
    // Zgodnie z poleceniem: U(x,0) = 1 dla x < 0, 0 dla x >= 0
    if (t <= 0.0) {
        return (x < 0.0) ? 1.0 : 0.0;
    }

    // Zgodnie ze wzorem: z = x/(2?(Dt))
    // gdzie D - wspolczynnik dyfuzji (zdefiniowany w params.h)
    const double denominator = 2.0 * std::sqrt(D * t);
    const double z = x / denominator;

    // Zgodnie ze wzorem: U(x,t) = 1/2 * erfc(z)
    return 0.5 * calerfpack::erfc_LD(z);
}

void calculate_analytical_solution(const double* x_grid, int Nx, double t, double* result) {
    // Krok 1: Sprawdzenie poprawnosci wskaznikow
    if (x_grid == nullptr || result == nullptr) {
        throw std::invalid_argument("Invalid pointer passed");
    }

    // Krok 2: Obliczenie rozwiazania dla kazdego punktu siatki
    for (int i = 0; i < Nx; ++i) {
        result[i] = analytical_solution(x_grid[i], t);
    }
}

void calculate_absolute_errors(const double* numerical, const double* analytical, int N, double* errors) {
    // Krok 1: Sprawdzenie poprawnosci wskaznikow
    if (numerical == nullptr || analytical == nullptr || errors == nullptr) {
        throw std::invalid_argument("Invalid pointer passed");
    }

    // Krok 2: Obliczenie błędów bezwzględnych
    for (int i = 0; i < N; ++i) {
        errors[i] = std::abs(numerical[i] - analytical[i]);
    }
}