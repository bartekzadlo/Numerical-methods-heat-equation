#ifndef PARAMS_H
#define PARAMS_H

#include <cmath>

// Stałe do symulacji
constexpr double D = 1.0;                 // współczynnik dyfuzji
constexpr double tmax = 2.0;              // maksymalny czas symulacji
const double a = 6.0 * std::sqrt(D * tmax); // zakres przestrzenny [-a, a] (≈8.4853)

// Parametry siatki
constexpr int Nx = 170;                   // liczba węzłów przestrzennych
const double h = (2.0 * a) / (Nx - 1);   // krok przestrzenny (≈0.099895)

// Warunki stabilności
constexpr double lambda_explicit = 0.4;   // Dla metody jawnej
const double dt_explicit = lambda_explicit * h * h / D; // krok czasowy jawnej (≈0.003992)

constexpr double lambda_implicit = 1.0;   // Dla metody niejawnej (Crank-Nicolson)
const double dt_implicit = lambda_implicit * h * h / D; // krok czasowy niejawnej (≈0.009979)

// Liczba kroków czasowych
const int Nt_explicit = static_cast<int>(tmax / dt_explicit) + 1; // ≈501
const int Nt_implicit = static_cast<int>(tmax / dt_implicit) + 1; // ≈201

#endif // PARAMS_H