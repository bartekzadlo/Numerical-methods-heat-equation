#pragma once
#include <cmath>

constexpr double D = 1.0;
constexpr double t_max = 2.0;
constexpr double a = 6.0 * std::sqrt(D * t_max); // a >= 6*sqrt(D*t_max)

constexpr double dx = 0.05;  // krok przestrzenny
constexpr double lambda_explicit = 0.4; // stabilny dla metody eksplicytnej
constexpr double dt_explicit = lambda_explicit * dx * dx / D;

constexpr double lambda_implicit = 1.0; // dla Crank-Nicolson
constexpr double dt_implicit = lambda_implicit * dx * dx / D;

constexpr int Nx = static_cast<int>(2 * a / dx) + 1;  // liczba punktów przestrzennych
constexpr int Nt_explicit = static_cast<int>(t_max / dt_explicit) + 1;
constexpr int Nt_implicit = static_cast<int>(t_max / dt_implicit) + 1;
