#pragma once

#include <cmath>

namespace params {
    constexpr double D = 1.0;               // Współczynnik dyfuzji
    constexpr double t_max = 2.0;           // Czas symulacji

    inline double calculate_a() {
        return 6.0 * std::sqrt(D * t_max); // a >= 6 * sqrt(D * t_max)
    }

    // Dla analizy błędów - h będzie zmieniane w main.cpp
    constexpr double lambda_explicit = 0.4; // Stabilny dla metody jawnej
    constexpr double lambda_implicit = 1.0;  // Dla Crank-Nicolson
}