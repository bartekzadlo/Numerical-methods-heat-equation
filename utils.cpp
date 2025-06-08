// utils.cpp
#include <vector>
#include "params.h"

// Inicjalizacja siatki przestrzennej
std::vector<double> init_space_grid() {
    std::vector<double> x(Nx);
    for (int i = 0; i < Nx; ++i) {
        x[i] = -a + i * h;
    }
    return x;
}

// Inicjalizacja warunku początkowego
std::vector<double> init_condition(const std::vector<double>& x) {
    std::vector<double> U0(Nx);
    for (int i = 0; i < Nx; ++i) {
        U0[i] = (x[i] < 0.0) ? 1.0 : 0.0;
    }
    return U0;
}