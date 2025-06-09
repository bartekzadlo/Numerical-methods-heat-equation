// utils.cpp
#include "params.h"

// Inicjalizacja siatki przestrzennej
void init_space_grid(double x[Nx]) {
    for (int i = 0; i < Nx; ++i) {
        x[i] = -a + i * h;
    }
}

// Inicjalizacja warunku początkowego
void init_condition(const double x[Nx], double U0[Nx]) {
    for (int i = 0; i < Nx; ++i) {
        U0[i] = (x[i] < 0.0) ? 1.0 : 0.0;
    }
}

// Inicjalizacja warunku brzegowego
void apply_boundary_conditions(double U[Nx]) {
    U[0] = 1.0;
    U[Nx - 1] = 0.0;
}

