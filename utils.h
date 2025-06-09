#ifndef UTILS_H
#define UTILS_H

// Zakładamy, że Nx jest zdefiniowane w params.h
#include "params.h"

// Funkcje inicjalizujące (zamiast zwracać wektory, modyfikują przekazane tablice)
void init_space_grid(double x[Nx]);
void init_condition(const double x[Nx], double U0[Nx]);

#endif // UTILS_H
