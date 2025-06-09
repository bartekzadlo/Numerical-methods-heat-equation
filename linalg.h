#ifndef LINALG_H
#define LINALG_H

#include <vector>

// Rozwiązuje układ równań z macierzą trójprzekątną metodą Thomasa
void thomas_solver(const std::vector<double>& a,
                   const std::vector<double>& b,
                   const std::vector<double>& c,
                   std::vector<double>& d);

// Rozkład LU macierzy kwadratowej A na L (dolną) i U (górną)
void lu_decompose(const std::vector<std::vector<double>>& A,
                  std::vector<std::vector<double>>& L,
                  std::vector<std::vector<double>>& U);

// Rozwiązuje układ równań LUx = b, gdzie L i U to macierze z rozkładu LU
std::vector<double> lu_solve(const std::vector<std::vector<double>>& L,
                             const std::vector<std::vector<double>>& U,
                             const std::vector<double>& b);

#endif // LINALG_H
