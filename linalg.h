#pragma once

#include <vector>

// Rozwiązuje układ trójdiagonalny Ax = d metodą Thomasa
// a - podprzekątna (n-1 elementów)
// b - przekątna (n elementów)
// c - nadprzekątna (n-1 elementów)
// d - wektor prawej strony, wynik w d (n elementów)

void thomas_solver(const std::vector<double>& a,
                   const std::vector<double>& b,
                   const std::vector<double>& c,
                   std::vector<double>& d);
