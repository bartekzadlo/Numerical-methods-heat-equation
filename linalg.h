#ifndef LINALG_H
#define LINALG_H

#include <vector>

void thomas_solver(const std::vector<double>& a, const std::vector<double>& b,
                   const std::vector<double>& c, std::vector<double>& d);
void lu_decompose(const std::vector<std::vector<double>>& A,
                  std::vector<std::vector<double>>& L,
                  std::vector<std::vector<double>>& U);
std::vector<double> lu_solve(const std::vector<std::vector<double>>& L,
                             const std::vector<std::vector<double>>& U,
                             const std::vector<double>& b);

#endif