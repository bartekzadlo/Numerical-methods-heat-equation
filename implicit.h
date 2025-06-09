#ifndef IMPLICIT_H
#define IMPLICIT_H

#include <vector>

void solve_implicit_thomas(std::vector<double>& u, double dx, double dt, int Nt, double D);
void solve_implicit_lu(std::vector<double>& u, double dx, double dt, int Nt, double D);

#endif