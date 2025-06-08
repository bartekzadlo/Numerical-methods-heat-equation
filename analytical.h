#ifndef ANALYTICAL_H
#define ANALYTICAL_H

#include "params.h"
#include "calerf.h"

/**
 * Oblicza rozwiazanie analityczne w punkcie (x,t)
 * zgodnie ze wzorem: U(x,t) = 1/2 * erfc(x/(2?(Dt)))
 *
 * @param x Polozenie przestrzenne
 * @param t Czas
 * @return Wartosc rozwiazania analitycznego
 */
double analytical_solution(double x, double t);

/**
 * Oblicza rozwiazanie analityczne dla calej siatki przestrzennej
 *
 * @param x_grid Tablica wspolrzednych przestrzennych
 * @param Nx Liczba punktow siatki
 * @param t Czas
 * @param result Tablica wynikowa (musi byc zaalokowana przed wywo�aniem)
 */
void calculate_analytical_solution(const double* x_grid, int Nx, double t, double* result);

/**
 * Oblicza bledy bezwzgledne miedzy rozwiazaniami numerycznymi i analitycznymi
 *
 * @param numerical Tablica rozwiazan numerycznych
 * @param analytical Tablica rozwiazan analitycznych
 * @param N Liczba punktow
 * @param errors Tablica bledow (musi byc zaalokowana przed wywolaniem)
 */
void calculate_absolute_errors(const double* numerical, const double* analytical, int N, double* errors);

#endif // ANALYTICAL_H