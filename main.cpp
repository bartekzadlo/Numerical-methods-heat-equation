#include <iostream>
#include <fstream>
#include "params.h"
#include "utils.h"
#include "analytical.h"

int main() {
    try {
        // 1. Inicjalizacja siatki przestrzennej - tablica statyczna
        double x[Nx];
        init_space_grid(x);

        // 2. Wybrane momenty czasu do testów (t=0, t=0.1, t=0.5, t=1.0, t=2.0)
        constexpr int Nt_points = 5;
        double time_points[Nt_points] = {0.0, 0.1, 0.5, 1.0, 2.0};

        // 3. Plik wyjściowy dla wyników
        std::ofstream out("analytical_results.csv");
        if (!out) {
            std::cerr << "Blad: Nie mozna otworzyc pliku do zapisu\n";
            return 1;
        }
        out << "x,t,U_analytical\n"; // Nagłówek CSV

        // 4. Obliczenia dla każdego czasu
        double U_analytical[Nx];
        for (int k = 0; k < Nt_points; ++k) {
            double t = time_points[k];
            calculate_analytical_solution(x, Nx, t, U_analytical);

            // Zapis wyników
            for (int i = 0; i < Nx; ++i) {
                out << x[i] << "," << t << "," << U_analytical[i] << "\n";
            }
        }

        std::cout << "Obliczenia analityczne zapisane do analytical_results.csv\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Blad: " << e.what() << "\n";
        return 1;
    }
}
