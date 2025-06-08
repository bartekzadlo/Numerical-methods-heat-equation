// main.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include "params.h"
#include "utils.h"
#include "analytical.h"

int main() {
    try {
        // 1. Inicjalizacja siatki przestrzennej
        std::vector<double> x = init_space_grid();

        // 2. Wybrane momenty czasu do testów (t=0, t=0.1, t=0.5, t=1.0, t=2.0)
        const std::vector<double> time_points = {0.0, 0.1, 0.5, 1.0, 2.0};

        // 3. Plik wyjściowy dla wyników
        std::ofstream out("analytical_results.csv");
        out << "x,t,U_analytical\n"; // Nagłówek CSV

        // 4. Obliczenia dla każdego czasu
        for (double t : time_points) {
            std::vector<double> U_analytical(x.size());
            calculate_analytical_solution(x.data(), x.size(), t, U_analytical.data());

            // Zapis wyników
            for (int i = 0; i < x.size(); ++i) {
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