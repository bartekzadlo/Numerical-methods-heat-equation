#pragma once

#include <vector>
#include <string>

// Zapisuje profil rozwiązania wraz z błędem
void save_profile(const std::string& filename,
                 const std::vector<double>& u,
                 double dx,
                 double time);

// Zapisuje maksymalny błąd w czasie
void save_max_error(const std::string& filename,
                   const std::vector<double>& max_errors,
                   double dt);

// Oblicza maksymalny błąd bezwzględny
double calculate_max_error(const std::vector<double>& u,
                         double dx,
                         double time);

// Funkcja zapisuje do pliku tekstowego wartości logarytmów kroku przestrzennego (log_dx)
// oraz odpowiadających im logarytmów maksymalnych błędów (log_err).
void save_log_error_plot(const std::string& filename,
                         const std::vector<double>& log_dx,
                         const std::vector<double>& log_err)