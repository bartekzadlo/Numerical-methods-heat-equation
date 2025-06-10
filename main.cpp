#include "analytical.h"
#include "params.h"
#include "explicit.h"
#include "implicit.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <functional>
#include "output.h"



void save_time_error_data(const std::string& filename,
                         const std::vector<double>& time_points,
                         const std::vector<double>& errors) {
    std::ofstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << "# t max_error\n";
    for (size_t i = 0; i < time_points.size(); ++i) {
        file << time_points[i] << " " << errors[i] << "\n";
    }
    file.close();
}

std::string generate_time_error_data(const std::string& methodName,
                                    std::vector<double> u_initial,
                                    double dx,
                                    double dt,
                                    const std::vector<double>& time_points,
                                    double D,
                                    std::function<void(std::vector<double>&, double, double, int, double)> solve_method) {

    if (u_initial.empty()) {
        return "";
    }

    std::vector<double> errors_over_time(time_points.size());
    std::vector<double> current_u_state;

    std::cout << "\n[" << methodName[0] << "] " << methodName << "...\n";

    for (size_t i = 0; i < time_points.size(); ++i) {
        double t_target = time_points[i];

        current_u_state = u_initial;

        int steps_to_reach_target = static_cast<int>(t_target / dt + 1e-9);
        double actual_t_reached = steps_to_reach_target * dt;

        solve_method(current_u_state, dx, dt, steps_to_reach_target, D);

        errors_over_time[i] = calculate_max_error(current_u_state, dx, actual_t_reached);

        std::cout << "t = " << actual_t_reached
                  << ", max error = " << errors_over_time[i] << "\n";
    }

    std::string plik_danych = "time_error_" + methodName + ".dat";
    save_time_error_data(plik_danych, time_points, errors_over_time);
    return plik_danych;
}

void run_time_error_analysis() {
    const double D = params::D;
    const double t_max = params::t_max;
    const double dx = 0.0249567;
    const int Nx = 681;
    const double dt_explicit = 0.000249135;
    const double dt_implicit = 0.000622837;
    const double a = params::calculate_a();

    const int num_time_points = static_cast<int>(t_max / 0.1) + 1;
    std::vector<double> time_points(num_time_points);
    for (int i = 0; i < num_time_points; ++i) {
        time_points[i] = i * 0.1;
    }

    std::vector<double> u_initial(Nx);
    double x_start = -a;
    for (int j = 0; j < Nx; ++j) {
        double x = x_start + j * dx;
        double u0 = (x < 0.0) ? 1.0 : (x > 0.0 ? 0.0 : 0.5);
        u_initial[j] = u0;
    }

    std::cout << "Rozpoczynam analizę błędów w czasie...\n";

    std::string explicit_file = generate_time_error_data(
        "explicit", u_initial, dx, dt_explicit, time_points, D, solve_explicit);

    std::string thomas_file = generate_time_error_data(
        "thomas", u_initial, dx, dt_implicit, time_points, D, solve_implicit_thomas);

    std::string lu_file = generate_time_error_data(
        "lu", u_initial, dx, dt_implicit, time_points, D, solve_implicit_lu);

    std::cout << "\nAnaliza zakończona. Dane zapisane do plików:\n";
    std::cout << " - " << explicit_file << "\n";
    std::cout << " - " << thomas_file << "\n";
    std::cout << " - " << lu_file << "\n";
}

int main() {
    try {
        auto total_start = std::chrono::high_resolution_clock::now();
        run_time_error_analysis();
        auto total_end = std::chrono::high_resolution_clock::now();

        std::cout << "\nTotal execution time: "
                 << std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start).count()
                 << " seconds\n";
    } catch (const std::exception& e) {
        std::cerr << "\nERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}