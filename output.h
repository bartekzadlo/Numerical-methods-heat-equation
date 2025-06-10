#ifndef OUTPUT_H
#define OUTPUT_H

#include <vector>
#include <string>

void save_log_error_data(const std::string& filename,
                        const std::vector<double>& h_values,
                        const std::vector<double>& errors);

double calculate_max_error(const std::vector<double>& u,
                         double dx,
                         double t);

void save_solutions_for_time(const std::string& filename,
                             const std::vector<double>& x,
                             const std::vector<double>& u_explicit,
                             const std::vector<double>& u_thomas,
                             const std::vector<double>& u_lu,
                             double t);

#endif
