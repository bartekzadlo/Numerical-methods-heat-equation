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

#endif