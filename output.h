#pragma once

#include <vector>
#include <string>

void save_profile(const std::string& filename, const std::vector<double>& u,
                  double dx, double time);

void save_error(const std::string& filename,
                const std::vector<double>& errors,
                double dt);
