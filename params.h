#pragma once
#include <cmath>
namespace params {
    constexpr double D = 1.0;
    constexpr double t_max = 2.0;

    inline double get_a() { return 6.0 * std::sqrt(D * t_max); }
}
