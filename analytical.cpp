#include "analytical.h"
#include "params.h"

double analytical_solution(double x, double t) {
    if (t <= 0.0) {
        return (x < 0.0) ? 1.0 : 0.0;
    }
    long double denom = 2.0L * std::sqrt(static_cast<long double>(params::D) * static_cast<long double>(t));
    long double z = static_cast<long double>(x) / denom;

    long double val = calerfpack::erfc_LD(z);
    return 0.5 * static_cast<double>(val);
}