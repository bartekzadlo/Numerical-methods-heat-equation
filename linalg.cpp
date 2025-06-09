#include "linalg.h"

void thomas_solver(const std::vector<double>& a,
                   const std::vector<double>& b,
                   const std::vector<double>& c,
                   std::vector<double>& d) {
    int n = static_cast<int>(b.size());
    std::vector<double> c_prime(n, 0.0);

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double m = b[i] - a[i - 1] * c_prime[i - 1];
        if (i < n - 1)
            c_prime[i] = c[i] / m;
        d[i] = (d[i] - a[i - 1] * d[i - 1]) / m;
    }

    // Back substitution
    for (int i = n - 2; i >= 0; --i) {
        d[i] = d[i] - c_prime[i] * d[i + 1];
    }
}
