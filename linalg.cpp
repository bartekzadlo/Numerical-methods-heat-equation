#include "linalg.h"
#include <stdexcept>

void thomas_solver(const std::vector<double>& a, const std::vector<double>& b,
                   const std::vector<double>& c, std::vector<double>& d) {
    const int n = b.size();
    std::vector<double> c_prime(n);
    std::vector<double> d_prime(n);

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double m = b[i] - a[i-1] * c_prime[i-1];
        c_prime[i] = c[i] / m;
        d_prime[i] = (d[i] - a[i-1] * d_prime[i-1]) / m;
    }

    d[n-1] = d_prime[n-1];
    for (int i = n-2; i >= 0; --i) {
        d[i] = d_prime[i] - c_prime[i] * d[i+1];
    }
}

void lu_decompose(const std::vector<std::vector<double>>& A,
                  std::vector<std::vector<double>>& L,
                  std::vector<std::vector<double>>& U) {
    const int n = A.size();
    L.assign(n, std::vector<double>(n, 0.0));
    U.assign(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
                sum += L[i][j] * U[j][k];
            U[i][k] = A[i][k] - sum;
        }

        L[i][i] = 1.0;
        for (int k = i+1; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
                sum += L[k][j] * U[j][i];
            L[k][i] = (A[k][i] - sum) / U[i][i];
        }
    }
}

std::vector<double> lu_solve(const std::vector<std::vector<double>>& L,
                             const std::vector<std::vector<double>>& U,
                             const std::vector<double>& b) {
    const int n = b.size();
    std::vector<double> y(n), x(n);

    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j)
            sum += L[i][j] * y[j];
        y[i] = b[i] - sum;
    }

    for (int i = n-1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i+1; j < n; ++j)
            sum += U[i][j] * x[j];
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}