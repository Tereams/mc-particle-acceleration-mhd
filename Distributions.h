#pragma once
#include <cmath>
#include "Random.h"

// Sample x from power-law between [xmin, xmax] with exponent a:
// pdf ~ x^{-a}  (MATLAB code uses: x = (xmin^(1-a) + (xmax^(1-a)-xmin^(1-a))*u)^(1/(1-a)))
inline double sample_powerlaw(Random& rng, double xmin, double xmax, double a) {
    double u = rng.uniform();
    double p1 = std::pow(xmin, 1.0 - a);
    double p2 = std::pow(xmax, 1.0 - a);
    return std::pow(p1 + (p2 - p1) * u, 1.0 / (1.0 - a));
}

// Truncated normal with hard cutoff [lo, hi]
inline double sample_trunc_normal(Random& rng, double mean, double sigma, double lo, double hi) {
    // Simple rejection sampling; matches MATLAB while-loop logic
    double x = mean;
    do {
        x = mean + sigma * rng.normal01();
    } while (x < lo || x > hi);
    return x;
}

// Random sign ±1
inline int sample_sign(Random& rng) {
    return (rng.uniform() > 0.5) ? 1 : -1;
}
