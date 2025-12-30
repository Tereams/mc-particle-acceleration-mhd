#pragma once
#include <random>

class Random {
public:
    Random(unsigned seed = 42) : gen(seed), uni(0.0, 1.0) {}

    double uniform() {
        return uni(gen);
    }

    double normal() {
        return norm(gen);
    }

private:
    std::mt19937 gen;
    std::uniform_real_distribution<double> uni;
    std::normal_distribution<double> norm;
};
