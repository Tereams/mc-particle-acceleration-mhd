#pragma once
#include "Particle.h"
#include "Random.h"

struct ScatterResult {
    double dW = 0.0;     // energy gain [erg]
    double dtScat = 0.0; // extra time [s] (MATLAB: dtScat)
};

class ScatteringModel {
public:
    virtual ~ScatteringModel() = default;
    virtual ScatterResult scatter(const Particle& p, Random& rng) = 0;
};
