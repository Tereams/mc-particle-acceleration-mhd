#pragma once
#include <array>

struct Particle {
    // position [cm]
    double x = 0, y = 0, z = 0;

    // kinematics
    double v = 0;     // speed [cm/s]
    double W = 0;     // kinetic energy [erg]
    double t = 0;     // time [s]

    // direction on lattice (MATLAB: pdir has one ±1 entry, others 0)
    std::array<int,3> pdir{1,0,0};

    // statistics
    int kkicks = 0;
    double sum_dW = 0.0;
    bool escaped = false;
};
