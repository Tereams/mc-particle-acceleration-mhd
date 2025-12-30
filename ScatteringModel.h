#pragma once

struct Particle {
    // position
    double x, y, z;

    // kinematics
    double v;     // speed [cm/s]
    double W;     // kinetic energy [erg]
    double t;     // time [s]

    // statistics
    int kkicks = 0;
    double sum_dW = 0.0;

    bool escaped = false;
};
