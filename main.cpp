#include <iostream>
#include <cmath>
#include "Particle.h"
#include "Random.h"
#include "UCSModel.h"

// ---- constants (CGS)
constexpr double c  = 2.99792458e10;
constexpr double qe = 4.803204e-10;
constexpr double mp = 1.67262e-24;

// ---- simulation parameters
constexpr double B0 = 100.0;     // Gauss
constexpr double n0 = 1e9;       // cm^-3
constexpr double L  = 1e9;       // box size [cm]

constexpr double kT0_eV = 100.0;
constexpr double erg2eV = 6.2415e11;

constexpr double Dr_min = 1e2;
constexpr double Dr_max = L;
constexpr double P_law_index = 1.2;

double sample_step(Random& rng) {
    double u = rng.uniform();
    return std::pow(
        std::pow(Dr_min, 1 - P_law_index) +
        (std::pow(Dr_max, 1 - P_law_index) - std::pow(Dr_min, 1 - P_law_index)) * u,
        1.0 / (1.0 - P_law_index)
    );
}

int main() {
    Random rng(123);

    // --- derived quantities
    double VA = B0 / std::sqrt(4 * M_PI * n0 * mp);

    double m = mp;
    double q = qe;

    double Wrest = m * c * c;

    // initial thermal velocity
    double Wth = kT0_eV / erg2eV;
    double v0 = std::sqrt(2 * Wth / m);

    // --- UCS model params (same as MATLAB)
    double BLim1 = 1e-5;
    double BLim2 = 1e2;
    double BPLa  = 5.0 / 3.0;

    double LeffLim1 = 1e2;
    double LeffLim2 = 1e5;

    double alpha1 = c * (LeffLim2 - LeffLim1) / (VA * (BLim2 - BLim1));
    double beta1  = LeffLim1 - BLim1 * (LeffLim2 - LeffLim1) / (BLim2 - BLim1);

    UCSModel model(q, VA, c, BLim1, BLim2, BPLa, alpha1, beta1);

    // ---- particle
    Particle p;
    p.x = 0.0;
    p.y = 0.0;
    p.z = 0.0;
    p.v = v0;
    p.W = Wth;
    p.t = 0.0;

    // ---- Monte Carlo loop
    const int Nsteps = 10000;

    for (int step = 0; step < Nsteps; ++step) {
        // random direction
        double costheta = 2.0 * rng.uniform() - 1.0;
        double sintheta = std::sqrt(1 - costheta * costheta);
        double phi = 2 * M_PI * rng.uniform();

        double dx = sintheta * std::cos(phi);
        double dy = sintheta * std::sin(phi);
        double dz = costheta;

        double dr = sample_step(rng);

        // position update
        p.x += dr * dx;
        p.y += dr * dy;
        p.z += dr * dz;

        // time update
        double tau = dr / std::abs(p.v);
        p.t += tau;

        // boundary check (open)
        if (std::abs(p.x) > L / 2 ||
            std::abs(p.y) > L / 2 ||
            std::abs(p.z) > L / 2) {
            p.escaped = true;
            break;
        }

        // scattering (UCS)
        double dW = model.compute_dW(p, rng);
        p.W += dW;
        p.sum_dW += dW;
        p.kkicks++;

        // update velocity
        double gamma = 1.0 + p.W / Wrest;
        if (gamma > 1.00001)
            p.v = c * std::sqrt(1 - 1 / (gamma * gamma));
        else
            p.v = std::sqrt(2 * p.W / m);
    }

    std::cout << "Escaped: " << p.escaped << "\n";
    std::cout << "Final time [s]: " << p.t << "\n";
    std::cout << "Final energy [eV]: " << p.W * erg2eV << "\n";
    std::cout << "Number of kicks: " << p.kkicks << "\n";

    return 0;
}
