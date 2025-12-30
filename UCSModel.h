#pragma once
#include "ScatteringModel.h"
#include <cmath>

class UCSModel : public ScatteringModel {
public:
    UCSModel(
        double q_,
        double VA_,
        double c_,
        double BLim1_,
        double BLim2_,
        double BPLa_,
        double alpha1_,
        double beta1_
    )
        : q(q_), VA(VA_), c(c_),
          BLim1(BLim1_), BLim2(BLim2_),
          BPLa(BPLa_), alpha1(alpha1_), beta1(beta1_) {}

    double compute_dW(const Particle& p, Random& rng) override {
        // --- sample Beff from power law
        double u = rng.uniform();
        double Beff = std::pow(
            std::pow(BLim1, 1 - BPLa) +
            (std::pow(BLim2, 1 - BPLa) - std::pow(BLim1, 1 - BPLa)) * u,
            1.0 / (1.0 - BPLa)
        );

        // Eeff = Beff * VA / c
        double Eeff = Beff * VA / c;

        // Leff = alpha * Eeff + beta
        double Leff = alpha1 * Eeff + beta1;

        // sign = +1 (reconnecting UCS)
        double dW = std::abs(q) * Eeff * Leff;
        return dW;
    }

private:
    double q, VA, c;
    double BLim1, BLim2, BPLa;
    double alpha1, beta1;
};
