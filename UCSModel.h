#pragma once
#include <cmath>
#include <algorithm>
#include "ScatteringModel.h"
#include "Distributions.h"

struct UCSParams {
    // Common
    double q = 0.0;     // statC
    double m = 0.0;     // g
    double c = 0.0;
    double VA = 0.0;

    // Switch
    bool enable_nonreconn = false; // Non_reconn_UCS
    double P_reconnect = 1.0;      // Probability to encounter reconnecting UCS

    // Reconnecting UCS: Beff power-law
    double BLim1 = 1e-5;
    double BLim2 = 1e2;
    double BPLa  = 5.0/3.0; // exponent in MATLAB sampling

    // Reconnecting: Leff depends on Eeff (linear), from lBdepend branch
    double alpha1 = 0.0;
    double beta1  = 0.0;

    // Non-reconnecting: Eeff distribution relative to Dreicer field ED
    double ED = 1.0;                 // Provide from MATLAB to match exactly
    bool lEplaw_non = true;
    double EPLa_non = 5.0/3.0;
    double EeffParamLim1_non = 1e2;  // times ED
    double EeffParamLim2_non = 1e4;  // times ED

    bool Eeff_normal = false;
    double Eeff_mean_coeff = 1e4;         // * ED
    double Eeef_sigma_coeff = (1.0/800.0) * 1e4; // * ED

    // Non-reconnecting: Leff
    double LeffLim1_non = 1e2;
    double LeffLim2_non = 1e5;
    bool Leff_plaw = false;          // if true, Leff ~ PL
    bool lLeffrand_non = false;      // if false, constant Leff=LeffLim1_non

    // If Leff is linear function of Eeff (as in MATLAB alpha_non/beta_non)
    // We'll compute (alpha_non, beta_non) from limits when needed.
};

class UCSModel : public ScatteringModel {
public:
    explicit UCSModel(UCSParams p) : P(p) {}

    ScatterResult scatter(const Particle& part, Random& rng) override {
        (void)part;

        double Eeff = 0.0;
        double Leff = 0.0;
        int sign1 = 1; // reconnecting fixed +1; non-reconnecting random ±1

        const bool do_nonreconn = P.enable_nonreconn && (rng.uniform() > P.P_reconnect);

        if (!do_nonreconn) {
            // ---- Reconnecting UCS (MATLAB: Beff -> Eeff=Beff*VA/c, Leff=alpha1*Eeff+beta1, sign1=+1)
            double Beff = sample_powerlaw(rng, P.BLim1, P.BLim2, P.BPLa);
            Eeff = Beff * P.VA / P.c;
            Leff = P.alpha1 * Eeff + P.beta1;
            sign1 = 1;
        } else {
            // ---- Non-reconnecting UCS (MATLAB: Eeff from PL/normal/uniform; Leff from linear/PL/const; sign1=±1)
            sign1 = sample_sign(rng);

            // Eeff limits
            const double EeffLim1 = P.EeffParamLim1_non * P.ED;
            const double EeffLim2 = P.EeffParamLim2_non * P.ED;

            if (P.lEplaw_non && !P.Eeff_normal) {
                Eeff = sample_powerlaw(rng, EeffLim1, EeffLim2, P.EPLa_non);
            } else if (P.Eeff_normal) {
                double mean  = P.Eeff_mean_coeff * P.ED;
                double sigma = P.Eeef_sigma_coeff * P.ED;
                Eeff = sample_trunc_normal(rng, mean, sigma, mean - 3*sigma, mean + 3*sigma);
            } else {
                // Uniform/constant fallback: set Eeff at mid or ED*coef; keep simplest
                Eeff = EeffLim2; // you can change to a dedicated Efactor_non*ED if you want 1:1 with MATLAB uniform branch
            }

            // Leff
            if (!P.lLeffrand_non) {
                Leff = P.LeffLim1_non; // constant
            } else if (P.Leff_plaw) {
                Leff = sample_powerlaw(rng, P.LeffLim1_non, P.LeffLim2_non, P.EPLa_non);
            } else {
                // Linear function of Eeff using endpoints
                // alpha_non = (LeffLim2-LeffLim1)/(EeffLim2-EeffLim1)
                // beta_non  = LeffLim1 - EeffLim1*alpha_non
                double alpha_non = (P.LeffLim2_non - P.LeffLim1_non) / (EeffLim2 - EeffLim1);
                double beta_non  = P.LeffLim1_non - EeffLim1 * alpha_non;
                Leff = alpha_non * Eeff + beta_non;
            }
        }

        ScatterResult out;
        out.dW = sign1 * std::abs(P.q) * Eeff * Leff;
        out.dtScat = 0.0; // MATLAB has idtCS*Leff/sqrt(2Wi/m). Keep 0 for now (idtCS=0).
        return out;
    }

private:
    UCSParams P;
};
