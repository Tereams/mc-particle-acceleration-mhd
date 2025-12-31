#pragma once
#include "ScatteringModel.h"

struct HybridParams {
    double PdW = 0.0;                 // probability to choose model 1 vs model 2 (see below)
    bool magnetic_field_decay = false;
    double PdW_init = 0.095;
    double t_decay = 1e-11;
};

class HybridModel12 : public ScatteringModel {
public:
    HybridModel12(ScatteringModel* model1, ScatteringModel* model2, HybridParams hp)
        : m1(model1), m2(model2), H(hp) {}

    ScatterResult scatter(const Particle& p, Random& rng) override {
        double PdW_eff = H.PdW;

        // MATLAB:
        // if magnetic_field_decay: PdW = PdW_init * log10(t/t_decay)
        if (H.magnetic_field_decay) {
            if (p.t > 0 && H.t_decay > 0) {
                PdW_eff = H.PdW_init * std::log10(p.t / H.t_decay);
                // keep PdW in [0,1]
                if (PdW_eff < 0) PdW_eff = 0;
                if (PdW_eff > 1) PdW_eff = 1;
            }
        }

        // MATLAB imodel==12:
        // idWChoice = (rand > PdW) + 1  -> PdW chooses model1 with prob PdW? (depends on convention)
        // Their code: (rand > PdW)+1 => if rand > PdW => 2 else 1
        // So model1 chosen with probability PdW, model2 with probability (1-PdW).
        if (rng.uniform() <= PdW_eff) return m1->scatter(p, rng);
        return m2->scatter(p, rng);
    }

private:
    ScatteringModel* m1;
    ScatteringModel* m2;
    HybridParams H;
};
