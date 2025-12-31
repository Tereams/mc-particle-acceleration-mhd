#pragma once
#include <cmath>
#include "ScatteringModel.h"

struct FermiParams {
    double V = 0.0;      // cloud speed [cm/s]
    double GammaV = 1.0; // Lorentz factor for cloud
    double c = 0.0;
    double Wrest = 0.0;  // rest energy [erg]
};

inline std::array<int,3> random_lattice_dir(Random& rng) {
    // one axis chosen uniformly, sign ±1
    int axis = static_cast<int>(rng.uniform(0.0, 3.0)); // 0,1,2
    int sgn = (rng.uniform() > 0.5) ? 1 : -1;
    std::array<int,3> d{0,0,0};
    d[axis] = sgn;
    return d;
}

class FermiModel : public ScatteringModel {
public:
    explicit FermiModel(FermiParams p) : P(p) {}

    ScatterResult scatter(const Particle& part, Random& rng) override {
        // MATLAB:
        // sdir random lattice direction
        // isgn = dot(pdir,sdir)
        // dW = 2*GammaV*WiTot*(V^2/c^2 - isgn*(V*v)/c^2)

        auto sdir = random_lattice_dir(rng);
        int isgn = part.pdir[0]*sdir[0] + part.pdir[1]*sdir[1] + part.pdir[2]*sdir[2];

        double Wi = part.W;
        double WiTot = Wi + P.Wrest;
        double c2 = P.c * P.c;

        ScatterResult out;
        out.dW = 2.0 * P.GammaV * WiTot * ( (P.V*P.V)/c2 - isgn * (P.V * part.v)/c2 );
        out.dtScat = 0.0;
        return out;
    }

private:
    FermiParams P;
};
