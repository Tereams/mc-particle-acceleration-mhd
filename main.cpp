#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "Particle.h"
#include "Random.h"
#include "Distributions.h"
#include "UCSModel.h"
#include "FermiModel.h"
#include "HybridModel.h"

// ---- constants (CGS)
constexpr double c  = 2.99792458e10;
constexpr double qe = 4.803204e-10;
constexpr double mp = 1.67262e-24;

constexpr double erg2eV = 6.2415e11;

// ---- transport parameters (match your MATLAB defaults)
constexpr double L  = 1e9;   // [cm]
constexpr double Dr_min = 1e2;
constexpr double Dr_max = L;
constexpr double P_law_index = 1.2;

// sample step length (same transform as MATLAB)
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

    // ---- User-chosen simulation scale
    const int nP = 10000;     // start smaller than 1e6 for C++ bring-up
    const int Nsteps = 5000;  // per particle

    // ---- plasma parameters (match MATLAB defaults)
    double B0 = 100.0; // Gauss
    double n0 = 1e9;   // cm^-3
    double VA = B0 / std::sqrt(4 * M_PI * n0 * mp);

    // particle: ion, Zion=1
    double m = mp;
    double q = qe;
    double Wrest = m * c * c;

    // initial energy: kT0=100 eV -> erg
    double kT0_eV = 100.0;
    double W0 = kT0_eV / erg2eV;
    double v0 = std::sqrt(2 * W0 / m);

    // ---- Choose model (match your intent)
    // imodel = 2 or 12 (2=UCS only, 12=Fermi+UCS)
    int imodel = 2;

    // ---- UCS params (reconnecting + non-reconnecting)
    UCSParams up;
    up.q = q; up.m = m; up.c = c; up.VA = VA;

    // Reconnecting UCS
    up.BLim1 = 1e-5;
    up.BLim2 = 1e2;
    up.BPLa  = 5.0/3.0;

    // lBdepend linear Leff(Eeff)
    double LeffLim1 = 1e2;
    double LeffLim2 = 1e5;
    up.alpha1 = c * (LeffLim2 - LeffLim1) / (VA * (up.BLim2 - up.BLim1));
    up.beta1  = LeffLim1 - up.BLim1 * (LeffLim2 - LeffLim1) / (up.BLim2 - up.BLim1);

    // Non-reconnecting UCS switches
    up.enable_nonreconn = true; // Non_reconn_UCS = 1
    up.P_reconnect = 1.0;       // set <1 to enable non-reconnecting branch

    // Provide ED from MATLAB if you want strict match
    up.ED = 1.0; // TODO: paste MATLAB ED value here for exact consistency

    up.lEplaw_non = true;
    up.Eeff_normal = false;
    up.EPLa_non = 5.0/3.0;
    up.EeffParamLim1_non = 1e2;
    up.EeffParamLim2_non = 1e4;

    up.LeffLim1_non = 1e2;
    up.LeffLim2_non = 1e5;
    up.lLeffrand_non = false; // if true -> varies; else constant
    up.Leff_plaw = false;

    UCSModel ucs(up);

    // ---- Fermi model params (for imodel=12)
    double Vfactor = 1.0;
    double V = Vfactor * VA;
    double GammaV = 1.0 / std::sqrt(1.0 - (V*V)/(c*c));

    FermiParams fp;
    fp.V = V; fp.GammaV = GammaV; fp.c = c; fp.Wrest = Wrest;
    FermiModel fermi(fp);

    // ---- Hybrid model 12 (Fermi + UCS)
    HybridParams hp;
    hp.PdW = 0.095; // MATLAB PdW
    hp.magnetic_field_decay = false;
    hp.PdW_init = 0.095;
    hp.t_decay = 1e-11;

    HybridModel12 hybrid12(&fermi, &ucs, hp);

    // ---- Output CSV
    std::ofstream out("results.csv");
    out << "ip,escaped,t_final,W_final_eV,kkicks,sum_dW_eV\n";

    int escaped_count = 0;

    for (int ip = 0; ip < nP; ++ip) {
        Particle p;
        p.x = 0; p.y = 0; p.z = 0;
        p.v = v0;
        p.W = W0;
        p.t = 0.0;
        p.pdir = {1,0,0};

        for (int step = 0; step < Nsteps; ++step) {
            // random direction in 3D (continuous for transport)
            double costheta = 2.0 * rng.uniform() - 1.0;
            double sintheta = std::sqrt(std::max(0.0, 1 - costheta * costheta));
            double phi = 2 * M_PI * rng.uniform();

            double dx = sintheta * std::cos(phi);
            double dy = sintheta * std::sin(phi);
            double dz = costheta;

            double dr = sample_step(rng);

            // position update
            p.x += dr * dx;
            p.y += dr * dy;
            p.z += dr * dz;

            // time update (event-driven)
            double tau = dr / std::abs(p.v);
            p.t += tau;

            // open boundary check
            if (std::abs(p.x) > L/2 || std::abs(p.y) > L/2 || std::abs(p.z) > L/2) {
                p.escaped = true;
                break;
            }

            // scattering event
            ScatterResult sr;
            if (imodel == 2) {
                sr = ucs.scatter(p, rng);
            } else if (imodel == 12) {
                sr = hybrid12.scatter(p, rng);
            } else {
                // keep minimal; you can add imodel==1,3,13 later
                sr = ucs.scatter(p, rng);
            }

            p.kkicks++;
            p.sum_dW += sr.dW;

            // update energy and time
            p.W += sr.dW;
            p.t += sr.dtScat;

            // velocity update (relativistic switch)
            double gamma = 1.0 + p.W / Wrest;
            if (gamma > 1.00001) p.v = c * std::sqrt(1 - 1 / (gamma * gamma));
            else p.v = std::sqrt(2 * p.W / m);

            // MATLAB: after scattering, pdir is re-randomized on lattice
            p.pdir = random_lattice_dir(rng);
        }

        if (p.escaped) escaped_count++;

        out << ip << ","
            << (p.escaped ? 1 : 0) << ","
            << p.t << ","
            << (p.W * erg2eV) << ","
            << p.kkicks << ","
            << (p.sum_dW * erg2eV)
            << "\n";
    }

    out.close();

    std::cout << "Done. Escaped: " << escaped_count << " / " << nP << "\n";
    std::cout << "Wrote results.csv\n";
    return 0;
}
