#include "LLKine.h"

#include <cmath>

using std::array;

using linal::Vect;
using linal::LVect;

const Vect LLKine::c_ebeam(0., 0., 1.);

// Cosine of polar angle of Lambda in the global frame
double LLKine::cosLambda(const LVect& pr, const LVect& pin) {
    return dot(c_ebeam, (pr+pin).Vec().unit());
}

// Angles in the Lambda frame
array<double, 2> LLKine::omegaLambda(LVect pr, const LVect& pin) {
    LVect lamP4(pr + pin);
    Vect ez = lamP4.Vec().unit();
    Vect ey = cross(ez, c_ebeam).unit();
    Vect boost = -lamP4.BoostVec();
    Vect pdir = pr.Boost(boost).Vec().unit();

    double costh = dot(pdir, ez);
    double sinphi = dot(pdir, ey) / sqrt(1. - costh*costh);

    return {costh, sinphi};
}

// Complete five phase space parameters
array<double, 5> LLKine::xi(const LVect& pr, const LVect& pin, const LVect& pbar, const LVect& pip) {
    auto om1 = omegaLambda(pr, pin);
    auto om2 = omegaLambda(pbar, pip);
    return {cosLambda(pr, pin), om1[0], om1[1], om2[0], om2[1]};
}
