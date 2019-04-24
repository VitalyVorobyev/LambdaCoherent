/** Reduced matrix element squared for J/psi -> [Lambda -> p pi-] [anti-Lambda -> pbar pi+]
 *
 *  The matrix elemnt takes into account electron polarization.
 *  The expressions are obtained by A. Grabovsky and A. Reznichenko (BINP, Novosibirk)
 * 
 *  Coded by V. Vorobyev (BINP)
 *
 *  Date: 23 April 2019
 *
 **/

#ifndef COHMSQLL_H__
#define COHMSQLL_H__

// libLinal
#include "lvect.h"

#include <array>

/** Fourth rank Levi-Civita tensor **/
class LeviCivita {
    std::array<int, 256> m_lc;
    static int levi();
    static int idx;

 public:
    LeviCivita();

    int operator()(int i, int j, int k, int l) const;
};

class CohMsqLL {
    /** J/psi mass **/
    static constexpr double m_mJpsi   = 3.09692;
    static constexpr double m_mLambda = 1.11568;
    static constexpr double m_mproton = 0.93827;
    static constexpr double m_mpion   = 0.13957;
    static constexpr double m_mLambdaSq = m_mLambda * m_mLambda;
    static constexpr double m_lp = 0.5 * (m_mLambdaSq + m_mproton*m_mproton - m_mpion*m_mpion);
    static constexpr double m_mJpsiSq = m_mJpsi * m_mJpsi;

    LeviCivita m_lc;
    static constexpr double m_eps = 1.e-6;
    double Vol(const linal::LVect& v1, const linal::LVect& v2,
               const linal::LVect& v3, const linal::LVect& v4) const;

    /** electrons polarization level **/
    const double m_xi;
    /** ratio of absolute values of the Lambda formfactors **/
    const double m_alpha;
    /** cos relative phase of the Lambda formfactors **/
    const double m_cosphi;
    /** sin relative phase of the Lambda formfactors **/
    const double m_sinphi;
    /** Lambda0 decay parameter **/
    const double m_alpha1;
    /** anti-Lambda0 decay parameter **/
    const double m_alpha2;

    /** Total (J/psi) momentum **/
    const linal::LVect m_P;
    /** Lambda and anti-Lambda momenta difference **/
    mutable linal::LVect m_Q;
    /** Positron momentum **/
    const linal::LVect m_kp;
    /** Electron momentum **/
    const linal::LVect m_kn;
    /** Proton momentum **/
    mutable linal::LVect m_l1;
    /** anti-Proton momentum **/
    mutable linal::LVect m_l2;
    /** Lambda momentum **/
    mutable linal::LVect m_p1;
    /** anti-Lambda momentum **/
    mutable linal::LVect m_p2;

    // Precomputed values
    /** 5*M**5 **/
    const double m_5mLam5;
    /** 4 * alpha * cos(phi) */
    const double m_4alpCosPhi;
    /** 4 * alpha * sin(phi) */
    const double m_4alpSinPhi;
    /** alpha**2 **/
    const double m_alphaSq;
    /** Q**2 **/
    mutable double m_Qsq;
    /** dot(kp, Q) **/
    mutable double m_kpQ;
    /** dot(kp, l1) **/
    mutable double m_kpl1;
    /** dot(kn, l1) **/
    mutable double m_knl1;
    /** dot(kp, l2) **/
    mutable double m_kpl2;
    /** dot(kn, l2) **/
    mutable double m_knl2;
    /** dot(l1, P) **/
    mutable double m_l1P;
    /** dot(l2, P) **/
    mutable double m_l2P;

    double a() const;
    double b() const;

    double aRR() const;
    double aRSSR() const;
    double aSS() const;

    double bRS() const;
    double bSR() const;
    double bSS() const;

 public:
    /** Constructor
        @args:
            xi - electrons polarization level
            alpha - ratio of absolute values of the Lambda formfactors
            phi - relative phase of the Lambda formfactors
            alpha1 - Lambda0 decay parameter
            alpha2 - anti-Lambda0 decay parameter
     **/
    CohMsqLL(double xi, double alpha, double phi, double alpha1, double alpha2);
    
    /** Calcualtes matrix element **/
    double operator()(const linal::LVect& p, const linal::LVect& pbar,
               const linal::LVect& pip, const linal::LVect& pin) const;
};

 #endif  // COHMSQLL_H__
