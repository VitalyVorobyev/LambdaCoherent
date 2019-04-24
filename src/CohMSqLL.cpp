#include "CohMSqLL.h"

#include <iostream>
#include <cmath>
#include <algorithm>  // std::generate

using linal::LVect;
using linal::dot;

using std::cout;
using std::endl;
using std::cerr;

int LeviCivita::idx;

int LeviCivita::levi() {
    int i =  idx &   3;        // 0b00000011
    int j = (idx &  12) >> 2;  // 0b00001100
    int k = (idx &  48) >> 4;  // 0b00110000
    int l = (idx & 192) >> 6;  // 0b11000000
    idx++;

    int result = (i-j) * (i-k) * (i-l) * (j-k) * (j-l) * (k-l) / 12;
    // cout << "levi: " << idx << " -> (" << i << ", " << j << ", " << k << ", " << l << "): " << result << endl;
    return -result;
}

LeviCivita::LeviCivita() {
    idx = 0;
    std::generate(m_lc.begin(), m_lc.end(), LeviCivita::levi);
}

int LeviCivita::operator()(int i, int j, int k, int l) const {
    return m_lc[i + (j << 2) + (k << 4) + (l << 6)];
}

double CohMsqLL::Vol(const LVect& v1, const LVect& v2, const LVect& v3, const LVect& v4) const {
    double res = 0;
    for (int i = 0; i < 4; i++) {
        if (fabs(v1.at(i)) < m_eps)
            continue;
        for (int j = 0; j < 4; j++) {
            if ((i == j) || (fabs(v2.at(j)) < m_eps))
                continue;
            for (int k = 0; k < 4; k++) {
                if ((k == i) || (k == j) || (fabs(v3.at(k)) < m_eps))
                    continue;
                for (int l = 0; l < 4; l++) {
                    if ((l == i) || (l == j) || (l == k) || (fabs(v4.at(l)) < m_eps))
                        continue;
                    res += m_lc(i, j, k, l) * v1.at(i) * v2.at(j) * v3.at(k) * v4.at(l);
                }
            }
        }
    }
    return res;
}

CohMsqLL::CohMsqLL(double xi, double alpha, double phi, double alpha1, double alpha2) :
    m_xi(xi), m_alpha(alpha), m_cosphi(cos(phi)), m_sinphi(sin(phi)),
    m_alpha1(alpha1), m_alpha2(alpha2),
    m_P(m_mJpsi, 0., 0., 0.),
    m_kp(0.5 * m_mJpsi, 0., 0., -0.5 * m_mJpsi),
    m_kn(0.5 * m_mJpsi, 0., 0.,  0.5 * m_mJpsi),
    m_5mLam5(5. * pow(m_mLambda, 5)),
    m_4alpCosPhi(4. * m_alpha * m_cosphi),
    m_4alpSinPhi(4. * m_alpha * m_sinphi),
    m_alphaSq(m_alpha*m_alpha)
    {}

double CohMsqLL::operator()(const LVect& p, const LVect& pbar, const LVect& pip, const LVect& pin) const {
    m_l1 = p;
    m_l2 = pbar;
    m_p1 = p + pin;
    m_p2 = pbar + pin;
    m_Q = m_p1 - m_p2;
    m_Qsq = m_Q.m2();
    m_kpQ = dot(m_kp, m_Q);
    m_knl1 = dot(m_kn, m_l1);
    m_knl2 = dot(m_kn, m_l2);
    m_kpl1 = dot(m_kp, m_l1);
    m_kpl2 = dot(m_kp, m_l2);
    m_l1P = dot(m_l1, m_P);
    m_l2P = dot(m_l2, m_P);

    double aa = a();
    double bb = b();
    cout << "a: " << aa << ", b: " << bb << endl;
    return aa + m_xi * bb;
}

double CohMsqLL::a() const {
    return aRR() + aRSSR() + m_alpha1 * m_alpha2 * aSS();
}

double CohMsqLL::b() const {
    return m_alpha1 * bRS() + m_alpha2 * bSR() + m_alpha1 * m_alpha2 * bSS();
}

double CohMsqLL::aRR() const {
    double dblkpQsq = 2. * m_kpQ * m_kpQ;
    return 2. * dblkpQsq + m_mJpsiSq * (4. * m_mLambdaSq + m_mJpsiSq) - 
        (dblkpQsq + 0.5 * m_Qsq * m_mJpsiSq) * (m_4alpCosPhi - 0.5 * m_alphaSq * m_Qsq / m_mLambdaSq);
}

double CohMsqLL::aRSSR() const {
    return m_4alpSinPhi * m_kpQ * Vol(m_kn - m_kp, m_alpha1*m_l1 + m_alpha2*m_l2, m_p1, m_p2);
}

double CohMsqLL::aSS() const {
    double l1p2 = dot(m_l1, m_p2);
    double l2p1 = dot(m_l2, m_p1);
    double l1l2 = dot(m_l1, m_l2);
    double l1Q  = dot(m_l1, m_Q);
    double l2Q  = dot(m_l2, m_Q);
    double lpsq = m_lp * m_lp;
    double dblkpQsq = 2. * m_kpQ * m_kpQ;
    double kpnl1 = m_kpl1 - m_knl1;
    double kpnl2 = m_kpl2 - m_knl2;
    double lpnl21 = kpnl2 - kpnl1;

    double aSS1 = (dblkpQsq + 0.5 * m_Qsq * m_mJpsiSq) * (
        0.5 * m_Qsq / m_mLambdaSq * lpsq + l1Q * l2Q - 0.5 * m_Qsq * l1l2
    );

    double aSS2 = m_mLambdaSq * m_mJpsiSq * (m_knl1*m_knl2 + m_kpl1*m_kpl2 - l1p2*l2p1 - 0.25*m_Qsq*l1l2) -
                  m_mLambdaSq * m_kpQ * (l1l2 * m_kpQ - 2.*m_knl2*m_kpl1 + 2.*m_knl1*m_kpl2) -
                  0.25 * lpsq * (2. * dblkpQsq + m_mJpsiSq);

    double lp2l1p2l2p1 = m_lp * (l1p2 + l2p1);
    double msqdiff = m_mJpsiSq - 2.*m_mLambdaSq;
    double aSS3 = dblkpQsq * (2.*m_mLambdaSq*l1l2 - lp2l1p2l2p1) + m_kpQ * (
                      -msqdiff * m_lp * lpnl21 + 2.*m_mLambdaSq * (l1p2*kpnl2 - l2p1*kpnl1)
                    ) + m_mJpsiSq * (msqdiff * lpsq +
                       m_mLambdaSq * (2.*l1p2*l2p1 + m_Qsq * l1l2 - 2. * lp2l1p2l2p1)
                  );

    return m_alphaSq * aSS1 + 4.*aSS2 + 0.5*m_4alpCosPhi*aSS3;
}

double CohMsqLL::bRS() const {
    double kpnl1 = m_kpl1 - m_knl1;
    return m_mJpsiSq * (
        4. * (m_lp * m_kpQ - m_mLambdaSq * kpnl1) +
        0.25 * m_4alpCosPhi * (2.*m_kpQ*(m_l1P - 2.*m_lp + m_Qsq*kpnl1))
    );
}

double CohMsqLL::bSR() const {
    double kpnl2 = m_kpl2 - m_knl2;
    return m_mJpsiSq * (
        -4. * (m_lp * m_kpQ + m_mLambdaSq * kpnl2) +
        0.25 * m_4alpCosPhi * (2.*m_kpQ*(2.*m_lp - m_l2P + m_Qsq*kpnl2))
    );
}

double CohMsqLL::bSS() const {
    auto p4comb = (2. * m_mLambdaSq * m_l2P - m_mJpsiSq * m_lp) * m_l1 -
                  (2. * m_mLambdaSq * m_l1P - m_mJpsiSq * m_lp) * m_l2;
    return 0.5 * m_4alpSinPhi * Vol(m_kn, m_kp, p4comb, m_Q);
}
