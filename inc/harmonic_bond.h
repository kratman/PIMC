
#pragma once

#include "physical_constants.h"
#include "particle.h"
#include "miscellaneous_constants.h"

struct BondParameters
{
    int at2; //atom ID 2
    double K; //Force constant (eV/Ang^2)
    double R; //Equil. distance (Ang)
};

inline double HarmEnergy(BondParameters& bond, Qpart& at1, Qpart& at2, int p)
{
    double r = sqrt(Dist2(at1.P[p],at2.P[p]));
    if ((r == 0.0) and (CoulSub == 1))
    {
        return HugeNum; //Avoid dividing by zero
    }
    double E = 0.5*bond.K*(r-bond.R)*(r-bond.R);
    //Coulomb subtract, if applicable
    if (CoulSub == 1)
    {
        E -= (C2eV*at1.q*at2.q/r);
    }
    return E;
};
