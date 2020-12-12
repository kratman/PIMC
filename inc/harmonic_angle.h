
#pragma once

#include "particle.h"
#include "physics_settings.h"

struct AngleParameters
{
    int at2; //atom ID 2
    int at3; //atom ID 3
    double K; //Force constant (eV/rad^2)
    double ang; //Equilibrium angle (Rad)
};

inline double AngEnergy(AngleParameters& ang, Qpart& at1, Qpart& at2, Qpart& at3, int p)
{
    double r1 = sqrt(Dist2(at1.P[p],at2.P[p]));
    double r2 = sqrt(Dist2(at1.P[p],at3.P[p]));
    double dx1 = at2.P[p].x-at1.P[p].x;
    double dx2 = at3.P[p].x-at1.P[p].x;
    double dy1 = at2.P[p].y-at1.P[p].y;
    double dy2 = at3.P[p].y-at1.P[p].y;
    double dz1 = at2.P[p].z-at1.P[p].z;
    double dz2 = at3.P[p].z-at1.P[p].z;
    double Theta = dx1*dx2+dy1*dy2+dz1*dz2;
    Theta /= r1*r2;
    Theta = acos(Theta);
    double E = 0.5*ang.K*(Theta-ang.ang)*(Theta-ang.ang);
    //Coulomb subtract, if applicable
    if (CoulSub == 1)
    {
        double r3 = sqrt(Dist2(at2.P[p],at3.P[p]));
        if (r3 == 0.0)
        {
            return HugeNum; //Avoid dividing by zero
        }
        E -= (C2eV*at2.q*at3.q/r3);
    }
    return E;
};
