
#pragma once

#include "physical_constants.h"

struct LennardJones
{
    double sig; //Sphere radius
    double eps; //Well depth
};

inline double LJEnergy(double s, double e, double r, double q1, double q2)
{
    double E = 0.0;
    if (r == 0.0)
    {
        return HugeNum; //Prevents LJ singularity
    }
    double r6 = s/r;
    r6 *= r6*r6;
    r6 *= r6;
    double r12 = r6*r6;
    //Note: r6 and r12 are inverse
    E = 4*e*(r12-r6)+(C2eV*q1*q2/r);
    return E;
}
