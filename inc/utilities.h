
#pragma once

#include "system.h"

double Dist2(ParticleCoordinates& a, ParticleCoordinates& b)
{
    //Displacements
    double dx = a.x-b.x;
    double dy = a.y-b.y;
    double dz = a.z-b.z;
    //PBC
    bool check = 1;
    while (check == 1)
    {
        check = 0;
        if (abs(dx) > (0.5*Lx))
        {
            dx = Lx-abs(dx);
            check = 1;
        }
        if (abs(dy) > (0.5*Ly))
        {
            dy = Ly-abs(dy);
            check = 1;
        }
        if (abs(dz) > (0.5*Lz))
        {
            dz = Lz-abs(dz);
            check = 1;
        }
    }
    //Squared radius
    double r2 = dx*dx+dy*dy+dz*dz;
    return r2;
}
