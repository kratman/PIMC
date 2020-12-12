
#pragma once

#include <vector>
#include <string>
#include "coordinates.h"
#include "system.h"

struct BondParameters;
struct AngleParameters;

struct Qpart
{
    std::string typ; //Atom type
    double m; //mass (amu)
    double q; //Charge (au)
    std::vector<double> rad; //Radius, electron only (Ang)
    int spin; //Spin, electron only
    double x; //x position (Ang)
    double y; //y position (Ang)
    double z; //z position (Ang)
    double Ep; //Temp. energy for parallel
    std::vector<ParticleCoordinates> P; //Bead coordinates
    std::vector<BondParameters> Bonds; //Harmonic bonds
    std::vector<AngleParameters> Angs; //Harmonic angles
};

void Get_Centroid(Qpart& part)
{
    double x=0,y=0,z=0;
    for (int i=0;i<Nbeads;i++)
    {
        x += part.P[i].x;
        y += part.P[i].y;
        z += part.P[i].z;
    }
    x /= Nbeads;
    y /= Nbeads;
    z /= Nbeads;
    part.x = x;
    part.y = y;
    part.z = z;
    return;
};
