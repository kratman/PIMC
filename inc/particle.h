
#pragma once

#include <vector>
#include <string>
#include "coordinates.h"
#include "harmonic_bond.h"
#include "harmonic_angle.h"

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
