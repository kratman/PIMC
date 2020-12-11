
#pragma once

const double k = 8.6173324e-5; //Boltzmann constant (eV)
const double hbar = 6.58211928e-16; //Reduced Planck Constant (eV)
const double m2Ang = 1.0e10; //Angstroms to meters
const double amu2kg = 1.660538921e-27; //Atomic mass units to kg
const double pi = 4*atan(1); //Pi
const double SI2eV = 1/(1.602176565e-19); //Convert SI to eV
const double ToeV = amu2kg*SI2eV/(m2Ang*m2Ang); //Convert to eV units
const double C2eV = m2Ang/(4*pi*SI2eV*8.854187817e-12); //Coulomb to eV
const double Masse = 9.10938291e-31; //Mass of an electron (kg)
const double BohrRad = 0.52917721092; //Bohr radius (Ang)
const double Har2eV = 27.21138386; //Hartrees to eV
const double atm2eV = SI2eV*1.01325e-25; //atmA^3 to eV
