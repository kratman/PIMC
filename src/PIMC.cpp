
#include <ctime>
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>
using namespace std;

//Compile options
const bool Debug = false; //Turn debugging on/off
const bool RCOM = true; //Remove center of mass
const bool CoulSub = true; //Turn Coulomb subtracting on/off
const bool Isotrop = true; //Force isotropic expansion
const bool qGrid = false; //Restrict point charge movements

//Compile time eFF options
const double Rad = 0.10; //Initial electron radius
const double rho = 1.0; //Parameter for eFF VB mixing
const double sbar = 1.0; //Parameter for eFF radius scaling
const double rbar = 1.0; //Parameter for eFF distance scaling
const double Charm = 0; //Force constant for eFF harmonic constraint
const bool Scale_eFF = 0; //Scale the eFF kinetic energy by rtNbeads
const double Scale_POW = 0.05; //rtNbeads = pow(P,Scale_POW)

//Compile-time move options
double step = 0.15; //Initial step size for the beads
const double StepMin = 0.075; //Minimum step size
const double StepMax = 1.0; //Maximum step size
const double Centratio= 10.0; //Scales 'step' for centroids
const int Acc_Check = 5000; //Eq steps before checking accratio
const double RadMin = 0.01; //Minimum electron radius
const double RadMax = 25.0; //Maximum electron radius

//Move Probabilities
//Note: These probabilities allow for multi-particle moves
double BeadProb = 0.55; //Probability to move a single bead
double CentProb = 0.55; //Probability to move a centroid
double ElBeadProb = 0.25; //Probability to move an electron bead
double ElCentProb = 0.25; //Probability to move an electron centroid
double RadProb = 0; //Probability to change electron radius
double SwapProb = 0.05; //Probability to swap spins
double FlipProb = 0.05; //Probability to flip a single spin
double VolProb = 0.05; //Volume change probability

//Physical Constants
const double k = 8.6173324e-5; //Boltzmann constant (eV)
const double hbar = 6.58211928e-16; //Reduced Planck Constant (eV)
const double hbarSI = 1.054571726e-34; //Reduced Planck Constant (SI)
const double kb = 0.69503476; //Boltzmann constant (cm-1)
const double kSI = 1.3806488e-23; //Boltzmann constant (SI)
const double m2Ang = 1.0e10; //Angstroms to meters
const double amu2kg = 1.660538921e-27; //Atomic mass units to kg
const double cs = 2.99792458e8; //Speed of light (m)
const double pi = 4*atan(1); //Pi
const double h = 2*pi*hbar; //Planck Constant (eV)
const double SI2eV = 1/(1.602176565e-19); //Convert SI to eV
const double ToeV = amu2kg*SI2eV/(m2Ang*m2Ang); //Convert to eV units
const double C2eV = m2Ang/(4*pi*SI2eV*8.854187817e-12); //Coulomb to eV
const double Masse = 9.10938291e-31; //Mass of an electron (kg)
const double BohrRad = 0.52917721092; //Bohr radius (Ang)
const double Har2eV = 27.21138386; //Hartrees to eV
const double atm2eV = SI2eV*1.01325e-25; //atmA^3 to eV

//Root constants
const double sqrt2 = sqrt(2); //Square root of 2
const double TwoRtSix = pow(2,(1.0/6.0)); //2**(1/6) for WCA
const double rtMe = 0.2*sqrt(amu2kg/Masse); //Increase electron step

//Global parameters
double Beta,Lx,Ly,Lz,LJcut,Press; //Parameters needed for functions
double rtNbeads; //Scale factor for eFF kinetic energy
int Ensemble = 0; //NVT=0, NPT=1
int Nbeads; //Number of time slices
const double HugeNum = 1e200; //Large number

//Data Structures
struct Coord
{
  double x; //x position
  double y; //y position
  double z; //z position
};

struct BondParam
{
  int at2; //atom ID 2
  double K; //Force constant (eV/Ang^2)
  double R; //Equil. distance (Ang)
};

struct AngParam
{
  int at2; //atom ID 2
  int at3; //atom ID 3
  double K; //Force constant (eV/rad^2)
  double ang; //Equilibrium angle (Rad)
};

struct Qpart
{
  string typ; //Atom type
  double m; //mass (amu)
  double q; //Charge (au)
  vector<double> rad; //Radius, electron only (Ang)
  int spin; //Spin, electron only
  double x; //x position (Ang)
  double y; //y position (Ang)
  double z; //z position (Ang)
  double Ep; //Temp. energy for parallel
  vector<Coord> P; //Bead coordinates
  vector<BondParam> Bonds; //Harmonic bonds
  vector<AngParam> Angs; //Harmonic angles
};

struct LJparam
{
  double sig; //Sphere radius
  double eps; //Well depth
};

double Dist2(Coord& a, Coord& b)
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
};

double LJEnergy(double s, double e, double r, double q1, double q2)
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
};

double SpringEnergy(double k, double r2)
{
  //General harmonic bond for PI rings
  double E = 0.5*k*r2;
  return E;
};

double RadSpring(vector<Qpart>& parts)
{
  //Harmonic bond to keep the radius in check
  double E = 0.0;
  double w0 = 1/(Beta*hbar);
  w0 *= Charm*w0*ToeV;
  #pragma omp parallel for
  for (int i=0;i<parts.size();i++)
  {
    parts[i].Ep = 0.0;
    double w = w0*parts[i].m*Nbeads;
    for (int j=0;j<Nbeads;j++)
    {
      //Weakly restrain the radius to zero
      parts[i].Ep += 0.5*w*parts[i].rad[j]*parts[i].rad[j];
    }
  }
  #pragma omp barrier
  for (int i=0;i<parts.size();i++)
  {
    E += parts[i].Ep;
  }
  return E;
};

double HarmEnergy(BondParam& bond, Qpart& at1, Qpart& at2, int p)
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

double AngEnergy(AngParam& ang, Qpart& at1, Qpart& at2, Qpart& at3, int p)
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

double Get_Espring(vector<Qpart>& parts)
{
  //Calculate total harmonic PI ring energy
  double E = 0.0;
  double w0 = 1/(Beta*hbar);
  w0 *= w0*ToeV;
  #pragma omp parallel for
  for (int i=0;i<parts.size();i++)
  {
    parts[i].Ep = 0.0;
    double w = w0*parts[i].m*Nbeads;
    for (int j=0;j<Nbeads;j++)
    {
      //Bead energy, one bond to avoid double counting
      int j2 = j-1;
      if (j2 == -1)
      {
        j2 = Nbeads-1; //Ring PBC
      }
      double dr2 = Dist2(parts[i].P[j],parts[i].P[j2]);
      parts[i].Ep += SpringEnergy(w,dr2);
    }
  }
  #pragma omp barrier
  for (int i=0;i<parts.size();i++)
  {
    E += parts[i].Ep;
  }
  return E;
};

double Get_Epot(vector<Qpart>& parts, map<string,LJparam>& LJmap)
{
  //Bonded and non-bonded potential for atoms
  double E = 0.0;
  #pragma omp parallel for
  for (int i=0;i<parts.size();i++)
  {
    parts[i].Ep = 0.0;
    for (int j=0;j<Nbeads;j++)
    {
      //Bond energy, Coulomb subtracted
      for (int l=0;l<parts[i].Bonds.size();l++)
      {
        int at2 = parts[i].Bonds[l].at2;
        if (at2 > i)
        {
          parts[i].Ep += HarmEnergy(parts[i].Bonds[l],parts[i],parts[at2],j);
        }
      }
      //Angle energy, Coulomb subtracted
      for (int l=0;l<parts[i].Angs.size();l++)
      {
        int at2 = parts[i].Angs[l].at2;
        int at3 = parts[i].Angs[l].at3;
        parts[i].Ep += AngEnergy(parts[i].Angs[l],parts[i],
        parts[at2],parts[at3],j);
      }
      //Non-Bonded energy
      for (int k=0;k<i;k++)
      {
        double dr2 = Dist2(parts[i].P[j],parts[k].P[j]);
        //LJ and Coulomb energy
        double sig = LJmap[parts[i].typ+parts[k].typ].sig;
        double eps = LJmap[parts[i].typ+parts[k].typ].eps;
        double rtest = LJcut*LJcut;
        if (dr2 < rtest)
        {
          double r = sqrt(dr2);
          parts[i].Ep += LJEnergy(sig,eps,r,parts[i].q,parts[k].q);
        }
      }
    }
  }
  #pragma omp barrier
  for (int i=0;i<parts.size();i++)
  {
    E += parts[i].Ep;
  }
  E /= Nbeads; //Removes double counting
  return E;
};

double EFFEnergy(Qpart& atom, Qpart& elec, int p)
{
  //Atom-lepton interactions
  double E = 0.0;
  double r = Dist2(atom.P[p],elec.P[p]);
  if (r <= LJcut*LJcut)
  {
    if (r == 0.0)
    {
      E = C2eV*atom.q*elec.q*sqrt(8/pi)/elec.rad[p];
    }
    else
    {
      r = sqrt(r);
      E = (C2eV*atom.q*elec.q/r)*erf(sqrt2*r/elec.rad[p]);
    }
  }
  return E;
};

double EFFCorr(Qpart& elec1, Qpart& elec2, int p)
{
  //Lepton-lepton interactions
  double E = 0.0;
  double r = Dist2(elec1.P[p],elec2.P[p]);
  //Electrostatic energy
  if (r <= LJcut*LJcut)
  {
    r = sqrt(r);
    if (r == 0.0)
    {
      E = C2eV*elec1.q*elec2.q*sqrt(8/pi);
      E /= sqrt(elec1.rad[p]*elec1.rad[p]+elec2.rad[p]*elec2.rad[p]);
      return HugeNum; //Escape to avoid singularities later
    }
    else
    {
      double radij = elec1.rad[p]*elec1.rad[p];
      radij += elec2.rad[p]*elec2.rad[p];
      radij = sqrt(radij);
      E = (C2eV*elec1.q*elec2.q/r);
      E *= erf(sqrt2*r/radij);
    }
    //Pauli repulsion
    if (elec1.typ == elec2.typ)
    {
      //Overlap
      double Sij = 2/((elec1.rad[p]/elec2.rad[p])+(elec2.rad[p]/elec1.rad[p]));
      Sij *= Sij*Sij;
      Sij = sqrt(Sij);
      double tmp = -1*rbar*rbar*r*r;
      tmp /= (elec1.rad[p]*elec1.rad[p]+elec2.rad[p]*elec2.rad[p]);
      tmp /= sbar*sbar;
      Sij *= exp(tmp);
      //Kinetic energy difference
      double Tij = 1/(elec1.rad[p]*elec1.rad[p]);
      Tij += 1/(elec2.rad[p]*elec2.rad[p]);
      Tij *= 3/(2*sbar*sbar);
      tmp = 6*sbar*sbar*(elec1.rad[p]*elec1.rad[p]+elec2.rad[p]*elec2.rad[p]);
      tmp -= 4*rbar*rbar*r*r;
      tmp /= sbar*sbar*(elec1.rad[p]*elec1.rad[p]+elec2.rad[p]*elec2.rad[p]);
      tmp /= sbar*sbar*(elec1.rad[p]*elec1.rad[p]+elec2.rad[p]*elec2.rad[p]);
      Tij -= tmp;
      Tij *= Har2eV*BohrRad*BohrRad;
      if (elec1.spin == elec2.spin)
      {
        //Symmetric VB spin-orbital
        double Etmp = Sij*Sij/(1-(Sij*Sij));
        Etmp += (1-rho)*Sij*Sij/(1+(Sij*Sij));
        E += Etmp*Tij;
      }
      else
      {
        //Antisymmetric VB spin orbital
        E += -1*rho*Sij*Sij*Tij/(1+(Sij*Sij));
      }
    }
  }
  return E;
};

double Get_EeFF(vector<Qpart>& parts, vector<Qpart>& elecs)
{
  //Total eFF interaction energy
  double E = 0;
  #pragma omp parallel for
  for (int i=0;i<parts.size();i++)
  {
    parts[i].Ep = 0;
    for (int j=0;j<elecs.size();j++)
    {
      for (int k=0;k<Nbeads;k++)
      {
        parts[i].Ep += EFFEnergy(parts[i],elecs[j],k);
      }
    }
  }
  #pragma omp barrier
  #pragma omp parallel for
  for (int i=0;i<elecs.size();i++)
  {
    elecs[i].Ep = 0.0;
    for (int k=0;k<Nbeads;k++)
    {
      //Lepton kinetic energy
      double Etmp = 3/(2*elecs[i].rad[k]*elecs[i].rad[k]);
      Etmp *= Har2eV*BohrRad*BohrRad;
      if (Scale_eFF == 1)
      {
        //Reduce kinetic energy as the beads increase
        Etmp /= rtNbeads;
      }
      elecs[i].Ep += Etmp;
      for (int j=0;j<i;j++)
      {
        //Lepton-lepton correlation
        elecs[i].Ep += EFFCorr(elecs[i],elecs[j],k);
      }
    }
  }
  #pragma omp barrier
  for (int i=0;i<parts.size();i++)
  {
    E += parts[i].Ep;
  }
  for (int i=0;i<elecs.size();i++)
  {
    E += elecs[i].Ep;
  }
  E /= Nbeads; //Removes double counting
  return E;
};

bool MCMove(vector<Qpart>& parts, vector<Qpart>& elecs,
 map<string,LJparam>& LJmap)
{
  bool acc = 0;
  //Copy parts
  vector<Qpart> parts2;
  vector<Qpart> elecs2;
  parts2 = parts;
  elecs2 = elecs;
  //Pick random move and apply PBC
  double randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum > (1-CentProb))
  {
    //Move a centroid
    int p = (rand()%parts2.size());
    double randx = (((double)rand())/((double)RAND_MAX));
    double randy = (((double)rand())/((double)RAND_MAX));
    double randz = (((double)rand())/((double)RAND_MAX));
    double dx = 2*(randx-0.5)*step*Centratio;
    double dy = 2*(randy-0.5)*step*Centratio;
    double dz = 2*(randz-0.5)*step*Centratio;
    double x = parts2[p].x+dx;
    double y = parts2[p].y+dy;
    double z = parts2[p].z+dz;
    bool check = 1;
    while (check == 1)
    {
      check = 0;
      if (x > Lx)
      {
        x -= Lx;
        check = 1;
      }
      if (x < 0.0)
      {
        x += Lx;
        check = 1;
      }
      if (y > Ly)
      {
        y -= Ly;
        check = 1;
      }
      if (y < 0.0)
      {
        y += Ly;
        check = 1;
      }
      if (z > Lz)
      {
        z -= Lz;
        check = 1;
      }
      if (z < 0.0)
      {
        z += Lz;
        check = 1;
      }
    }
    parts2[p].x = x;
    parts2[p].y = y;
    parts2[p].z = z;
    #pragma omp parallel for
    for (int i=0;i<Nbeads;i++)
    {
      double xp = parts2[p].P[i].x+dx;
      double yp = parts2[p].P[i].y+dy;
      double zp = parts2[p].P[i].z+dz;
      bool check = 1;
      while (check == 1)
      {
        check = 0;
        if (xp > Lx)
        {
          xp -= Lx;
          check = 1;
        }
        if (xp < 0.0)
        {
          xp += Lx;
          check = 1;
        }
        if (yp > Ly)
        {
          yp -= Ly;
          check = 1;
        }
        if (yp < 0.0)
        {
          yp += Ly;
          check = 1;
        }
        if (zp > Lz)
        {
          zp -= Lz;
          check = 1;
        }
        if (zp < 0.0)
        {
          zp += Lz;
          check = 1;
        }
      }
      parts2[p].P[i].x = xp;
      parts2[p].P[i].y = yp;
      parts2[p].P[i].z = zp;
    }
    #pragma omp barrier
  }
  if (randnum < BeadProb)
  {
    //Move a single bead
    int p = (rand()%parts2.size());
    int p2 = (rand()%Nbeads);
    double randx = (((double)rand())/((double)RAND_MAX));
    double randy = (((double)rand())/((double)RAND_MAX));
    double randz = (((double)rand())/((double)RAND_MAX));
    double dx = 2*(randx-0.5)*step;
    double dy = 2*(randy-0.5)*step;
    double dz = 2*(randz-0.5)*step;
    parts2[p].P[p2].x += dx;
    bool check = 1;
    while (check == 1)
    {
      check = 0;
      if (parts2[p].P[p2].x > Lx)
      {
        parts2[p].P[p2].x -= Lx;
        check = 1;
      }
      if (parts2[p].P[p2].x < 0.0)
      {
        parts2[p].P[p2].x += Lx;
        check = 1;
      }
    }
    parts2[p].P[p2].y += dy;
    check = 1;
    while (check == 1)
    {
      check = 0;
      if (parts2[p].P[p2].y > Ly)
      {
        parts2[p].P[p2].y -= Ly;
        check = 1;
      }
      if (parts2[p].P[p2].y < 0.0)
      {
        parts2[p].P[p2].y += Ly;
        check = 1;
      }
    }
    parts2[p].P[p2].z += dz;
    check = 1;
    while (check == 1)
    {
      check = 0;
      if (parts2[p].P[p2].z > Lz)
      {
        parts2[p].P[p2].z -= Lz;
        check = 1;
      }
      if (parts2[p].P[p2].z < 0.0)
      {
        parts2[p].P[p2].z += Lz;
        check = 1;
      }
    }
  }
  //Pick random electron move and apply PBC
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum < ElCentProb)
  {
    //Move an electron centroid
    int p = (rand()%elecs2.size());
    double randx = (((double)rand())/((double)RAND_MAX));
    double randy = (((double)rand())/((double)RAND_MAX));
    double randz = (((double)rand())/((double)RAND_MAX));
    double dx = 2*(randx-0.5)*step*Centratio;
    double dy = 2*(randy-0.5)*step*Centratio;
    double dz = 2*(randz-0.5)*step*Centratio;
    double x = elecs2[p].x+dx;
    double y = elecs2[p].y+dy;
    double z = elecs2[p].z+dz;
    bool check = 1;
    while (check == 1)
    {
      check = 0;
      if (x > Lx)
      {
        x -= Lx;
        check = 1;
      }
      if (x < 0.0)
      {
        x += Lx;
        check = 1;
      }
      if (y > Ly)
      {
        y -= Ly;
        check = 1;
      }
      if (y < 0.0)
      {
        y += Ly;
        check = 1;
      }
      if (z > Lz)
      {
        z -= Lz;
        check = 1;
      }
      if (z < 0.0)
      {
        z += Lz;
        check = 1;
      }
    }
    elecs2[p].x = x;
    elecs2[p].y = y;
    elecs2[p].z = z;
    #pragma omp parallel for
    for (int i=0;i<Nbeads;i++)
    {
      double xp = elecs2[p].P[i].x+dx;
      double yp = elecs2[p].P[i].y+dy;
      double zp = elecs2[p].P[i].z+dz;
      bool check = 1;
      while (check == 1)
      {
        check = 0;
        if (xp > Lx)
        {
          xp -= Lx;
          check = 1;
        }
        if (xp < 0.0)
        {
          xp += Lx;
          check = 1;
        }
        if (yp > Ly)
        {
          yp -= Ly;
          check = 1;
        }
        if (yp < 0.0)
        {
          yp += Ly;
          check = 1;
        }
        if (zp > Lz)
        {
          zp -= Lz;
          check = 1;
        }
        if (zp < 0.0)
        {
          zp += Lz;
          check = 1;
        }
      }
      elecs2[p].P[i].x = xp;
      elecs2[p].P[i].y = yp;
      elecs2[p].P[i].z = zp;
    }
    #pragma omp barrier
  }
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum < ElBeadProb)
  {
    //Move a single electron bead
    int p = (rand()%elecs2.size());
    int p2 = (rand()%Nbeads);
    double randx = (((double)rand())/((double)RAND_MAX));
    double randy = (((double)rand())/((double)RAND_MAX));
    double randz = (((double)rand())/((double)RAND_MAX));
    double dx = 2*(randx-0.5)*step;
    double dy = 2*(randy-0.5)*step;
    double dz = 2*(randz-0.5)*step;
    elecs2[p].P[p2].x += dx;
    elecs2[p].P[p2].y += dy;
    elecs2[p].P[p2].z += dz;
    bool check = 1;
    while (check == 1)
    {
      check = 0;
      if (elecs2[p].P[p2].x > Lx)
      {
        elecs2[p].P[p2].x -= Lx;
        check = 1;
      }
      if (elecs2[p].P[p2].x < 0.0)
      {
        elecs2[p].P[p2].x += Lx;
        check = 1;
      }
      if (elecs2[p].P[p2].y > Ly)
      {
        elecs2[p].P[p2].y -= Ly;
        check = 1;
      }
      if (elecs2[p].P[p2].y < 0.0)
      {
        elecs2[p].P[p2].y += Ly;
        check = 1;
      }
      if (elecs2[p].P[p2].z > Lz)
      {
        elecs2[p].P[p2].z -= Lz;
        check = 1;
      }
      if (elecs2[p].P[p2].z < 0.0)
      {
        elecs2[p].P[p2].z += Lz;
        check = 1;
      }
    }
  }
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum < RadProb)
  {
    //Change electron radius
    int p = (rand()%elecs2.size());
    int p2 = (rand()%Nbeads);
    randnum = (((double)rand())/((double)RAND_MAX));
    double radnew = elecs2[p].rad[p2];
    radnew *= 0.90+(randnum*0.20);
    if (radnew > RadMax)
    {
      radnew = RadMax;
    }
    if (radnew < RadMin)
    {
      radnew = RadMin;
    }
    elecs2[p].rad[p2] = radnew;
  }
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum < SwapProb)
  {
    int p = (rand()%elecs2.size());
    int p2 = (rand()%elecs2.size());
    int tmp = elecs2[p].spin;
    int tmp2 = elecs2[p2].spin;
    elecs2[p].spin = tmp2;
    elecs2[p2].spin = tmp;
  }
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum < FlipProb)
  {
    int p = (rand()%elecs2.size());
    elecs2[p].spin *= -1;
  }
  //Calculate energies
  double Eold = 0;
  if (parts.size() > 0)
  {
    Eold += Get_Epot(parts,LJmap);
    Eold += Get_Espring(parts);
  }
  if (elecs.size() > 0)
  {
    Eold += Get_EeFF(parts,elecs);
    Eold += Get_Espring(elecs)+RadSpring(elecs);
  }
  double Enew = 0;
  double Lxtmp = Lx;
  double Lytmp = Ly;
  double Lztmp = Lz;
  randnum = (((double)rand())/((double)RAND_MAX));
  if (randnum < VolProb)
  {
    if (Isotrop == 0)
    {
      randnum = (((double)rand())/((double)RAND_MAX));
      double Lmin = 0.90*Lx;
      double Lmax = 1.10*Lx;
      Lx = Lmin+randnum*(Lmax-Lmin);
      randnum = (((double)rand())/((double)RAND_MAX));
      Lmin = 0.90*Ly;
      Lmax = 1.10*Ly;
      Ly = Lmin+randnum*(Lmax-Lmin);
      randnum = (((double)rand())/((double)RAND_MAX));
      Lmin = 0.90*Lz;
      Lmax = 1.10*Lz;
      Lz = Lmin+randnum*(Lmax-Lmin);
      if (Lx < 2.1*LJcut)
      {
        Lx = 2.1*LJcut;
      }
      if (Ly < 2.1*LJcut)
      {
        Ly = 2.1*LJcut;
      }
      if (Lz < 2.1*LJcut)
      {
        Lz = 2.1*LJcut;
      }
    }
    if (Isotrop == 1)
    {
      randnum = (((double)rand())/((double)RAND_MAX));
      double expan = 0.9+randnum*0.20;
      Lx *= expan;
      Ly *= expan;
      Lz *= expan;
      bool bad_size = 0;
      if (Lx < 2.1*LJcut)
      {
        bad_size = 1;
      }
      if (Ly < 2.1*LJcut)
      {
        bad_size = 1;
      }
      if (Lz < 2.1*LJcut)
      {
        bad_size = 1;
      }
      if (bad_size == 1)
      {
        Lx = Lxtmp;
        Ly = Lytmp;
        Lz = Lztmp;
      }
    }
    #pragma omp parallel for
    for (int i=0;i<parts2.size();i++)
    {
      for (int j=0;j<Nbeads;j++)
      {
        parts2[i].P[j].x *= Lx/Lxtmp;
        parts2[i].P[j].y *= Ly/Lytmp;
        parts2[i].P[j].z *= Lz/Lztmp;
      }
      parts2[i].x *= Lx/Lxtmp;
      parts2[i].y *= Ly/Lytmp;
      parts2[i].z *= Lz/Lztmp;
    }
    #pragma omp barrier
    #pragma omp parallel for
    for (int i=0;i<elecs2.size();i++)
    {
      for (int j=0;j<Nbeads;j++)
      {
        elecs2[i].P[j].x *= Lx/Lxtmp;
        elecs2[i].P[j].y *= Ly/Lytmp;
        elecs2[i].P[j].z *= Lz/Lztmp;
      }
      elecs2[i].x *= Lx/Lxtmp;
      elecs2[i].y *= Ly/Lytmp;
      elecs2[i].z *= Lz/Lztmp;
    }
    #pragma omp barrier
    Eold += Press*Lxtmp*Lytmp*Lztmp*atm2eV;
    Enew += Press*Lx*Ly*Lz*atm2eV;
  }
  if (parts.size() > 0)
  {
    Enew += Get_Epot(parts2,LJmap);
    Enew += Get_Espring(parts2);
  }
  if (elecs.size() > 0)
  {
    Enew += Get_EeFF(parts2,elecs2);
    Enew += Get_Espring(elecs2)+RadSpring(elecs2);
  }
  //Accept or reject
  double dE = Enew-Eold;
  if ((dE == 0.0) and (Debug == 1))
  {
    cout << "Warning two structures have identical energies." << '\n';
  }
  double Prob = exp(-1*dE*Beta);
  randnum = (((double)rand())/((double)RAND_MAX));
  if ((dE <= 0) or (randnum < Prob))
  {
    parts = parts2;
    elecs = elecs2;
    acc = 1;
  }
  else
  {
    Lx = Lxtmp;
    Ly = Lytmp;
    Lz = Lztmp;
  }
  return acc;
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

Coord Get_COM(vector<Qpart>& parts, vector<Qpart>& elecs)
{
  double x=0,y=0,z=0,M=0;
  Coord com;
  //Gather all atom centroids
  #pragma omp parallel for
  for (int i=0;i<parts.size();i++)
  {
    Get_Centroid(parts[i]);
  }
  #pragma omp barrier
  //Gather all electron centroids
  #pragma omp parallel for
  for (int i=0;i<elecs.size();i++)
  {
    Get_Centroid(elecs[i]);
  }
  #pragma omp barrier
  //Calculate COM
  for (int i=0;i<parts.size();i++)
  {
    //Atoms
    x += parts[i].m*parts[i].x;
    y += parts[i].m*parts[i].y;
    z += parts[i].m*parts[i].z;
    M += parts[i].m;
  }
  for (int i=0;i<elecs.size();i++)
  {
    //Electrons
    x += elecs[i].m*elecs[i].x;
    y += elecs[i].m*elecs[i].y;
    z += elecs[i].m*elecs[i].z;
    M += elecs[i].m;
  }
  if ((parts.size() > 0) or (elecs.size() > 0))
  {
    com.x = x/M;
    com.y = y/M;
    com.z = z/M;
  }
  else
  {
    com.x = Lx/2;
    com.y = Ly/2;
    com.z = Lz/2;
  }
  return com;
};

void Get_Range(vector<Qpart>& parts, Coord& Ls, Coord& Hs)
{
  //Find the min and max values of the positions
  double minx = 1000000000.0;
  double miny = 1000000000.0;
  double minz = 1000000000.0;
  double maxx = -1000000000.0;
  double maxy = -1000000000.0;
  double maxz = -1000000000.0;
  for (int i=0;i<parts.size();i++)
  {
    Get_Centroid(parts[i]);
    if (parts[i].x < minx)
    {
      minx = parts[i].x;
    }
    if (parts[i].x > maxx)
    {
      maxx = parts[i].x;
    }
    if (parts[i].y < miny)
    {
      miny = parts[i].y;
    }
    if (parts[i].y > maxy)
    {
      maxy = parts[i].y;
    }
    if (parts[i].z < minz)
    {
      minz = parts[i].z;
    }
    if (parts[i].z > maxz)
    {
      maxz = parts[i].z;
    }
  }
  Ls.x = minx;
  Ls.y = miny;
  Ls.z = minz;
  Hs.x = maxx;
  Hs.y = maxy;
  Hs.z = maxz;
  if (parts.size() == 0)
  {
    Ls.x = (Lx/2)-(LJcut/2);
    Ls.y = (Ly/2)-(LJcut/2);
    Ls.z = (Lz/2)-(LJcut/2);
    Hs.x = (Lx/2)+(LJcut/2);
    Hs.y = (Ly/2)+(LJcut/2);
    Hs.z = (Lz/2)+(LJcut/2);
  }
  return;
};

void Remove_COM(vector<Qpart>& parts, vector<Qpart>& elecs)
{
  Coord com = Get_COM(parts,elecs);
  //Subtract COM for atoms
  #pragma omp parallel for
  for (int i=0;i<parts.size();i++)
  {
    double x = parts[i].x-com.x-0.5*Lx;
    double y = parts[i].y-com.y-0.5*Ly;
    double z = parts[i].z-com.z-0.5*Lz;
    bool check = 1;
    while (check == 1)
    {
      check = 0;
      if (x > Lx)
      {
        x -= Lx;
        check = 1;
      }
      if (x < 0.0)
      {
        x += Lx;
        check = 1;
      }
      if (y > Ly)
      {
        y -= Ly;
        check = 1;
      }
      if (y < 0.0)
      {
        y += Ly;
        check = 1;
      }
      if (z > Lz)
      {
        z -= Lz;
        check = 1;
      }
      if (z < 0.0)
      {
        z += Lz;
        check = 1;
      }
    }
    parts[i].x = x;
    parts[i].y = y;
    parts[i].z = z;
    for (int j=0;j<Nbeads;j++)
    {
      x = parts[i].P[j].x-com.x-0.5*Lx;
      y = parts[i].P[j].y-com.y-0.5*Ly;
      z = parts[i].P[j].z-com.z-0.5*Lz;
      check = 1;
      while (check == 1)
      {
        check = 0;
        if (x > Lx)
        {
          x -= Lx;
          check = 1;
        }
        if (x < 0.0)
        {
          x += Lx;
          check = 1;
        }
        if (y > Ly)
        {
          y -= Ly;
          check = 1;
        }
        if (y < 0.0)
        {
          y += Ly;
          check = 1;
        }
        if (z > Lz)
        {
          z -= Lz;
          check = 1;
        }
        if (z < 0.0)
        {
          z += Lz;
          check = 1;
        }
      }
      parts[i].P[j].x = x;
      parts[i].P[j].y = y;
      parts[i].P[j].z = z;
    }
  }
  #pragma omp barrier
  //Repeat for electrons
  #pragma omp parallel for
  for (int i=0;i<elecs.size();i++)
  {
    double x = elecs[i].x-com.x-0.5*Lx;
    double y = elecs[i].y-com.y-0.5*Ly;
    double z = elecs[i].z-com.z-0.5*Lz;
    bool check = 1;
    while (check == 1)
    {
      check = 0;
      if (x > Lx)
      {
        x -= Lx;
        check = 1;
      }
      if (x < 0.0)
      {
        x += Lx;
        check = 1;
      }
      if (y > Ly)
      {
        y -= Ly;
        check = 1;
      }
      if (y < 0.0)
      {
        y += Ly;
        check = 1;
      }
      if (z > Lz)
      {
        z -= Lz;
        check = 1;
      }
      if (z < 0.0)
      {
        z += Lz;
        check = 1;
      }
    }
    elecs[i].x = x;
    elecs[i].y = y;
    elecs[i].z = z;
    for (int j=0;j<Nbeads;j++)
    {
      x = elecs[i].P[j].x-com.x-0.5*Lx;
      y = elecs[i].P[j].y-com.y-0.5*Ly;
      z = elecs[i].P[j].z-com.z-0.5*Lz;
      check = 1;
      while (check == 1)
      {
        check = 0;
        if (x > Lx)
        {
          x -= Lx;
          check = 1;
        }
        if (x < 0.0)
        {
          x += Lx;
          check = 1;
        }
        if (y > Ly)
        {
          y -= Ly;
          check = 1;
        }
        if (y < 0.0)
        {
          y += Ly;
          check = 1;
        }
        if (z > Lz)
        {
          z -= Lz;
          check = 1;
        }
        if (z < 0.0)
        {
          z += Lz;
          check = 1;
        }
      }
      elecs[i].P[j].x = x;
      elecs[i].P[j].y = y;
      elecs[i].P[j].z = z;
    }
  }
  #pragma omp barrier
  return;
};

void Print_traj(vector<Qpart>& parts, vector<Qpart>& elecs,
 fstream& traj, string mode)
{
  if (RCOM == 1)
  {
    Remove_COM(parts,elecs);
  }
  if ((mode == "All") or (mode == "all"))
  {
    int Ntot = 0;
    if (parts.size() > 0)
    {
      Ntot += Nbeads*parts.size();
    }
    if (elecs.size() > 0)
    {
      Ntot += Nbeads*elecs.size();
    }
    traj << Ntot << '\n' << '\n';
    for (int i=0;i<parts.size();i++)
    {
      for (int j=0;j<Nbeads;j++)
      {
        traj << parts[i].typ << " ";
        traj << parts[i].P[j].x << " ";
        traj << parts[i].P[j].y << " ";
        traj << parts[i].P[j].z << '\n';
      }
    }
    for (int i=0;i<elecs.size();i++)
    {
      for (int j=0;j<Nbeads;j++)
      {
        traj << elecs[i].typ << " ";
        traj << elecs[i].P[j].x << " ";
        traj << elecs[i].P[j].y << " ";
        traj << elecs[i].P[j].z << '\n';
      }
    }
  }
  if ((mode == "COM") or (mode == "com"))
  {
    traj << parts.size()+elecs.size() << '\n' << '\n';
    for (int i=0;i<parts.size();i++)
    {
      traj << parts[i].typ << " ";
      traj << parts[i].x << " ";
      traj << parts[i].y << " ";
      traj << parts[i].z << '\n';
    }
    for (int i=0;i<elecs.size();i++)
    {
      traj << elecs[i].typ << " ";
      traj << elecs[i].x << " ";
      traj << elecs[i].y << " ";
      traj << elecs[i].z << '\n';
    }
  }
  traj.flush();
  return;
};


//##  Main Code ##//
int main()
{
  cout << '\n';
  //Initialize dynamic parameters
  srand((unsigned)time(0)); //Serial only random numbers
  string dummy,PrintMode,SpinMode;
  string filename1 = "pimc.param"; //Input file
  string filename2 = "traj.xyz"; //Output trajectory
  string filename3 = "pimc.pot"; //Potential file
  string filename4; //xyz input
  string filename5; //Connectivity input
  fstream paramfile,trajfile,potfile,xyzfile,bondfile;
  int Natoms,Nsteps,Ntyps,Neq,Nprint,Npre,Nelec,Npos;
  double accratio,SumE,SumE2,VolAvg,Ek;
  vector<Qpart> Atoms;
  vector<string> Types;
  map<string,LJparam> LJparams;
  paramfile.open(filename1.c_str(),ios_base::in);
  trajfile.open(filename2.c_str(),ios_base::out);
  potfile.open(filename3.c_str(),ios_base::in);
  //End of section

  //Read Input
  if (Debug == 1)
  {
    cout << "Reading input..." << '\n';
  }
  paramfile >> dummy;
  paramfile >> dummy >> dummy;
  if (dummy == "NVT")
  {
    Ensemble = 0;
  }
  if (dummy == "NPT")
  {
    Ensemble = 1;
  }
  paramfile >> dummy >> Beta;
  paramfile >> dummy >> Press;
  paramfile >> dummy >> Neq;
  paramfile >> dummy >> Nsteps;
  paramfile >> dummy >> Nbeads;
  paramfile >> dummy >> accratio;
  paramfile >> dummy >> Nprint;
  paramfile >> dummy >> PrintMode;
  paramfile >> dummy >> filename4;
  paramfile >> dummy >> Lx >> Ly >> Lz;
  paramfile >> dummy >> LJcut;
  paramfile >> dummy >> filename5;
  paramfile >> dummy >> SpinMode;
  paramfile >> dummy >> Npre;
  paramfile >> dummy >> Nelec;
  paramfile >> dummy >> Npos;
  xyzfile.open(filename4.c_str(),ios_base::in);
  xyzfile >> Natoms;
  for (int i=0;i<Natoms;i++)
  {
    string typ;
    double x,y,z;
    xyzfile >> typ >> x >> y >> z;
    Qpart tmp;
    tmp.typ = typ;
    tmp.x = x;
    tmp.y = y;
    tmp.z = z;
    tmp.spin = 0;
    for (int j=0;j<Nbeads;j++)
    {
      tmp.rad.push_back(0);
      double randx = (((double)rand())/((double)RAND_MAX));
      double randy = (((double)rand())/((double)RAND_MAX));
      double randz = (((double)rand())/((double)RAND_MAX));
      if (j == 0)
      {
        randx = 0.5;
        randy = 0.5;
        randz = 0.5;
      }
      if (qGrid == 1)
      {
        randx = 0.5;
        randy = 0.5;
        randz = 0.5;
      }
      Coord temp;
      //Set random bead displacements
      temp.x = x+(randx-0.5)*step;
      temp.y = y+(randy-0.5)*step;
      temp.z = z+(randz-0.5)*step;
      tmp.P.push_back(temp);
    }
    Get_Centroid(tmp);
    Atoms.push_back(tmp);
  }
  potfile >> dummy;
  potfile >> dummy >> Ntyps;
  potfile >> dummy; //Types
  for (int i=0;i<Ntyps;i++)
  {
    string tmp;
    potfile >> tmp;
    Types.push_back(tmp);
  }
  potfile >> dummy; //Masses
  for (int i=0;i<Ntyps;i++)
  {
    double tmp;
    potfile >> tmp;
    for (int j=0;j<Natoms;j++)
    {
      if (Types[i] == Atoms[j].typ)
      {
        Atoms[j].m = tmp;
      }
    }
  }
  potfile >> dummy; //Charges
  for (int i=0;i<Ntyps;i++)
  {
    double tmp;
    potfile >> tmp;
    for (int j=0;j<Natoms;j++)
    {
      if (Types[i] == Atoms[j].typ)
      {
        Atoms[j].q = tmp;
      }
    }
  }
  potfile >> dummy; //Lennard-Jones
  for (int i=0;i<Ntyps;i++)
  {
    double sig,eps;
    potfile >> eps >> sig;
    LJparam tmp;
    tmp.sig = sig;
    tmp.eps = eps;
    LJparams[Types[i]+Types[i]] = tmp;
  }
  for (int i=0;i<Ntyps;i++)
  {
    for (int j=0;j<Ntyps;j++)
    {
      if (i != j)
      {
        double sig,eps;
        sig = LJparams[Types[i]+Types[i]].sig;
        sig += LJparams[Types[j]+Types[j]].sig;
        sig *= 0.5;
        eps = LJparams[Types[i]+Types[i]].eps;
        eps *= LJparams[Types[j]+Types[j]].eps;
        eps = sqrt(eps);
        LJparam tmp;
        tmp.sig = sig;
        tmp.eps = eps;
        LJparams[Types[i]+Types[j]] = tmp;
        LJparams[Types[j]+Types[i]] = tmp;
      }
    }
  }
  if (filename5 != "None")
  {
    if (filename5 != "none")
    {
      //Collect bonding info
      bondfile.open(filename5.c_str(),ios_base::in);
      int Nbonds;
      bondfile >> dummy >> Nbonds;
      for (int i=0;i<Nbonds;i++)
      {
        BondParam tmp,tmp2;
        bondfile >> tmp.at2 >> tmp2.at2 >> tmp.K >> tmp.R;
        tmp2.K = tmp.K;
        tmp2.R = tmp.R;
        if ((tmp.at2 > Atoms.size()) or (tmp.at2 > Atoms.size()))
        {
          cout << "Bond to non-existent atom!!!!" << '\n';
          cout << "Check connectivity file." << '\n';
          return 0;
        }
        Atoms[tmp2.at2].Bonds.push_back(tmp);
        Atoms[tmp.at2].Bonds.push_back(tmp2);
      }
      //Collect angle info
      bondfile >> dummy >> Nbonds;
      for (int i=0;i<Nbonds;i++)
      {
        AngParam tmp;
        int cent;
        bondfile >> cent >> tmp.at2 >> tmp.at3 >> tmp.K >> tmp.ang;
        tmp.ang *= pi/180.0;
        if ((tmp.at2 > Atoms.size()) or (tmp.at2 > Atoms.size())
        or (cent > Atoms.size()))
        {
          cout << "Angle with non-existent atom!!!!" << '\n';
          cout << "Check connectivity file." << '\n';
          return 0;
        }
        Atoms[cent].Angs.push_back(tmp);
      }
    }
  }
  vector<Qpart> Elecs;
  if (Nelec > 0)
  {
    int spinct = -1;
    Coord Lows;
    Coord Highs;
    Get_Range(Atoms,Lows,Highs);
    //Create electron
    for (int i=0;i<Nelec;i++)
    {
      spinct *= -1; //Make spins alternate
      Qpart tmp;
      tmp.spin = spinct;
      tmp.m = Masse/amu2kg;
      tmp.typ = "e";
      tmp.q = -1;
      tmp.Ep = 0;
      double randx = (((double)rand())/((double)RAND_MAX));
      double randy = (((double)rand())/((double)RAND_MAX));
      double randz = (((double)rand())/((double)RAND_MAX));
      tmp.x = Lows.x+randx*(Highs.x-Lows.x);
      tmp.y = Lows.y+randy*(Highs.y-Lows.y);
      tmp.z = Lows.z+randz*(Highs.z-Lows.z);
      //Set random bead displacements
      for (int j=0;j<Nbeads;j++)
      {
        tmp.rad.push_back(Rad);
        Coord Temp;
        randx = (((double)rand())/((double)RAND_MAX));
        randy = (((double)rand())/((double)RAND_MAX));
        randz = (((double)rand())/((double)RAND_MAX));
        if (j == 0)
        {
          randx = 0.5;
          randy = 0.5;
          randz = 0.5;
        }
        Temp.x = tmp.x+(randx-0.5)*step;
        Temp.y = tmp.y+(randy-0.5)*step;
        Temp.z = tmp.z+(randz-0.5)*step;
        tmp.P.push_back(Temp);
      }
      Get_Centroid(tmp);
      Elecs.push_back(tmp);
    }
  }
  if (Npos > 0)
  {
    int spinct = -1;
    Coord Lows;
    Coord Highs;
    Get_Range(Atoms,Lows,Highs);
    //Create positron
    for (int i=0;i<Npos;i++)
    {
      spinct *= -1; //Make spins alternate
      Qpart tmp;
      tmp.spin = spinct;
      tmp.m = Masse/amu2kg;
      tmp.typ = "p";
      tmp.q = 1;
      tmp.Ep = 0;
      double randx = (((double)rand())/((double)RAND_MAX));
      double randy = (((double)rand())/((double)RAND_MAX));
      double randz = (((double)rand())/((double)RAND_MAX));
      tmp.x = Lows.x+randx*(Highs.x-Lows.x);
      tmp.y = Lows.y+randy*(Highs.y-Lows.y);
      tmp.z = Lows.z+randz*(Highs.z-Lows.z);
      //Set random bead displacements
      for (int j=0;j<Nbeads;j++)
      {
        tmp.rad.push_back(Rad);
        Coord Temp;
        randx = (((double)rand())/((double)RAND_MAX));
        randy = (((double)rand())/((double)RAND_MAX));
        randz = (((double)rand())/((double)RAND_MAX));
        if (j == 0)
        {
          randx = 0.5;
          randy = 0.5;
          randz = 0.5;
        }
        Temp.x = tmp.x+(randx-0.5)*step;
        Temp.y = tmp.y+(randy-0.5)*step;
        Temp.z = tmp.z+(randz-0.5)*step;
        tmp.P.push_back(Temp);
      }
      Get_Centroid(tmp);
      Elecs.push_back(tmp);
    }
  }
  //End of section

  //Print input for error checking
  cout << "Setting up simulation..." << '\n';
  cout << '\n';
  cout << "Atoms: " << Natoms << '\n';
  cout << "Electrons: " << Nelec << '\n';
  cout << "Positrons: " << Npos << '\n';
  cout << "Beads: " << Nbeads << '\n';
  cout << '\n';
  if ((Npre > 0) and (Elecs.size() > 0))
  {
    cout << "Pre-equilibration steps for electrons: ";
    cout << Npre << '\n';
  }
  cout << "Equilibration steps: " << Neq << '\n';
  cout << "Steps for production run: " << Nsteps << '\n';

  //Adjust probabilities
  if (Atoms.size() == 0)
  {
    //Remove atom moves
    BeadProb = 0.0;
    CentProb = 0.0;
    ElCentProb = ElCentProb/(ElCentProb+ElBeadProb);
    ElBeadProb = 1.0;
  }
  if (Elecs.size() == 0)
  {
    //Remove electron moves
    ElCentProb = 0.0;
    ElBeadProb = 0.0;
    RadProb = 0.0;
    FlipProb = 0.0;
    SwapProb = 0.0;
  }
  if ((Atoms.size() == 0) and (Elecs.size() == 1))
  {
    //Remove electron centroid moves
    ElCentProb = 0.0;
    ElBeadProb = 1.0;
  }
  if ((Atoms.size() == 1) and (Elecs.size() == 0))
  {
    //Remove atom centroid moves
    CentProb = 0.0;
    BeadProb = 1.0;
  }
  if (Elecs.size() <= 2)
  {
    SwapProb = 0.0;
  }
  if (Elecs.size() <= 1)
  {
    FlipProb = 0.0;
  }
  if ((SpinMode == "LPIMC") or (SpinMode == "lpimc")
  or (SpinMode == "LOW") or (SpinMode == "low"))
  {
    FlipProb = 0;
    cout << "Mode: ";
    if (Ensemble == 0)
    {
      cout << "NVT";
    }
    if (Ensemble == 1)
    {
      cout << "NPT";
    }
    if (Elecs.size() > 0)
    {
      cout << " Low-Spin";
    }
    cout << " PIMC" << '\n';
  }
  else
  {
    cout << "Mode: ";
    if (Ensemble == 0)
    {
      cout << "NVT";
    }
    if (Ensemble == 1)
    {
      cout << "NPT";
    }
    if (Elecs.size() > 0)
    {
      cout << " High-Spin";
    }
    cout << " PIMC" << '\n';
  }
  if (Ensemble == 0)
  {
    VolProb = 0.0;
  }
  if (qGrid == 1)
  {
    CentProb = 0.0;
    BeadProb = 0.0;
  }

  //Run simulations
  cout << '\n';
  Beta = 1/(k*Beta); //Invert temperature
  rtNbeads = pow(Nbeads,Scale_POW);
  SumE = 0;
  SumE2 = 0;
  VolAvg = 0;
  Ek = 0.0;
  if ((Atoms.size() > 0) and (qGrid == 0))
  {
    //Add atom kinetic energy
    Ek += 3*Atoms.size()*Nbeads/(2*Beta);
  }
  if (Elecs.size() > 0)
  {
    Ek += 3*Elecs.size()*Nbeads/(2*Beta);
  }
  int Nct = 0; //Step counter
  int ct = 0; //Secondary counter
  double Nacc = 0;
  double Nrej = 0;
  bool acc;
  if (RCOM == 1)
  {
    Remove_COM(Atoms,Elecs);
  }
  if ((Npre > 0) and (Elecs.size() > 0))
  {
    cout << "Starting electron pre-equilibration..." << '\n';
    double tmp1 = BeadProb;
    if (Atoms.size() > 0)
    {
      BeadProb = 0.10;
    }
    if (qGrid == 1)
    {
      BeadProb = 0;
    }
    double tmp2 = CentProb;
    CentProb = 0.0;
    double tmp3 = ElBeadProb;
    ElBeadProb = 1.0;
    double tmp4 = ElCentProb;
    if ((Atoms.size() > 0) and (Elecs.size() > 1))
    {
      ElCentProb = 0.50;
    }
    while (Nct <= Npre)
    {
      acc = MCMove(Atoms,Elecs,LJparams);
      if (acc)
      {
        Nct += 1;
      }
    }
    BeadProb = tmp1;
    CentProb = tmp2;
    ElBeadProb = tmp3;
    ElCentProb = tmp4;
  }
  cout << "Starting equilibration..." << '\n';
  if (Debug == 1)
  {
    Print_traj(Atoms,Elecs,trajfile,PrintMode);
  }
  Nct = 0;
  while (Nct <= Neq) //Equilibration
  {
    if(ct == Acc_Check)
    {
      if (Debug == 1)
      {
        cout << "Step: ";
        cout << Nct;
        cout << " Energy: ";
        double Veff = Ek; //Kinetic energy
        if (Atoms.size() > 0)
        {
          //Add atom-atom interactions
          Veff += Get_Epot(Atoms,LJparams);
          Veff -= Get_Espring(Atoms);
        }
        if (Elecs.size() > 0)
        {
          //Add atom-lepton and lepton-lepton interaction
          Veff += Get_EeFF(Atoms,Elecs);
          Veff -= Get_Espring(Elecs)+RadSpring(Elecs);
        }
        if (Ensemble == 1)
        {
          Veff += Press*Lx*Ly*Lz*atm2eV;
        }
        cout << Veff;
        cout << " Accept ratio: ";
        cout << (Nacc/(Nrej+Nacc));
        cout << " Step size: ";
        cout << step << '\n';
        Print_traj(Atoms,Elecs,trajfile,PrintMode);
      }
      if ((Nacc/(Nrej+Nacc)) > accratio)
      {
        step *= 1.10;
      }
      if ((Nacc/(Nrej+Nacc)) < accratio)
      {
        step *= 0.91;
      }
      if (step < StepMin)
      {
        step = StepMin;
      }
      if (step > StepMax)
      {
        step = StepMax;
      }
      ct = 0;
      Nacc = 0;
      Nrej = 0;
    }
    ct += 1;
    acc = MCMove(Atoms,Elecs,LJparams);
    if (acc)
    {
      Nct += 1;
      Nacc += 1;
    }
    else
    {
      Nrej += 1;
    }
  }
  Nct = 0;
  Nacc = 0;
  Nrej = 0;
  ct = 0;
  cout << "Starting production run..." << '\n';
  Print_traj(Atoms,Elecs,trajfile,PrintMode);
  while (Nct < Nsteps)
  {
    acc = MCMove(Atoms,Elecs,LJparams);
    if (acc)
    {
      Nct += 1;
      ct += 1;
      Nacc += 1;
      double Et = Ek;
      if (Atoms.size() > 0)
      {
        Et += Get_Epot(Atoms,LJparams);
        Et -= Get_Espring(Atoms);
      }
      if (Elecs.size() > 0)
      {
        Et += Get_EeFF(Atoms,Elecs);
        Et -= Get_Espring(Elecs)+RadSpring(Elecs);
      }
      if (Ensemble == 1)
      {
        Et += Press*Lx*Ly*Lz*atm2eV;
      }
      VolAvg += Lx*Ly*Lz;
      SumE += Et;
      SumE2 += Et*Et;
      if (ct == Nprint)
      {
        Print_traj(Atoms,Elecs,trajfile,PrintMode);
        ct = 0;
      }
    }
    else
    {
      Nrej += 1;
    }
  }
  SumE /= Nsteps;
  SumE2 /= Nsteps;
  VolAvg /= Nsteps;

  //Print output
  cout << '\n';
  cout << "Temperature: ";
  cout << 1.0/(k*Beta);
  cout << " Volume: ";
  cout << VolAvg;
  cout << '\n';
  cout << "Average Energy: ";
  cout << SumE;
  cout << " Variance: ";
  cout << (SumE2-SumE*SumE);
  cout << '\n';
  cout << "Acceptance ratio: ";
  cout << (Nacc/(Nrej+Nacc));
  cout << " Step size: ";
  cout << step;
  cout << '\n' << '\n';
  //End of section

  //Clean up and quit
  if (Debug == 1)
  {
    cout << "Cleaning up..." << '\n';
  }
  paramfile.close();
  trajfile.close();
  potfile.close();
  xyzfile.close();
  bondfile.close();
  cout << "Done." << '\n' << '\n';
  return 0;
};
