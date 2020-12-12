
#include <ctime>
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <map>

#include "coordinates.h"
#include "harmonic_bond.h"
#include "harmonic_angle.h"
#include "lennard_jones.h"
#include "electron.h"
#include "move_settings.h"
#include "particle.h"
#include "physics_settings.h"
#include "path_integral.h"
#include "system.h"
#include "monte_carlo.h"
#include "output.h"

int main()
{
  std::cout << '\n';
  //Initialize dynamic parameters
  srand((unsigned)time(0)); //Serial only random numbers
  std::string dummy,PrintMode,SpinMode;
  std::string filename1 = "pimc.param"; //Input file
  std::string filename2 = "traj.xyz"; //Output trajectory
  std::string filename3 = "pimc.pot"; //Potential file
  std::string filename4; //xyz input
  std::string filename5; //Connectivity input
  std::fstream paramfile,trajfile,potfile,xyzfile,bondfile;
  int Natoms,Nsteps,Ntyps,Neq,Nprint,Npre,Nelec,Npos;
  double accratio,SumE,SumE2,VolAvg,Ek;
  std::vector<Qpart> Atoms;
  std::vector<std::string> Types;
  std::map<std::string,LennardJones> LJparams;
  paramfile.open(filename1.c_str(),std::ios_base::in);
  trajfile.open(filename2.c_str(),std::ios_base::out);
  potfile.open(filename3.c_str(),std::ios_base::in);
  //End of section

  //Read Input
  if (Debug)
  {
    std::cout << "Reading input..." << '\n';
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
  xyzfile.open(filename4.c_str(),std::ios_base::in);
  xyzfile >> Natoms;
  for (int i=0;i<Natoms;i++)
  {
    std::string typ;
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
      if (qGrid)
      {
        randx = 0.5;
        randy = 0.5;
        randz = 0.5;
      }
      ParticleCoordinates temp;
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
    std::string tmp;
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
    LennardJones tmp;
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
        LennardJones tmp;
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
      bondfile.open(filename5.c_str(),std::ios_base::in);
      int Nbonds;
      bondfile >> dummy >> Nbonds;
      for (int i=0;i<Nbonds;i++)
      {
        BondParameters tmp,tmp2;
        bondfile >> tmp.at2 >> tmp2.at2 >> tmp.K >> tmp.R;
        tmp2.K = tmp.K;
        tmp2.R = tmp.R;
        if ((tmp.at2 > Atoms.size()) or (tmp.at2 > Atoms.size()))
        {
          std::cout << "Bond to non-existent atom!!!!" << '\n';
          std::cout << "Check connectivity file." << '\n';
          return 0;
        }
        Atoms[tmp2.at2].Bonds.push_back(tmp);
        Atoms[tmp.at2].Bonds.push_back(tmp2);
      }
      //Collect angle info
      bondfile >> dummy >> Nbonds;
      for (int i=0;i<Nbonds;i++)
      {
        AngleParameters tmp;
        int cent;
        bondfile >> cent >> tmp.at2 >> tmp.at3 >> tmp.K >> tmp.ang;
        tmp.ang *= pi/180.0;
        if ((tmp.at2 > Atoms.size()) or (tmp.at2 > Atoms.size())
        or (cent > Atoms.size()))
        {
          std::cout << "Angle with non-existent atom!!!!" << '\n';
          std::cout << "Check connectivity file." << '\n';
          return 0;
        }
        Atoms[cent].Angs.push_back(tmp);
      }
    }
  }
  std::vector<Qpart> Elecs;
  if (Nelec > 0)
  {
    int spinct = -1;
    ParticleCoordinates Lows;
    ParticleCoordinates Highs;
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
        ParticleCoordinates Temp;
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
    ParticleCoordinates Lows;
    ParticleCoordinates Highs;
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
        ParticleCoordinates Temp;
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
  std::cout << "Setting up simulation..." << '\n';
  std::cout << '\n';
  std::cout << "Atoms: " << Natoms << '\n';
  std::cout << "Electrons: " << Nelec << '\n';
  std::cout << "Positrons: " << Npos << '\n';
  std::cout << "Beads: " << Nbeads << '\n';
  std::cout << '\n';
  if ((Npre > 0) and (Elecs.size() > 0))
  {
    std::cout << "Pre-equilibration steps for electrons: ";
    std::cout << Npre << '\n';
  }
  std::cout << "Equilibration steps: " << Neq << '\n';
  std::cout << "Steps for production run: " << Nsteps << '\n';

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
    std::cout << "Mode: ";
    if (Ensemble == 0)
    {
      std::cout << "NVT";
    }
    if (Ensemble == 1)
    {
      std::cout << "NPT";
    }
    if (Elecs.size() > 0)
    {
      std::cout << " Low-Spin";
    }
    std::cout << " PIMC" << '\n';
  }
  else
  {
    std::cout << "Mode: ";
    if (Ensemble == 0)
    {
      std::cout << "NVT";
    }
    if (Ensemble == 1)
    {
      std::cout << "NPT";
    }
    if (Elecs.size() > 0)
    {
      std::cout << " High-Spin";
    }
    std::cout << " PIMC" << '\n';
  }
  if (Ensemble == 0)
  {
    VolProb = 0.0;
  }
  if (qGrid)
  {
    CentProb = 0.0;
    BeadProb = 0.0;
  }

  //Run simulations
  std::cout << '\n';
  Beta = 1/(k*Beta); //Invert temperature
  rtNbeads = pow(Nbeads,Scale_POW);
  SumE = 0;
  SumE2 = 0;
  VolAvg = 0;
  Ek = 0.0;
  if ((Atoms.size() > 0) and (!qGrid))
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
  if (RCOM)
  {
    Remove_COM(Atoms,Elecs);
  }
  if ((Npre > 0) and (Elecs.size() > 0))
  {
    std::cout << "Starting electron pre-equilibration..." << '\n';
    double tmp1 = BeadProb;
    if (Atoms.size() > 0)
    {
      BeadProb = 0.10;
    }
    if (qGrid)
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
  std::cout << "Starting equilibration..." << '\n';
  if (Debug)
  {
    Print_traj(Atoms,Elecs,trajfile,PrintMode);
  }
  Nct = 0;
  while (Nct <= Neq) //Equilibration
  {
    if(ct == Acc_Check)
    {
      if (Debug)
      {
        std::cout << "Step: ";
        std::cout << Nct;
        std::cout << " Energy: ";
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
        std::cout << Veff;
        std::cout << " Accept ratio: ";
        std::cout << (Nacc/(Nrej+Nacc));
        std::cout << " Step size: ";
        std::cout << step << '\n';
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
  std::cout << "Starting production run..." << '\n';
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
  std::cout << '\n';
  std::cout << "Temperature: ";
  std::cout << 1.0/(k*Beta);
  std::cout << " Volume: ";
  std::cout << VolAvg;
  std::cout << '\n';
  std::cout << "Average Energy: ";
  std::cout << SumE;
  std::cout << " Variance: ";
  std::cout << (SumE2-SumE*SumE);
  std::cout << '\n';
  std::cout << "Acceptance ratio: ";
  std::cout << (Nacc/(Nrej+Nacc));
  std::cout << " Step size: ";
  std::cout << step;
  std::cout << '\n' << '\n';
  //End of section

  //Clean up and quit
  if (Debug)
  {
    std::cout << "Cleaning up..." << '\n';
  }
  paramfile.close();
  trajfile.close();
  potfile.close();
  xyzfile.close();
  bondfile.close();
  std::cout << "Done." << '\n' << '\n';
  return 0;
};
