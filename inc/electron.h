
#pragma once

#include <vector>
#include "physical_constants.h"

const double Rad = 0.10; //Initial electron radius
const double rho = 1.0; //Parameter for eFF VB mixing
const double sbar = 1.0; //Parameter for eFF radius scaling
const double rbar = 1.0; //Parameter for eFF distance scaling
const double Charm = 0; //Force constant for eFF harmonic constraint
const bool Scale_eFF = false; //Scale the eFF kinetic energy by rtNbeads
const double Scale_POW = 0.05; //rtNbeads = pow(P,Scale_POW)

inline double EFFEnergy(Qpart& atom, Qpart& elec, int p)
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

inline double EFFCorr(Qpart& elec1, Qpart& elec2, int p)
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

inline double RadSpring(std::vector<Qpart>& parts)
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
