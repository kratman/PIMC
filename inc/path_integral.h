
#pragma once

#include "physical_constants.h"

inline double SpringEnergy(double k, double r2)
{
    //General harmonic bond for PI rings
    double E = 0.5*k*r2;
    return E;
};

inline double Get_Espring(vector<Qpart>& parts)
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

inline double Get_Epot(vector<Qpart>& parts, map<string,LennardJones>& LJmap)
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
}
