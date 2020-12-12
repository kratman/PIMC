
#pragma once

#include "center_of_mass.h"

inline void Print_traj(std::vector<Qpart>& parts, std::vector<Qpart>& elecs,
                std::fstream& traj, std::string mode)
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
