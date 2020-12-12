
#pragma once

#include <omp.h>

ParticleCoordinates Get_COM(std::vector<Qpart>& parts, std::vector<Qpart>& elecs)
{
    double x=0,y=0,z=0,M=0;
    ParticleCoordinates com;
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

void Get_Range(std::vector<Qpart>& parts, ParticleCoordinates& Ls, ParticleCoordinates& Hs)
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

void Remove_COM(std::vector<Qpart>& parts, std::vector<Qpart>& elecs)
{
    ParticleCoordinates com = Get_COM(parts, elecs);
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
