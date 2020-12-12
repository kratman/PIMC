
#pragma once

inline bool MCMove(std::vector<Qpart>& parts, std::vector<Qpart>& elecs,
                   std::map<std::string,LennardJones>& LJmap)
{
    bool acc = 0;
    //Copy parts
    std::vector<Qpart> parts2;
    std::vector<Qpart> elecs2;
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
        std::cout << "Warning two structures have identical energies." << '\n';
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
