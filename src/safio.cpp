#include "safio.h"
#include "string_uitls.h"
#include <vector>

void Safio::load(std::string safio_file)
{
    settings = *this;
    std::ifstream safio_input;
    std::string filename = safio_file + ".input";
    safio_input.open(filename);

    filename = safio_file + ".data";
    out_file.open(filename);

    filename = safio_file + ".crys";
    crystal_file.open(filename);

    filename = safio_file + ".dbug";
    debug_file.open(filename);
    debug_file << "Loading Info From: " << safio_file << "\n\n";

    if (safio_input.is_open())
    {
        std::string line;
        //This is line index,
        //matching first line of file
        int n = 1;
        //This is used for offsets for
        //values which affect line number
        int o = 0;

        while(getline(safio_input,line))
        {
            debug_file << line << '\n';

            //Skip blank lines
            if(line=="")
                continue;

            //Splits line for parsing.
            std::vector<std::string> args = split(line);

            if(n==1)
            {
                E0 = atof(args[0].c_str());
                THETA0 = atof(args[1].c_str());
                PHI0 = atof(args[2].c_str());
                MASS = atof(args[3].c_str());
                SYMION = &args[4][0];
                // Ensure phi is in correct range
                while (PHI0 > 180)
                    PHI0 -= 360;
                while (PHI0 < -180)
                    PHI0 += 360;
            }
            if(n==2)
            {
                EMIN = atof(args[0].c_str());
                EMAX = atof(args[1].c_str());
                ESIZE = atof(args[2].c_str());
                ASIZE = atof(args[3].c_str());
            }
            if(n==3)
            {
                NDTECT = atoi(args[0].c_str());
            }
            if(n==4)
            {
                DTECTPAR = to_double_array(args, 0, 3);
            }
            if(n==5)
            {
                DELLOW = atof(args[0].c_str());
                DELT0 = atof(args[1].c_str());
            }
            if(n==6)
            {
                DEMAX = atof(args[0].c_str());
                DEMIN = atof(args[1].c_str());
                ABSERR = atof(args[2].c_str());
            }
            if(n==7)
            {
                NPART = atoi(args[0].c_str());
            }
            if(n==8)
            {
                RECOIL = args[0]=="t";
            }
            if(n==9)
            {
                Z1 = atof(args[0].c_str());
            }
            if(n==10)
            {
                MAX_STEPS = atoi(args[0].c_str());
            }
            if(n==11)
            {
                R_MAX = atof(args[0].c_str());
                rr_max = R_MAX * R_MAX;
                DR_MIN_TAB = atof(args[1].c_str());
            }
            if(n==12)
            {
                ZMIN = atof(args[0].c_str());
                ZSTEP = atof(args[1].c_str());
            }
            if(n==13)
            {
                if(o<=0)
                {
                    MAXDIV = atoi(args[0].c_str());
                    MINDIV = atoi(args[1].c_str());
                    //If this is the case, we have more arguments
                    if(MAXDIV == MINDIV && MAXDIV == 1)
                    {
                        o = 4;
                    }
                }
                else
                {
                    if(o==3)
                    {
                        NUMCHA = atoi(args[0].c_str());
                    }
                    else if(o==2)
                    {
                        XSTART = atof(args[0].c_str());
                        XSTEP = atof(args[1].c_str());
                        XSTOP = atof(args[2].c_str());
                    }
                    else if(o==1)
                    {
                        YSTART = atof(args[0].c_str());
                        YSTEP = atof(args[1].c_str());
                        YSTOP = atof(args[2].c_str());
                    }
                }
            }
            if(n==14)
            {
                NWRITX = atoi(args[0].c_str());
                NWRITY = atoi(args[1].c_str());
            }
            if(n==15)
            {
                RAX = atof(args[0].c_str());
                RAY = atof(args[1].c_str());
            }
            if(n==16)
            {
                NPAR = atoi(args[0].c_str());
                IPOT = atoi(args[1].c_str());
            }
            if(n==17)
            {
                POTPAR = to_double_array(args, 0, NPAR - 1);
            }
            if(n==18)
            {
                NIMPAR = atoi(args[0].c_str());
                IIMPOT = atoi(args[1].c_str());
            }
            if(n==19)
            {
                if(o<=0)
                {
                    PIMPAR = to_double_array(args, 0, NIMPAR - 1);
                    //If this is the case, we have more arguments
                    if(IIMPOT == 2)
                    {
                        //This will be decremented after leaving this if
                        o = 4;
                    }
                }
                else
                {
                    if(o==6)
                    {
                        NBZ = atoi(args[0].c_str());
                        TOL = atof(args[1].c_str());
                    }
                    else if(o==5)
                    {
                        ZMAX = to_double_array(args, 0, NBZ - 1);
                    }
                    else if(o==4)
                    {
                        NZ = to_double_array(args, 0, NBZ - 1);
                    }
                    else if(o==3)
                    {
                        NBG = atoi(args[0].c_str());
                        GTOL = atof(args[1].c_str());
                    }
                    else if(o==2)
                    {
                        GMAX = to_double_array(args, 0, NBG - 1);
                    }
                    else if(o==1)
                    {
                        NG = to_double_array(args, 0, NBG - 1);
                    }
                }
            }
            if(n==20)
            {
                TEMP = atof(args[0].c_str());
                SEED = atof(args[1].c_str());
                NITER = atoi(args[2].c_str());
            }
            if(n==21)
            {
                IMAGE = args[0]=="t";
            }
            if(n==22)
            {
                SENRGY = atof(args[0].c_str());
                BDIST = atof(args[1].c_str());
            }
            if(n==23)
            {
                AX = atof(args[0].c_str());
                AY = atof(args[1].c_str());
                AZ = atof(args[2].c_str());
            }
            if(n==24)
            {
                if(o<=0)
                {
                    NBASIS = atoi(args[0].c_str());
                    o = NBASIS + 1;
                }
                else
                {
                    Site s;
                    s[0] = atof(args[0].c_str());
                    s[1] = atof(args[1].c_str());
                    s[2] = atof(args[2].c_str());
                    s.index =  atoi(args[3].c_str());
                    BASIS.push_back(s);
                }
            }
            if(n==25)
            {
                if(o<=0)
                {
                    NTYPES = atoi(args[0].c_str());
                    o = NTYPES * 2  + 1;
                }
                else
                {
                    if(o%2==0)
                    {
                        Atom a;
                        a.mass = atof(args[0].c_str());
                        a.charge = atof(args[1].c_str());
                        a.symbol = args[2];
                        a.index = o/2;
                        ATOMS.push_back(a);
                    }
                    else
                    {
                        Atom a = ATOMS.back();
                        a.spring[0] = atof(args[0].c_str());
                        a.spring[1] = atof(args[1].c_str());
                        a.spring[2] = atof(args[2].c_str());
                    }
                }
            }
            if(n==26)
            {
                CORR = args[0]=="t";
                ATOMK = atof(args[1].c_str());
                RNEIGH = atof(args[2].c_str());
            }
            if(n==27)
            {
                face = to_double_array(args, 0, 2);
            }
            // Decrement our sub-line first.
            o--;
            // Only increment this if not in a sub-line section.
            if (o <= 0)
                n++;

        }
        safio_input.close();
        debug_file << "\nLoaded SAFIO" << '\n';

        //Then we need these files.
        if(NUMCHA == 1)
        {
            traj_file.open(safio_file+".traj");
            xyz_file.open(safio_file+".xyz");
            debug_file << "Initializing for single shot mode." << '\n';
        }
    }
    else
    {
        std::cout << "Error opening Safio File" << '\n';
    }
    return;
}
