#include "safio.h"
#include "string_utils.h"
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

        while (getline(safio_input, line))
        {
            debug_file << line << '\n';

            //Skip blank lines
            if (line == "")
                continue;

            //Splits line for parsing.
            std::vector<std::string> args = split(line);

            if (n == 1)
            {
                E0 = atof(args[0].c_str());
                THETA0 = atof(args[1].c_str());
                PHI0 = atof(args[2].c_str());
                MASS = atof(args[3].c_str());
                SYMION = &args[4][0];

                //initialize the atom for this ion.
                ion.mass = MASS;
                ion.symbol = SYMION;
                ion.index = 0;
                //a.charge = TODO lookup

                // Ensure phi is in correct range
                while (PHI0 > 180)
                    PHI0 -= 360;
                while (PHI0 < -180)
                    PHI0 += 360;
            }
            if (n == 2)
            {
                EMIN = atof(args[0].c_str());
                EMAX = atof(args[1].c_str());
                ESIZE = atof(args[2].c_str());
                ASIZE = atof(args[3].c_str());
            }
            if (n == 3)
            {
                NDTECT = atoi(args[0].c_str());
            }
            if (n == 4)
            {
                DTECTPAR = to_double_array(args, 0, 3);
            }
            if (n == 5)
            {
                DELLOW = atof(args[0].c_str());
                DELT0 = atof(args[1].c_str());
            }
            if (n == 6)
            {
                DEMAX = atof(args[0].c_str());
                DEMIN = atof(args[1].c_str());
                ABSERR = atof(args[2].c_str());
            }
            if (n == 7)
            {
                NPART = atoi(args[0].c_str());
            }
            if (n == 8)
            {
                RECOIL = args[0] == "t";
            }
            if (n == 9)
            {
                Z1 = atof(args[0].c_str());
            }
            if (n == 10)
            {
                MAX_STEPS = atoi(args[0].c_str());
            }
            if (n == 11)
            {
                R_MAX = atof(args[0].c_str());
                rr_max = R_MAX * R_MAX;
                DR_MIN_TAB = atof(args[1].c_str());
            }
            if (n == 12)
            {
                ZMIN = atof(args[0].c_str());
                ZSTEP = atof(args[1].c_str());
            }
            if (n == 13)
            {
                if (o <= 0)
                {
                    DIST_SEARCH = atoi(args[0].c_str());
                    FAILED_DE = atoi(args[1].c_str());

                    //Since we re-purposed the above, we always
                    //have the 4 successive arguments.
                    o = 4;
                }
                else
                {
                    if (o == 3)
                    {
                        NUMCHA = atoi(args[0].c_str());
                    }
                    else if (o == 2)
                    {
                        XSTART = atof(args[0].c_str());
                        XSTEP = atof(args[1].c_str());
                        XSTOP = atof(args[2].c_str());
                    }
                    else if (o == 1)
                    {
                        YSTART = atof(args[0].c_str());
                        YSTEP = atof(args[1].c_str());
                        YSTOP = atof(args[2].c_str());
                    }
                }
            }
            if (n == 14)
            {
                SCAT_FLAG = atoi(args[0].c_str());
                SCAT_TYPE = atoi(args[1].c_str());
            }
            if (n == 15)
            {
                RAX = atof(args[0].c_str());
                RAY = atof(args[1].c_str());
            }
            if (n == 16)
            {
                NPAR = atoi(args[0].c_str());
                IPOT = atoi(args[1].c_str());
            }
            if (n == 17)
            {
                POTPAR = to_double_array(args, 0, NPAR - 1);
            }
            if (n == 18)
            {
                NIMPAR = atoi(args[0].c_str());
                IIMPOT = atoi(args[1].c_str());
            }
            if (n == 19)
            {
                if (o <= 0)
                {
                    PIMPAR = to_double_array(args, 0, NIMPAR - 1);
                    //If this is the case, we have more arguments
                    if (IIMPOT == 2)
                    {
                        //This will be decremented after leaving this if
                        o = 4;
                    }
                }
                else
                {
                    if (o == 6)
                    {
                        NBZ = atoi(args[0].c_str());
                        TOL = atof(args[1].c_str());
                    }
                    else if (o == 5)
                    {
                        ZMAX = to_double_array(args, 0, NBZ - 1);
                    }
                    else if (o == 4)
                    {
                        NZ = to_double_array(args, 0, NBZ - 1);
                    }
                    else if (o == 3)
                    {
                        NBG = atoi(args[0].c_str());
                        GTOL = atof(args[1].c_str());
                    }
                    else if (o == 2)
                    {
                        GMAX = to_double_array(args, 0, NBG - 1);
                    }
                    else if (o == 1)
                    {
                        NG = to_double_array(args, 0, NBG - 1);
                    }
                }
            }
            if (n == 20)
            {
                TEMP = atof(args[0].c_str());
                SEED = atof(args[1].c_str());
                NITER = atoi(args[2].c_str());
            }
            if (n == 21)
            {
                IMAGE = args[0] == "t";
            }
            if (n == 22)
            {
                SENRGY = atof(args[0].c_str());
                BDIST = atof(args[1].c_str());
            }
            if (n == 23)
            {
                AX = atof(args[0].c_str());
                AY = atof(args[1].c_str());
                AZ = atof(args[2].c_str());
            }
            if (n == 24)
            {
                if (o <= 0)
                {
                    NBASIS = atoi(args[0].c_str());
                    o = NBASIS + 1;
                }
                else
                {
                    Site s;
                    s.r_0[0] = atof(args[0].c_str());
                    s.r_0[1] = atof(args[1].c_str());
                    s.r_0[2] = atof(args[2].c_str());
                    s.index = atoi(args[3].c_str());
                    BASIS.push_back(s);
                }
            }
            if (n == 25)
            {
                if (o <= 0)
                {
                    NTYPES = atoi(args[0].c_str());
                    o = NTYPES * 2 + 1;
                }
                else
                {
                    if (o % 2 == 0)
                    {
                        Atom a;
                        a.mass = atof(args[0].c_str());
                        a.charge = atof(args[1].c_str());
                        a.symbol = args[2];
                        a.index = o / 2;
                        ATOMS.push_back(a);
                    }
                    else
                    {
                        Atom &a = ATOMS.back();
                        a.spring[0] = atof(args[0].c_str());
                        a.spring[1] = atof(args[1].c_str());
                        a.spring[2] = atof(args[2].c_str());
                    }
                }
            }
            if (n == 26)
            {
                CORR = args[0] == "t";
                ATOMK = atof(args[1].c_str());
                RNEIGH = atof(args[2].c_str());
            }
            if (n == 27)
            {
                face = to_double_array(args, 0, 2);
                load_crystal = args.size() > 3 and args[3] == "t"; 
                if(load_crystal)
                {
                    loaded_face = to_double_array(args, 4,6);
                }
                else
                {
                    //Default basis should be defined in 001 direction
                    loaded_face = new double[3];
                    loaded_face[0] = 0;
                    loaded_face[1] = 0;
                    loaded_face[2] = 1;
                }
                
            }
            // Decrement our sub-line first.
            o--;
            // Only increment this if not in a sub-line section.
            if (o <= 0)
                n++;

        }

        //Populate basis atoms indecies
        for(int i = 0; i<NBASIS; i++)
        {
            BASIS[i].atom = &ATOMS[BASIS[i].index-1];
        }


        safio_input.close();
        debug_file << "\nLoaded SAFIO" << '\n';

        //Then we need these files.
        if (NUMCHA == 1)
        {
            traj_file.open(safio_file + ".traj");
            xyz_file.open(safio_file + ".xyz");
            debug_file << "Initializing for single shot mode." << '\n';
        }
    }
    else
    {
        std::cout << "Error opening Safio File" << '\n';
    }
    return;
}

double zeros[3] = { 0,0,0 };
void Site::reset()
{
    //Reset positions and momenta
    std::copy(r_0, r_0 + 3, r);
    std::copy(r_0, r_0 + 3, r_t);
    std::copy(p_0, p_0 + 3, p);
    //Reset forces
    std::copy(std::begin(zeros), std::end(zeros), dp_dt);
    std::copy(std::begin(zeros), std::end(zeros), dp_dt_t);
}

void Site::write_info()
{
    debug_file << "Atom: " << atom->symbol << std::endl;
    debug_file << "r  : " << r[0] << " " << r[1] << " " << r[2] << std::endl;
    debug_file << "p  : " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    debug_file << "r_t: " << r_t[0] << " " << r_t[1] << " " << r_t[2] << std::endl;
    debug_file << "F: " << dp_dt[0] << " " << dp_dt[1] << " " << dp_dt[2] << std::endl;
    debug_file << "F_t: " << dp_dt_t[0] << " " << dp_dt_t[1] << " " << dp_dt_t[2] << std::endl;
}
