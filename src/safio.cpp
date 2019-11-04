#include "safio.h"
#include "string_utils.h"
#include "temps.h"
#include <vector>

void Safio::load(std::map<std::string, ArgValue>& args)
{
    //This should have been included to args if not originally present.
    std::string safio_file = args["-i"].as_string();

    settings = *this;
    std::ifstream safio_input;
    std::string filename = safio_file + ".input";
    safio_input.open(filename);

    //check if we have alternate output file name
    //All other streams after this will use the output file name instead.
    if(args["-o"])
    {
        safio_file = args["-o"].as_string();
        std::cout << "Output files: " << safio_file << std::endl;
    }

    //Only open these if flagged, 
    //this allows modules to use safio, but not open the files
    if(args["-f"])
    {
        filename = safio_file + ".data";
        out_file.open(filename);

        filename = safio_file + ".crys";
        crystal_file.open(filename);
    }

    //Open the debug file, this is alway used.
    filename = safio_file + ".dbug";
    debug_file.open(filename);
    //This one needs to specify the old file name.
    debug_file << "Loading Info From: " << args["-i"].as_string() << "\n\n";

    //Number of parameters for potentials.
    //This is used to know how many values to look for in a line
    int npar = -1;

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
            //Skip blank lines
            if (line == "")
                continue;
            //Allow having comment lines in the file
            //Comment lines start with a #
            if(starts_with(line, "#")) continue;

            //Print the input line to debug file.
            debug_file << line << '\n';

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
                detector_type = atoi(args[0].c_str());
            }
            if (n == 4)
            {
                detect_parameters = to_double_array(args, 0, 3);
            }
            if (n == 5)
            {
                DELLOW = atof(args[0].c_str());
                DELT0 = atof(args[1].c_str());
            }
            if (n == 6)
            {
                error_exponent = atof(args[0].c_str());
                DEMIN = atof(args[1].c_str());
                error_scale = atof(args[2].c_str());
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
                npar = atoi(args[0].c_str());
                binary_potential_type = atoi(args[1].c_str());
            }
            if (n == 17)
            {
                binary_potential_parameters = to_double_array(args, 0, npar - 1);
            }
            if (n == 18)
            {
                npar = atoi(args[0].c_str());
                image_potential_type = atoi(args[1].c_str());
            }
            if (n == 19)
            {
                if (o <= 0)
                {
                    image_parameters = to_double_array(args, 0, npar - 1);
                    //If this is the case, we have more arguments
                    if (image_potential_type == 2)
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
                ion_index = atoi(args[2].c_str());
            }
            if (n == 21)
            {
                use_image = args[0] == "t";
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
            if(n == 28)
            {
                F_a = atof(args[0].c_str());
                F_b = atof(args[1].c_str());
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

        //check arguments for specifically enabling single-shot mode
        if(args["-s"].as_bool())
        {
            NUMCHA = 1;
            SCAT_FLAG = 666;
            //Set x-start
            if(args["-x"])
            {
                XSTART = args["-x"].as_double();
            }
            //Set y-start
            if(args["-y"])
            {
                YSTART = args["-y"].as_double();
            }
            //If this is true, it will output only nearish to xyz
            SCAT_TYPE = args["-r"].as_bool();

            debug_file << "After commandline args: " << '\n';
            debug_file << "NUMCHA: " << NUMCHA << '\n';
            debug_file << "SCAT_FLAG: " << SCAT_FLAG << '\n';
            debug_file << "XSTART: " << XSTART << '\n';
            debug_file << "YSTART: " << YSTART << '\n';
            debug_file << "SCAT_TYPE: " << SCAT_TYPE << '\n';
        }

        safio_input.close();
        debug_file << "\nLoaded SAFIO" << '\n';
        std::cout << "\nLoaded SAFIO" << '\n';

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
