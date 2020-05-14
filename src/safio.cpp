#include "safio.h"

#include "string_utils.h" //ArgValue, to_double_array, etc

void findAndReplaceAll(std::string &data, std::string toSearch, std::string replaceStr)
{
    // Get the first occurrence
    size_t pos = data.find(toSearch);

    // Repeat till end is reached
    while (pos != std::string::npos)
    {
        // Replace this occurrence of Sub String
        data.replace(pos, toSearch.size(), replaceStr);
        // Get the next occurrence from the current position
        pos = data.find(toSearch, pos + replaceStr.size());
    }
}

void Safio::load(std::map<std::string, ArgValue> &prog_args)
{
    //This should have been included to prog_args if not originally present.
    input_name = prog_args["-i"].as_string();
    output_name = input_name;

    settings = *this;
    std::ifstream safio_input;
    std::string filename = output_name + ".input";
    safio_input.open(filename);

    //check if we have alternate output file name
    //All other streams after this will use the output file name instead.
    if (prog_args["-o"])
    {
        output_name = prog_args["-o"].as_string();
    }
    std::cout << "Output files: " << output_name << std::endl;

    //Only open these if flagged,
    //this allows modules to use safio, but not open the files
    if (prog_args["-f"])
    {
        filename = output_name + ".data";
        out_file.open(filename);

        filename = output_name + ".crys";
        crystal_file.open(filename);
    }

    //Open the debug file, this is alway used.
    filename = output_name + ".dbug";
    debug_file.open(filename);
    //This one needs to specify the old file name.
    debug_file << "Loading Info From: " << prog_args["-i"].as_string() << "\n\n";

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
            findAndReplaceAll(line, "\r", "");
            findAndReplaceAll(line, "\n", "");
            //Skip blank lines
            if (line == "")
                continue;
            if (starts_with(line, " "))
                continue;
            if (starts_with(line, "\t"))
                continue;
            //Allow having comment lines in the file
            //Comment lines start with a #
            if (starts_with(line, "#"))
                continue;

            //Print the input line to debug file.
            debug_file << line << '\n';

            //Splits line for parsing.
            std::vector<std::string> line_args = split(line);

            if (n == 1)
            {
                E0 = atof(line_args[0].c_str());
                THETA0 = atof(line_args[1].c_str());
                PHI0 = atof(line_args[2].c_str());
                MASS = atof(line_args[3].c_str());
                SYMION = &line_args[4][0];

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
                EMIN = atof(line_args[0].c_str());
                EMAX = atof(line_args[1].c_str());
                ESIZE = atof(line_args[2].c_str());
                ASIZE = atof(line_args[3].c_str());
            }
            if (n == 3)
            {
                detector_type = atoi(line_args[0].c_str());
                if (line_args.size() > 1 && line_args[1] == "f")
                    save_errored = false;
            }
            if (n == 4)
            {
                detect_parameters = to_double_array(line_args, 0, line_args.size() - 1);
            }
            if (n == 5)
            {
                DELLOW = atof(line_args[0].c_str());
                DELT0 = atof(line_args[1].c_str());
            }
            if (n == 6)
            {
                error_exponent = atof(line_args[0].c_str());
                DEMIN = atof(line_args[1].c_str());
                error_scale = atof(line_args[2].c_str());
            }
            if (n == 7)
            {
                NPART = atoi(line_args[0].c_str());
            }
            if (n == 8)
            {
                RECOIL = line_args[0] == "t";
            }
            if (n == 9)
            {
                Z1 = atof(line_args[0].c_str());
            }
            if (n == 10)
            {
                MAX_STEPS = atoi(line_args[0].c_str());
            }
            if (n == 11)
            {
                R_MAX = atof(line_args[0].c_str());
                rr_max = R_MAX * R_MAX;
                DR_MIN_TAB = atof(line_args[1].c_str());
            }
            if (n == 12)
            {
                ZMIN = atof(line_args[0].c_str());
                ZSTEP = atof(line_args[1].c_str());
            }
            if (n == 13)
            {
                if (o <= 0)
                {
                    DIST_SEARCH = atoi(line_args[0].c_str());
                    FAILED_DE = atoi(line_args[1].c_str());

                    //Since we re-purposed the above, we always
                    //have the 4 successive arguments.
                    o = 4;
                }
                else
                {
                    if (o == 3)
                    {
                        NUMCHA = atoi(line_args[0].c_str());
                        singleshot = NUMCHA == 1;
                    }
                    else if (o == 2)
                    {
                        XSTART = atof(line_args[0].c_str());
                        XSTEP = atof(line_args[1].c_str());
                        XSTOP = atof(line_args[2].c_str());

                        if (line_args.size() > 3)
                        {
                            x_mask_points = to_double_array(line_args, 3, line_args.size() - 1);
                            n_x_mask = line_args.size() - 4;
                        }
                    }
                    else if (o == 1)
                    {
                        YSTART = atof(line_args[0].c_str());
                        YSTEP = atof(line_args[1].c_str());
                        YSTOP = atof(line_args[2].c_str());

                        if (line_args.size() > 3)
                        {
                            y_mask_points = to_double_array(line_args, 3, line_args.size() - 1);
                            n_y_mask = line_args.size() - 4;
                        }
                    }
                }
            }
            if (n == 14)
            {
                SCAT_FLAG = atoi(line_args[0].c_str());
                SCAT_TYPE = atoi(line_args[1].c_str());

                montecarlo = SCAT_TYPE == 666;
                gridscat = SCAT_TYPE == 777;
                chainscat = SCAT_TYPE == 888;
                adaptivegrid = SCAT_TYPE < 100;

                // 100 bifurcations is currently computationally infeasable, so
                // this will never be a valid choice for this, at least until
                // computers get many, many orders of magnitude better, when
                // that happens, here is where this needs to be changed!
            }
            if (n == 15)
            {
                RAX = atof(line_args[0].c_str());
                RAY = atof(line_args[1].c_str());
            }
            if (n == 16)
            {
                npar = atoi(line_args[0].c_str());
                binary_potential_type = atoi(line_args[1].c_str());
                lattice_potential_start = npar;
                lattice_potential_type = line_args.size() > 2 ? atoi(line_args[2].c_str()) : 0;
                if (prog_args["--latflag"])
                    lattice_potential_type = prog_args["--latflag"].as_int();
                useAtomSpings = lattice_potential_type & 1;
                useLennardJones = lattice_potential_type & 2;
                useEinsteinSprings = !(useAtomSpings || useLennardJones);
                rigidBounds = lattice_potential_type & 4;
                dynamicNeighbours = lattice_potential_type & 8;
                saveSputter = (montecarlo or gridscat) and (lattice_potential_type & 16);
            }
            if (n == 17)
            {
                binary_potential_parameters = to_double_array(line_args, 0, line_args.size() - 1);
            }
            if (n == 18)
            {
                npar = atoi(line_args[0].c_str());
                image_potential_type = atoi(line_args[1].c_str());
            }
            if (n == 19)
            {
                if (o <= 0)
                {
                    image_parameters = to_double_array(line_args, 0, npar - 1);
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
                        NBZ = atoi(line_args[0].c_str());
                        TOL = atof(line_args[1].c_str());
                    }
                    else if (o == 5)
                    {
                        ZMAX = to_double_array(line_args, 0, NBZ - 1);
                    }
                    else if (o == 4)
                    {
                        NZ = to_double_array(line_args, 0, NBZ - 1);
                    }
                    else if (o == 3)
                    {
                        NBG = atoi(line_args[0].c_str());
                        GTOL = atof(line_args[1].c_str());
                    }
                    else if (o == 2)
                    {
                        GMAX = to_double_array(line_args, 0, NBG - 1);
                    }
                    else if (o == 1)
                    {
                        NG = to_double_array(line_args, 0, NBG - 1);
                    }
                }
            }
            if (n == 20)
            {
                TEMP = atof(line_args[0].c_str());
                SEED = atof(line_args[1].c_str());
                ion_index = atoi(line_args[2].c_str());
            }
            if (n == 21)
            {
                use_image = line_args[0] == "t";
            }
            if (n == 22)
            {
                SENRGY = atof(line_args[0].c_str());
                BDIST = atof(line_args[1].c_str());
            }
            if (n == 23)
            {
                AX = atof(line_args[0].c_str());
                AY = atof(line_args[1].c_str());
                AZ = atof(line_args[2].c_str());
            }
            if (n == 24)
            {
                if (o <= 0)
                {
                    NBASIS = atoi(line_args[0].c_str());
                    o = NBASIS + 1;
                }
                else
                {
                    Site s;
                    s.r_0[0] = atof(line_args[0].c_str());
                    s.r_0[1] = atof(line_args[1].c_str());
                    s.r_0[2] = atof(line_args[2].c_str());
                    s.index = atoi(line_args[3].c_str());
                    BASIS.push_back(s);
                }
            }
            if (n == 25)
            {
                if (o <= 0)
                {
                    NTYPES = atoi(line_args[0].c_str());
                    o = NTYPES * 2 + 1;
                }
                else
                {
                    if (o % 2 == 0)
                    {
                        Atom a;
                        a.mass = atof(line_args[0].c_str());
                        a.charge = atof(line_args[1].c_str());
                        a.symbol = line_args[2];
                        a.index = o / 2;
                        ATOMS.push_back(a);
                    }
                    else
                    {
                        Atom &a = ATOMS.back();
                        a.spring[0] = atof(line_args[0].c_str());
                        a.spring[1] = atof(line_args[1].c_str());
                        a.spring[2] = atof(line_args[2].c_str());
                    }
                }
            }
            if (n == 26)
            {
                CORR = line_args[0] == "t";
                ATOMK = atof(line_args[1].c_str());
                max_spring_V = atof(line_args[2].c_str());
                neighbour_count = line_args.size() > 3 ? atoi(line_args[3].c_str()) : 1;
            }
            if (n == 27)
            {
                face = to_double_array(line_args, 0, 2);
                load_crystal = line_args.size() > 3 and line_args[3] == "t";
                if (load_crystal)
                {
                    loaded_face = to_double_array(line_args, 4, 6);
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
            if (n == 28)
            {
                F_a = atof(line_args[0].c_str());
                F_b = atof(line_args[1].c_str());
            }
            // Decrement our sub-line first.
            o--;
            // Only increment this if not in a sub-line section.
            if (o <= 0)
                n++;
        }

        //Populate basis atoms indecies
        for (int i = 0; i < NBASIS; i++)
        {
            BASIS[i].atom = &ATOMS[BASIS[i].index - 1];
        }

        if (prog_args["-n"])
        {
            NUMCHA = prog_args["-n"].as_int();
        }

        if (prog_args["-e"])
        {
            E0 = prog_args["-e"].as_double();
            debug_file << "Override of E0: " << E0 << '\n';
        }

        if (prog_args["-t"])
        {
            TEMP = prog_args["-t"].as_double();
            debug_file << "Override of Temperature: " << TEMP << '\n';
        }

        //check arguments for specifically enabling single-shot mode
        if (prog_args["-s"].as_bool())
        {
            NUMCHA = 1;
            singleshot = true;
            // We override the saving of sputter here
            saveSputter = (lattice_potential_type & 16);
            if (SCAT_TYPE)
                SCAT_TYPE = 666;
            //Set x-start
            if (prog_args["-x"])
            {
                XSTART = prog_args["-x"].as_double();
            }
            //Set y-start
            if (prog_args["-y"])
            {
                YSTART = prog_args["-y"].as_double();
            }
            //If this is true, it will output only nearish to xyz
            SCAT_TYPE = prog_args["-r"].as_bool();

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

        if (saveSputter)
        {
            sptr_file.open(output_name + ".sptr");
        }

        //Then we need these files.
        if (singleshot)
        {
            SCAT_FLAG = 666;
            if (SCAT_TYPE)
                SCAT_TYPE = 666;
            traj_file.open(output_name + ".traj");
            xyz_file.open(output_name + ".xyz");
            debug_file << "Initializing for single shot mode." << '\n';
        }
    }
    else
    {
        std::cout << "Error opening Safio File" << '\n';
    }
    return;
}
