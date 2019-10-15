#include "safio.h"
#include <vector>

void safio::load()
{
    settings = *this;
    std::string safio_file = "safari.input";
    char* file = &safio_file[0];
    std::ifstream input;
    input.open(file);
    if (input.is_open())
    {
        while(getline(input,safio_file))
        {
            std::cout << safio_file << '\n';
            break;
        }
        file = &safio_file[0];
        input.close();
        std::ifstream safio_input;
        std::cout << file << std::endl;
        safio_input.open(file);

        std::string output = safio_file;
        output.replace(output.end()-5, output.end(), "data");
        out_file.open(output);
        output.replace(output.end()-4, output.end(), "dbug");
        debug_file.open(output);

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
                int i = 0;

                //Skip blank lines
                if(line=="")
                    continue;

                std::istringstream iss(line);

                //TODO decide on max number or args.
                std::string args[10];

                //Split the string into an array.
                do
                {
                    std::string subs;
                    iss >> subs;
                    //cout << "Substring: " << subs << endl;
                    args[i] = subs;
                    i++;
                }
                while (iss);

                if(n==1)
                {
                    E0 = atof(args[0].c_str());
                    THETA0 = atof(args[1].c_str());
                    PHI0 = atof(args[2].c_str());
                    MASS = atof(args[3].c_str());
                    SYMION = args[4];
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
                    DTECTPAR.push_back(atof(args[0].c_str()));
                    DTECTPAR.push_back(atof(args[1].c_str()));
                    DTECTPAR.push_back(atof(args[2].c_str()));
                    DTECTPAR.push_back(atof(args[3].c_str()));
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
                    NTAB = atoi(args[0].c_str());
                }
                if(n==11)
                {
                    RRMIN = atof(args[0].c_str());
                    RRSTEP = atof(args[1].c_str());
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
                    FAX = atof(args[0].c_str());
                    FAY = atof(args[1].c_str());
                }
                if(n==16)
                {
                    NPAR = atoi(args[0].c_str());
                    IPOT = atoi(args[1].c_str());
                }
                if(n==17)
                {
                    for (int i = 0; i < NPAR; i++)
                    {
                        POTPAR.push_back(atof(args[i].c_str()));
                    }
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
                        for (int i = 0; i < NIMPAR; i++)
                        {
                            PIMPAR.push_back(atof(args[i].c_str()));
                        }
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
                            for (int i = 0; i < NBZ; i++)
                            {
                                ZMAX.push_back(atof(args[i].c_str()));
                            }
                        }
                        else if(o==4)
                        {
                            for (int i = 0; i < NBZ; i++)
                            {
                                NZ.push_back(atof(args[i].c_str()));
                            }
                        }
                        else if(o==3)
                        {
                            NBG = atoi(args[0].c_str());
                            GTOL = atof(args[1].c_str());
                        }
                        else if(o==2)
                        {
                            for (int i = 0; i < NBG; i++)
                            {
                                GMAX.push_back(atof(args[i].c_str()));
                            }
                        }
                        else if(o==1)
                        {
                            for (int i = 0; i < NBG; i++)
                            {
                                NG.push_back(atof(args[i].c_str()));
                            }
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
                    if(i>2)
                        AZ = atof(args[2].c_str());
                    else
                        AZ = AX;
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
                        site s;
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
                            atom a;
                            a.mass = atof(args[0].c_str());
                            a.charge = atof(args[1].c_str());
                            a.symbol = args[2];
                            a.index = o/2;
                            ATOMS.push_back(a);
                        }
                        else
                        {
                            atom a = ATOMS.back();
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
                    face[0] = atof(args[0].c_str());
                    face[1] = atof(args[1].c_str());
                    face[2] = atof(args[2].c_str());
                }
                // Decrement our sub-line first.
                o--;
                // Only increment this if not in a sub-line section.
                if (o <= 0)
                    n++;

            }
            safio_input.close();
        }
        else
        {
            std::cout << "Error opening Safio File" << '\n';
        }
    }
    else
    {
        std::cout << "Error opening Safari.input" << '\n';
    }
    return;
}
