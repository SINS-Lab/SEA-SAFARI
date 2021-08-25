#include "potentials.h"

#include "safio.h"  // settings
#include "safari.h" // exit_fail

#include <math.h> // exp, sqrt, etc

/**
 * This contains sets of 
 * 24εσ^6 and 48εσ^12
 */
double **L_J_params;
double **L_J_dV_dr_cache;
double **Vr_r_cache = NULL;
double **dVr_dr_cache = NULL;
double r_max;
double dr_min;
double z_max;
double z_min;
int n_rmax = -1;
int num_atoms = 0;

std::vector<double> **Vr_r_all_cache = NULL;
std::vector<double> **dVr_dr_all_cache = NULL;

double dVr_dr_init(double r, int n)
{
    if (settings.binary_potential_type == 1)
    {
        double a, b, c, d;

        //The potpars are in groups of 4,
        //in order of the listed basis atoms
        //The first atom is listed as 1.
        int index = (n - 1) * 4;
        a = settings.binary_potential_parameters[index];
        b = settings.binary_potential_parameters[index + 1];
        c = settings.binary_potential_parameters[index + 2];
        d = settings.binary_potential_parameters[index + 3];
        return -b * a * exp(-b * r) - d * c * exp(-d * r);
    }
    else
    {
        exit_fail("ERROR WITH dVr_dr");
    }
    return 0;
}

double Vr_r_init(double r, int n)
{
    if (settings.binary_potential_type == 1)
    {
        double a, b, c, d;
        //The potpars are in groups of 4,
        //in order of the listed basis atoms
        //The first atom is listed as 1.
        int index = (n - 1) * 4;
        a = settings.binary_potential_parameters[index];
        b = settings.binary_potential_parameters[index + 1];
        c = settings.binary_potential_parameters[index + 2];
        d = settings.binary_potential_parameters[index + 3];
        return a * exp(-b * r) + c * exp(-d * r);
    }
    else
    {
        exit_fail("ERROR WITH Vr_r");
    }
    return 0;
}

double L_J_dV_dr_init(double r, int a, int b)
{
    //For now we just have a, not b
    double A = L_J_params[a - 1][0];
    double B = L_J_params[a - 1][1];
    return A / pow(r, 7) - B / pow(r, 13);
}

void print_pots()
{
    dr_min = settings.DR_MIN_TAB;
    r_max = settings.R_MAX;
    if (dr_min == 0)
        return;
    if (n_rmax == -1)
        n_rmax = r_max / dr_min;

    debug_file << "Saving Potentials" << std::endl;

    std::ofstream pots_file;
    std::string filename = settings.output_name + "_generated.pots";
    pots_file.open(filename);
    // Print out a header, as well as info on ranges in the file
    pots_file << "# SAFARI Potentials and Forces File\n#\n# Table Range in file (min, step, max): " << dr_min << "-" << dr_min << "-" << r_max << "\n#\n# Atom1\tAtom2\tV_r\t-dV_dr\n" << std::flush;
    //Start at 1, as these guys are not defined for r=0 anyway.
    for (int i = 1; i < n_rmax; i++)
    {
        double r = dr_min * i;
        for (int n = 0; n < settings.NTYPES; n++)
        {
            double F = -dVr_dr(r, 0, n + 1);
            double V = Vr_r(r, 0, n + 1);
            Atom *a = settings.ATOMS[n];
            char buffer[200];
            sprintf(buffer, "%s\t%s\t%.5f\t%.5f\n",
                    settings.ion.symbol.c_str(), a->symbol.c_str(), V, F);
            pots_file << buffer << std::flush;
            if (settings.useLennardJones)
            {
                F = -dVr_dr(r, n + 1, n + 1);
                V = Vr_r(r, n + 1, n + 1);
                sprintf(buffer, "%s\t%s\t%.5f\t%.5f\n",
                        a->symbol.c_str(), a->symbol.c_str(), V, F);
                pots_file << buffer << std::flush;
            }
        }
    }
    pots_file.close();
}

void init_potentials()
{
    // Initialize these here, as they are needed even if we load in caches
    dr_min = settings.DR_MIN_TAB;
    r_max = settings.R_MAX;
    if (dr_min == 0)
        return;
    if (n_rmax == -1)
        n_rmax = r_max / dr_min;

    // Already was initialized earlier, we don't want to do this again.
    if (Vr_r_cache != NULL)
    {
        return;
    }

    //Initialize Lennard Jones potentials here
    if (settings.useLennardJones)
    {
        debug_file << "Initializing Lennard Jones Potentials" << std::endl;
        //We should have been given 2 parameters for each lattice atom.
        //TODO later swap this to in the atoms, so we have all pairs
        int start = settings.lattice_potential_start;
        L_J_dV_dr_cache = new double *[settings.NTYPES];
        L_J_params = new double *[settings.NTYPES];
        //TODO later, account for 2 params for all pairs of atoms somehow.
        for (int n = 0; n < settings.NTYPES; n++)
        {
            L_J_dV_dr_cache[n] = new double[n_rmax];
            L_J_dV_dr_cache[n][0] = 0;
            L_J_params[n] = new double[2];
            // Epsilon
            double epsilon = settings.binary_potential_parameters[start + 2 * n + 0];
            // Sigma
            double sigma = settings.binary_potential_parameters[start + 2 * n + 1];
            L_J_params[n][0] = -24 * epsilon * pow(sigma, 6);
            L_J_params[n][1] = -48 * epsilon * pow(sigma, 12);
            debug_file << "A: " << L_J_params[n][0];
            debug_file << ", B: " << L_J_params[n][1] << std::endl;
        }
        debug_file << "Initialized Lennard Jones Potentials" << std::endl;
    }

    Vr_r_cache = new double *[settings.NTYPES];
    dVr_dr_cache = new double *[settings.NTYPES];

    for (int n = 0; n < settings.NTYPES; n++)
    {
        Vr_r_cache[n] = new double[n_rmax];
        Vr_r_cache[n][0] = 0;
        dVr_dr_cache[n] = new double[n_rmax];
        dVr_dr_cache[n][0] = 0;
    }

    //Start at 1, as these guys are not defined for r=0 anyway.
    for (int i = 1; i < n_rmax; i++)
    {
        double r = dr_min * i;
        for (int n = 0; n < settings.NTYPES; n++)
        {
            dVr_dr_cache[n][i] = dVr_dr_init(r, n + 1);
            Vr_r_cache[n][i] = Vr_r_init(r, n + 1);
            if (settings.useLennardJones)
            {
                L_J_dV_dr_cache[n][i] = L_J_dV_dr_init(r, n + 1, n + 1);
            }
        }
    }
    print_pots();
}

double interp_r(double r, double *arr)
{
    if (r >= r_max)
        return 0;
    //Index for before r
    int i_bef = (int)(r / dr_min);
    //Index for after r
    int i_aft = i_bef + 1;
    //Anything out of bounds of table is 0.
    if (i_aft >= n_rmax)
        return 0;
    //Fraction of way between before and after
    double frac = (r - dr_min * i_bef) / dr_min;
    //Value before
    double bef = arr[i_bef];
    //Value after
    double aft = arr[i_aft];
    //Linearly interpolate between the two.
    return (1.0 - frac) * bef + frac * aft;
}

double Vr_r(double r, int n)
{
    if (dr_min == 0)
        return Vr_r_init(r, n);
    return interp_r(r, Vr_r_cache[n - 1]);
}

double dVr_dr(double r, int n)
{
    if (dr_min == 0)
        return dVr_dr_init(r, n);
    return interp_r(r, dVr_dr_cache[n - 1]);
}

double Vr_r(double r, int a, int b)
{
    if (Vr_r_all_cache == NULL)
        return Vr_r(r, b);
    return interp_r(r, Vr_r_all_cache[a + num_atoms * b]->data());
}

double dVr_dr(double r, int a, int b)
{
    if (dVr_dr_all_cache == NULL)
        return dVr_dr(r, b);
    return interp_r(r, dVr_dr_all_cache[a + num_atoms * b]->data());
}

double L_J_dV_dr(double r, int a, int b)
{
    if (dVr_dr_all_cache != NULL)
        return dVr_dr(r, a, b);

    if (dr_min == 0)
        return L_J_dV_dr_init(r, a, b);
    //TODO accound for b.
    return interp_r(r, L_J_dV_dr_cache[a - 1]);
}

double Vi_z(double z, int q)
{
    if (settings.image_potential_type == 0)
    {
        //This shouldn't be used for below the surface.
        if (z < settings.Z1)
            return 0;

        //Type 0: Only apply to the entry and exit trajectories.
        //Does not actually do anything during the numerical integration
        double z_min = settings.image_parameters[0];
        double v_min = settings.image_parameters[1];
        if (z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25 * eqsr * q;
            double eq_v = eq / v_min;
            return -eq / sqrt(dz * dz + eq_v * eq_v);
        }
        else
        {
            return -q * v_min;
        }
    }
    else if (settings.image_potential_type == 1)
    {
        //Type 1: Saturated image potential
        //Saturates based on the given z_min and v_min
        double z_min = settings.image_parameters[0];
        double v_min = settings.image_parameters[1];
        if (z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25 * eqsr * q;
            double eq_v = eq / v_min;
            return -eq / sqrt(dz * dz + eq_v * eq_v);
        }
        else
        {
            return -q * v_min;
        }
    }
    else
    {
        exit_fail("ERROR WITH Vi_z");
    }
    return 0;
}

double dVi_dz(double z, int q)
{
    if (settings.image_potential_type == 0)
    {
        //Type 0: Only apply to the entry and exit trajectories.
        //Does not actually do anything during the numerical integration
        return 0;
    }
    if (settings.image_potential_type == 1)
    {
        double z_min = settings.image_parameters[0];
        double v_min = settings.image_parameters[1];
        if (z > z_min)
        {
            double dz = z - z_min;
            double eq = 0.25 * eqsr * q;
            double eq_v = eq / v_min;
            return eq * dz / pow(dz * dz + eq_v * eq_v, 1.5);
        }
        return 0;
    }
    else
    {
        exit_fail("ERROR WITH dVi_dz");
    }
    return 0;
}

struct Pots
{
    std::string A;
    std::map<std::string, std::vector<double> *> V;
    std::map<std::string, std::vector<double> *> F;
};

void Atom::init_pots(std::string &filename)
{
    debug_file << "Loading provided potential tables" << std::endl;
    std::cout << "Loading provided potential tables" << std::endl;

    std::ifstream input;
    filename = filename + ".pots";
    input.open(filename);
    debug_file << "Opened File "<< filename << std::endl;

    std::map<std::string, Atom *> atoms;
    std::map<std::string, Pots *> pots;

    // Initialize the maps

    // Add the ion first
    atoms[settings.ion.symbol] = &settings.ion;
    pots[settings.ion.symbol] = new Pots();
    pots[settings.ion.symbol]->A = settings.ion.symbol;
    
    for (auto a : settings.ATOMS)
    {
        atoms[a->symbol] = a;
        pots[a->symbol] = new Pots();
        pots[a->symbol]->A = a->symbol;
    }

    // Cleanup the maps so no duplicates (ie A-B = B-A)
    for(auto const& a: pots)
    {
        Pots *potA = pots[a.first];
        for(auto const& b: pots)
        {
            Pots *potB = pots[b.first];
            if(!potA->V[b.first])
            {
                potA->V[b.first] = new std::vector<double>();
                pots[a.first]->F[b.first] = new std::vector<double>();

                // Put a 0 for the first entry here,
                // our first distance index is 1, as these functions
                // are not defined at 0 distance anyway
                pots[a.first]->V[b.first]->push_back(0);
                pots[a.first]->F[b.first]->push_back(0);
            }
            potB->V[potA->A] = potA->V[b.first];
            potB->F[potA->A] = potA->F[b.first];
        }
    }

    if (input.is_open())
    {
        std::string line;
        debug_file << "Reading Potentials" << std::endl;
        while (getline(input, line))
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

            //Splits line for parsing.
            std::vector<std::string> line_args = split(line);

            // First argument is the atom A, second atom B
            std::string A = line_args[0];
            std::string B = line_args[1];

            Pots *pot = pots[A];
            pot->V[B]->push_back(atof(line_args[2].c_str()));
            pot->F[B]->push_back(-atof(line_args[3].c_str()));
        }
    }
    // Arrays are now populated with some tables, ideally they are ordered such that the entire
    // table for each atom is set up correctly, and no redundant tables, ie only one table
    // for each pair.
    debug_file << "Processing Potentials" << std::endl;

    num_atoms = settings.ATOMS.size() + 1;

    Vr_r_cache = new double *[0];
    dVr_dr_cache = new double *[0];

    Vr_r_all_cache = new std::vector<double> *[num_atoms * num_atoms];
    dVr_dr_all_cache = new std::vector<double> *[num_atoms * num_atoms];

    int m = 0;

    // Starts at 1, as index 0 is the ion
    int n = 1;
    Pots *potI = pots[settings.ion.symbol];

    Vr_r_all_cache[0] = potI->V[settings.ion.symbol];
    dVr_dr_all_cache[0] = potI->F[settings.ion.symbol];

    for (auto b : settings.ATOMS)
    {
        Vr_r_all_cache[m + n * num_atoms] = potI->V[b->symbol];
        dVr_dr_all_cache[m + n * num_atoms] = potI->F[b->symbol];
        Vr_r_all_cache[n + m * num_atoms] = potI->V[b->symbol];
        dVr_dr_all_cache[n + m * num_atoms] = potI->F[b->symbol];
        n++;
    }
    m = 1;
    
    for (auto a : settings.ATOMS)
    {
        Pots *potA = pots[a->symbol];
        // Starts at 1, as index 0 is the ion
        n = 1;
        for (auto b : settings.ATOMS)
        {
            Vr_r_all_cache[m + n * num_atoms] = potA->V[b->symbol];
            dVr_dr_all_cache[m + n * num_atoms] = potA->F[b->symbol];

            Vr_r_all_cache[n + m * num_atoms] = potA->V[b->symbol];
            dVr_dr_all_cache[n + m * num_atoms] = potA->F[b->symbol];
            n++;
        }
        m++;
    }
    print_pots();
}

double electron_density(Lattice *lattice, Ion &ion, bool predicted)
{
    //TODO see https://github.com/SINS-Lab/SAFARI/blob/10305e6f9ee597e89a6df7acfede3554371f42bc/src/hameqinel.f
    // it has calculations for electron density.
    double z0 = 0.755;
    double z = predicted ? ion.r_t[2] : ion.r[2];
    return z < 0.2 ? 0.087 : 0.1104 * std::exp(-z / z0);
}

double apply_friction(Lattice *lattice, Ion &ion, double *F, double dt, bool predicted)
{
    double vx = ion.p[0] * ion.atom->mass_inv;
    double vy = ion.p[1] * ion.atom->mass_inv;
    double vz = ion.p[2] * ion.atom->mass_inv;

    double v_sq = vx * vx + vy * vy + vz * vz;
    double v = sqrt(v_sq);
    //No velocity, no friction.
    if (v == 0)
        return 0;

    //Components of the velocity direction.
    vx /= v;
    vy /= v;
    vz /= v;

    double density = electron_density(lattice, ion, predicted);
    //No electron density, no friction
    if (density == 0)
        return 0;

    //assume magnitude of friction is:
    //Av + Bv^2, where v is magnitude of velocity.
    //This is then scaled by electron density for the location
    double fric = (settings.F_a * v + settings.F_b * v_sq) * density;
    //Direction of friction is opposite of velocity.
    F[0] -= vx * fric;
    F[1] -= vy * fric;
    F[2] -= vz * fric;

    //TODO return effective "potential" from this interaction.
    return 0;
}