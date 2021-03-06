#include "xyz.h"
#include "vec_math.h"
#include "string_utils.h"
#include <cstdio>
#include <string>

#define THREADCOUNT 6

class Particle
{
public:
    std::string atom;
    Vec3d position;
    Vec3d momentum;
    Vec3d velocity;
    double mass;
    int id;
    int nearby;

    void fromXYZ(XYZ_Single input, int index)
    {
        atom = input.atoms[index];
        double* vals = input.values[index];
        position.set(vals[0],vals[1],vals[2]);
        momentum.set(vals[3],vals[4],vals[5]);
        mass = vals[6];
        velocity = momentum/mass;
        id = vals[7];
        nearby = vals[8];
    }
};

std::string to_string(double num)
{
    char* arr = new char[8];
    sprintf(arr, "%.4f", num);
    std::string output = arr;
    return output;
}

XYZ_Single from_Particles(double time, std::vector<Particle> pset)
{
    XYZ_Single xyz_single;
    xyz_single.number = pset.size();
    xyz_single.comment = to_string(time);
    xyz_single.max_row_index = 9;
    xyz_single.atoms = new std::string[xyz_single.number];
    xyz_single.values = new double*[xyz_single.number];
    #pragma omp parallel for num_threads(THREADCOUNT)
    for(int i = 0; i<xyz_single.number; i++)
    {
        Particle &p = pset[i];
        xyz_single.values[i] = new double[9];
        xyz_single.atoms[i] = p.atom;
        xyz_single.values[i][0] = p.position[0];
        xyz_single.values[i][1] = p.position[1];
        xyz_single.values[i][2] = p.position[2];
        xyz_single.values[i][3] = p.momentum[0];
        xyz_single.values[i][4] = p.momentum[1];
        xyz_single.values[i][5] = p.momentum[2];
        xyz_single.values[i][6] = p.mass;
        xyz_single.values[i][7] = p.id;
        xyz_single.values[i][8] = p.nearby;
    }
    return xyz_single;
}

Vec3d linear_interpolate(double t_rel, double dt, Vec3d& prev, Vec3d& next)
{
    Vec3d interp;
    for(int i = 0; i<3; i++)
    {
        double dx_dt = (next[i] - prev[i]) / dt;
        double dx = dx_dt * t_rel;
        interp[i] = prev[i] + dx;
    }
    return interp;
}

std::vector<Particle> interpolate_states(double time, std::vector<double> times, 
                                std::vector<std::vector<Particle>> particles)
{
    int next_index = times.size()-1;
    //Find the index for the time after our current time.
    for(int i = 0; i<=next_index; i++)
    {
        double t = times[i];
        if(t > time) 
        {
            next_index = i;
            break;
        }
    }
    //If next time is 0, return first frame.
    if(next_index == 0) return particles[0];
    std::vector<Particle> new_particles;
    int prev_index = next_index - 1;
    double t_prev = times[prev_index];
    double t_next = times[next_index];

    double dt = t_next - t_prev;
    double t_rel = time - t_prev;

    int num = particles[0].size();
    new_particles.resize(num);
    //Otherwise, lets interpolate
    #pragma omp parallel for num_threads(THREADCOUNT)
    for(int i = 0; i<num; i++)
    {
        Particle& prev = particles[prev_index][i];
        Particle& next = particles[next_index][i];

        //Copy previous into interp
        Particle interp = prev;
        interp.position = linear_interpolate(t_rel, dt, prev.position, next.position);
        interp.momentum = linear_interpolate(t_rel, dt, prev.momentum, next.momentum);
        interp.velocity = linear_interpolate(t_rel, dt, prev.velocity, next.velocity);
        new_particles[i] = interp;
    }
    return new_particles;
}

void apply_colour(std::vector<Particle>& pset, std::string& colour)
{
    int num = pset.size();
    if(colour == "velocity")
    {
        //TODO make this instead by some multiple of KT from a safio file.
        double min_E = 0.05;
        for(int i = 0; i<num; i++)
        {
            Particle&p = pset[i];
            //Ignore the ion
            if(p.id==0) continue;
            double E = (p.momentum[0]*p.momentum[0]
                       +p.momentum[1]*p.momentum[1]
                       +p.momentum[2]*p.momentum[2])
                       /(2*p.mass);
            if(E > min_E) p.atom = "X";
        }
    }
    else if(colour == "nearby")
    {
        //Just flag it as an X if it is nearby.
        for(int i = 0; i<num; i++)
        {
            Particle&p = pset[i];
            //Ignore the ion
            if(p.id==0) continue;
            if(p.nearby) p.atom = "X";
        }
    }
    else if(colour!="")
    {
        std::cout << "Colour should be nearby or velocity" << std::endl;
        exit(EXIT_FAILURE);
    }
}

XYZ smooth(const XYZ& original, std::string colour)
{
    //first, lets convert the xyz into a 2d vector of particles.
    std::vector<std::vector<Particle>> particles;
    //Also store vector of the times.
    std::vector<double> times;
    for(XYZ_Single single: original.xyzs)
    {
        std::vector<Particle> pset;
        double* timearr = to_double_array(split(single.comment), 0, 0);
        times.push_back(timearr[0]);
        for(int i = 0; i<single.number; i++)
        {
            Particle p;
            p.fromXYZ(single, i);
            pset.push_back(p);
        }
        particles.push_back(pset);
    }
    //End of our time step is whatever the last time was.
    double end_time = times[times.size()-1];
    
    //Now we construct a new xyz in even time steps
    XYZ new_xyz;

    int SIZE = std::max(1000, (int)(original.xyzs.size() + 5));

    new_xyz.xyzs.resize(SIZE);
    //Note that we are doing this with a fixed 1000 frames.
    #pragma omp parallel for num_threads(THREADCOUNT)
    for(int i = 0; i<SIZE; i++)
    {
        double time = i * end_time/SIZE;
        std::vector<Particle> pset = interpolate_states(time, times, particles);
        apply_colour(pset, colour);
        XYZ_Single single = from_Particles(time, pset);
        new_xyz.xyzs[i] = single;
    }
    return new_xyz;
}

int main(int argc,char* argv[])
{
    std::map<std::string, ArgValue> args = get_arguments(argc, argv);
    std::string safio_file = args["-i"].as_string();

    std::string in_file = args["-i"].as_string();
    std::string out_file = args["-o"].as_string();
    std::string colour = args["-c"].as_string();

    if(in_file == "" || out_file == "") 
    {
        std::cout << "Usage: -i [inputfile] -o [outputfile] | optional: -c [nearest|velocity]" << std::endl;
        return 1;
    }

    XYZ xyz;
    std::ifstream input;
    std::cout << "Opening input" << std::endl;
    input.open(in_file);
    std::cout << "Loading xyz" << std::endl;
    xyz.load(input);
    std::cout << "Closing input, frames loaded: " <<xyz.xyzs.size() << std::endl;
    input.close();

    std::ofstream output;
    std::cout << "Opening output" << std::endl;
    output.open(out_file);
    std::cout << "Smoothing xyz" << std::endl;
    xyz = smooth(xyz, colour);
    std::cout << "Saving xyz, frames saved: " <<xyz.xyzs.size() << std::endl;
    xyz.save(output);
    std::cout << "Closing output" << std::endl;
    output.close();

    return 0;
}