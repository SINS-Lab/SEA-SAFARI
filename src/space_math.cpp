#include "space_math.h"
#include <cmath>

#include "safio.h"
#include <iostream>

//Note, this implementation is somehow rather slow, until
//Whatever -O3 does to it. Maybe it should be looked into.
int to_hash(double x, double y, double z)
{
    int i = (int)(x/CELL_SIZE + 512);
    int j = (int)(y/CELL_SIZE + 512);
    int k = (int)(z/CELL_SIZE + 512);
    return i + (j << 10) + (k << 20);
}

Mat3d make_rot_matrix(Vec3d a, Vec3d b)
{
    Mat3d I;
    I.identity();
    a = a.normalize();
    b = b.normalize();

    Vec3d v = a.cross(b);
    double c = a*b;
    c = 1.0/(1.0+c);

    Mat3d V;

    V[0] = 0;
    V[1] = -v[2];
    V[2] = v[1];

    V[3] = v[2];
    V[4] = 0;
    V[5] = -v[0];

    V[6] = -v[1];
    V[7] = v[0];
    V[8] = 0;

    return I + V + V*V*c;
}

void index_to_loc(int radius, int index, int diffSq, int diffCb, Vec3d &location)
{
    location[0] = 0;
    location[1] = 0; 
    location[2] = 0;
    int layerSize = (2 * radius + 1) * (2 * radius + 1);
    if (index == 0)
    {
        location[0] = -radius;
        location[1] = -radius;
        location[2] = radius;
        return;
    }
    // Fill z
    if (index >= layerSize && index < diffCb - layerSize)
    {
        int temp = (index - layerSize) / (diffSq) + 1;
        temp -= radius;
        temp = temp > radius ? radius : temp < -radius ? -radius : temp;
        location[2] = temp;
    }
    else
    {
        location[2] = index > layerSize ? -radius : radius;
    }
    
    // Fill x and y
    if (!(location[2] == radius || location[2] == -radius))
    {
        //Fill x
        int temp = (index) % diffSq;
        if (temp < diffSq / 2)
        {
            if (temp < radius)
            {
                location[0] = temp;
            }
            else if (temp > diffSq / 2 - radius)
            {
                location[0] = -(temp - (diffSq / 2));
            }
            else location[0] = radius;
        }
        else if (temp > diffSq / 2)
        {
            temp -= diffSq / 2;
            if (temp < radius)
            {
                location[0] = -temp;
            }
            else if (temp > diffSq / 2 - radius)
            {
                location[0] = (temp - (diffSq / 2));
            }
            else location[0] = -radius;
        }

        //Fill y
        temp = (index) % diffSq;
        temp = (temp + 2 * radius - 1) % diffSq + 1;
        if (temp < diffSq / 2)
        {
            if (temp < radius)
            {
                location[1] = temp;
            }
            else if (temp > diffSq / 2 - radius)
            {
                location[1] = -(temp - (diffSq / 2));
            }
            else location[1] = radius;
        }
        else if (temp > diffSq / 2)
        {
            temp -= diffSq / 2;
            if (temp < radius)
            {
                location[1] = -temp;
            }
            else if (temp > diffSq / 2 - radius)
            {
                location[1] = (temp - (diffSq / 2));
            }
            else location[1] = -radius;
        }
    }
    else
    {
        int temp = (index % layerSize);
        temp = temp % (2 * radius + 1);
        temp -= radius;
        location[0] = temp;

        temp = (index % layerSize) / (2 * radius + 1);
        temp -= radius;
        location[1] = temp;
    }
}

//If not populated, do so.
bool populated = false;

void init_lookup()
{
    populated = false;
    for(int i = 0; i<3375; i++)
    {
        Vec3d loc;
        index_to_loc(i, loc);

        space_mask[i][0] = loc[0];
        space_mask[i][1] = loc[1];
        space_mask[i][2] = loc[2];
    }
    populated = true;
}

void index_to_loc(int index, Vec3d &location)
{
    if(populated && index < 3375)
    {
        location.set(space_mask[index]);
    }
    else if (index > 0)
    {
        double ind = index;
        //Radius of cube we are on.
        int radius = floor(ceil(cbrt(ind))/2.0);
        //This is a special case where this formula doesn't work
        if(index==1) radius = 1;

        //Area of a face of the current cube.
        int current_area = (2*radius + 1) *  (2*radius + 1); 
        //Total volume of current cube.
        int current_volume = current_area * (2*radius + 1);

        //previous cube's radius
        int radius_p = radius - 1;
        //Area of a face of the previous cube.
        int prev_area = (2*radius_p + 1) *  (2*radius_p + 1);
        //Volume of previous cube.
        int prev_volume = prev_area * (2*radius_p + 1);

        //Fill from here.
        index_to_loc(radius, index - prev_volume, 
                             current_area - prev_area, 
                             current_volume - prev_volume, location);
    }
    else
    {
        //Origin is 0
        location[0] = 0;
        location[1] = 0; 
        location[2] = 0;
    }
}