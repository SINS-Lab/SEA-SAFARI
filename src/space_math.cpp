#include "space_math.h"
#include <cmath>

#include "safio.h"
#include <iostream>

double cube_root(double val)
{
    return pow(val, 1.0/3.0);
}

/**
 * This fills the location based on the given shell and index on the current shell.
 * 
 * diffSq is the difference in areas between this shell and previous.
 * diffCb is the difference in volumes between this shell and previous.
 */
void index_to_loc(int radius, int index, int diffSq, int diffCb, Vec3d &location)
{
    location[0] = 0;
    location[1] = 0; 
    location[2] = 0;
    int layerSize = (2 * radius + 1) * (2 * radius + 1);
    if (index == 0)
    {
        location[0] = -radius;
        location[1] = radius;
        location[2] = -radius;
        return;
    }
    // Fill z
    {
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
    }
    // Fill x
    if (!(location[2  ] == radius || location[2  ] == -radius))
    {
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
    }
    else
    {
        int temp = (index % layerSize);
        temp = temp % (2 * radius + 1);
        temp -= radius;
        location[0] = temp;
    }
    // Fill y
    if (!(location[2  ] == radius || location[2  ] == -radius))
    {
        int temp = (index) % diffSq;
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
        int temp = (index % layerSize) / (2 * radius + 1);
        temp -= radius;
        location[1] = temp;
    }
}

/**
 * This converts the index given to a location in space.
 * This is done such that the location given increases in
 * radius from 0, as cubic shells.
 */
void index_to_loc(int index, Vec3d &location)
{
    if (index > 0)
    {
        //Radius of cube we are on.
        int radius = floor(ceil(cube_root(index)/2.0));

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