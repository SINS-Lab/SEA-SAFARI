#pragma once
#include "vec_math.h"

#define CELL_SIZE 5.0

#define N_CUBE_MASK 29791 //This is 31^3, radius 15

extern double space_mask_cube[N_CUBE_MASK][3];

//Populates space_mask.
void init_lookup();

struct Point
{
    double x;
    double y;
};

class Mask
{
public:
    Point points[256];
    int num = 0;
    bool inside(Point &p);
    bool inside(double x, double y)
    {
        Point p;
        p.x = x;
        p.y = y;
        return inside(p);
    }
};

/**
 * This fills the location based on the given
 * shell and index on the current shell.
 * 
 * @param diffSq - the difference in areas 
 * between this shell and previous.
 * 
 * @param diffCb - the difference in volumes 
 * between this shell and previous.
 */
void index_to_loc(int radius, int index, int diffSq, int diffCb, Vec3d &location);

/**
 * This converts the index given to a location in space.
 * This is done such that the location given increases in
 * radius from 0, as cubic shells.
 * 
 * The radius of these cubes is 1 unit, and, if using
 * radii less ((than N_CUBE_MASK^1/3)-1)/2, will be looked up in space_mask
 * 
 * @param index - index of this cell, 0 being the origin
 * @param location - vector to fill with the location.
 */
void index_to_loc(int index, Vec3d &location);

/**
 * Makes a matrix which will rotate the given direction
 * onto the given axis.
 * 
 * @param direction - the reference direction
 * @param axis - the vector to rotate direction onto.
 */
Mat3d make_rot_matrix(Vec3d direction, Vec3d axis);

/**
 * This converts the given location into an index.
 * All locations within ion-axis-aligned cubes, of size
 * 5x5x5 A^3 will map to the same hash, which is
 * guarenteed to be > 0, so long as the location
 * is within 512A of the origin.
 * 
 * This allows for quick lookup of Cells of sites in 
 * an int-indexed map.
 */
int to_hash(double x, double y, double z);