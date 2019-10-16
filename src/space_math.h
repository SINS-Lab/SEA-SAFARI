#ifndef SPACE_MATH_H_INCLUDED
#define SPACE_MATH_H_INCLUDED
#include "vec_math.h"

double cube_root(double value);

void index_to_loc(int radius, int index, int diffSq, int diffCb, vec3d &location);

void index_to_loc(int index, vec3d &location);

#endif // SPACE_MATH_H_INCLUDED