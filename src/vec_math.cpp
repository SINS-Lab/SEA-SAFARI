#include "vec_math.h"
#include "safio.h"
#include <vector>
#include <math.h>

int to_hash(double x, double y, double z)
{
    int i = (int)(x/settings.AX + 512);
    int j = (int)(y/settings.AZ + 512);
    int k = (int)(z/settings.AY + 512);
    return i + (j << 10) + (k << 20);
}

void mat3d::identity()
{
    /*
    m = {1,0,0,
         0,1,0,
         0,0,1};*/
    m[0] = 1;
    m[1] = 0;
    m[2] = 0;

    m[3] = 0;
    m[4] = 1;
    m[5] = 0;

    m[6] = 0;
    m[7] = 0;
    m[8] = 1;
}

mat3d mat3d::invert()
{
    mat3d A = *this;
    // computes the inverse of a matrix m
    double det = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
                 A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
                 A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));

    double invdet = 1 / det;

    mat3d minv; // inverse of matrix m
    minv(0, 0) = (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) * invdet;
    minv(0, 1) = (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2)) * invdet;
    minv(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * invdet;
    minv(1, 0) = (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)) * invdet;
    minv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * invdet;
    minv(1, 2) = (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) * invdet;
    minv(2, 0) = (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1)) * invdet;
    minv(2, 1) = (A(2, 0) * A(0, 1) - A(0, 0) * A(2, 1)) * invdet;
    minv(2, 2) = (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) * invdet;

    return minv;
}

mat3d mat3d::operator*(double d)
{
    mat3d M;
    for(int i = 0; i<9; i++)
    {
        M[i] = m[i] * d;
    }
    return M;
}

mat3d mat3d::operator+(mat3d d)
{
    mat3d M;
    for(int i = 0; i<9; i++)
    {
        M[i] = m[i] + d[i];
    }
    return M;
}

mat3d mat3d::operator*(mat3d B)
{
    mat3d R;
    mat3d A = *this;
    for(int i = 0; i<9; i++)
    {
        int r = i/3;
        int c = i%3;
        R(r,c) = A(r,0)*B(0,c)+
                 A(r,1)*B(1,c)+
                 A(r,2)*B(2,c);
    }
    return R;
}

vec3d mat3d::operator*(vec3d v)
{
    vec3d x;
    x[0] = v[0]*m[0] + v[1]*m[1] + v[2]*m[2];
    x[1] = v[0]*m[3] + v[1]*m[4] + v[2]*m[5];
    x[2] = v[0]*m[6] + v[1]*m[7] + v[2]*m[8];
    return x;
}

double vec3d::norm_sq()
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

double vec3d::norm()
{
    return sqrt(norm_sq());
}

vec3d vec3d::normalize()
{
    vec3d n;
    double length = norm();
    if(length>0)
    {
        n.set(v);
        n = n * (1.0d/length);
    }
    else
    {
        n.set(1,0,0);
    }
    return n;
}

double vec3d::operator*(vec3d b)
{
    return v[0]*b[0] + v[1]*b[1] + v[2]*b[2];
}

vec3d vec3d::operator+(vec3d b)
{
    vec3d sum;
    sum[0] = v[0] + b[0];
    sum[1] = v[1] + b[1];
    sum[2] = v[2] + b[2];
    return sum;
}

vec3d vec3d::operator-(vec3d b)
{
    vec3d diff;
    diff[0] = v[0] - b[0];
    diff[1] = v[1] - b[1];
    diff[2] = v[2] - b[2];
    return diff;
}

vec3d vec3d::operator-(double b[])
{
    vec3d diff;
    diff[0] = v[0] - b[0];
    diff[1] = v[1] - b[1];
    diff[2] = v[2] - b[2];
    return diff;
}

vec3d vec3d::operator/(double b)
{
    vec3d mult;
    mult[0] = v[0]/b;
    mult[1] = v[1]/b;
    mult[2] = v[2]/b;
    return mult;
}

vec3d vec3d::operator*(double b)
{
    vec3d mult;
    mult[0] = v[0]*b;
    mult[1] = v[1]*b;
    mult[2] = v[2]*b;
    return mult;
}

void vec3d::set(double x, double y, double z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

void vec3d::set(double x[])
{
    v[0] = x[0];
    v[1] = x[1];
    v[2] = x[2];
}

void vec3d::scale(double x)
{
    v[0] *= x;
    v[1] *= x;
    v[2] *= x;
}

vec3d vec3d::cross(vec3d b)
{
    vec3d c;
    c[0] = v[1]*b[2] - v[2]*b[1];
    c[1] = v[2]*b[0] - v[0]*b[2];
    c[2] = v[0]*b[1] - v[1]*b[0];
    return c;
}

mat3d make_rot_matrix(vec3d a, vec3d b)
{
    mat3d I;
    I.identity();
    a = a.normalize();
    b = b.normalize();

    vec3d v = a.cross(b);
    double c = a*b;
    c = 1.0d/(1.0d+c);

    mat3d V;

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
