#include "vec_math.h"
#include "safio.h"
#include <vector>
#include <math.h>

double sqr(double *V)
{
    return V[0]*V[0]+V[1]*V[1]+V[2]*V[2];
}

void Mat3d::identity()
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

Mat3d Mat3d::invert()
{
    Mat3d &A = *this;
    // computes the determinant of A
    double det = A(0,0) * (A(1,1) * A(2,2) - A(2,1) * A(1,2)) -
                 A(0,1) * (A(1,0) * A(2,2) - A(1,2) * A(2,0)) +
                 A(0,2) * (A(1,0) * A(2,1) - A(1,1) * A(2,0));

    double invdet = 1 / det;

    Mat3d minv; // inverse of matrix A
    minv(0,0) = (A(1,1) * A(2,2) - A(2,1) * A(1,2)) * invdet;
    minv(0,1) = (A(0,2) * A(2,1) - A(0,1) * A(2,2)) * invdet;
    minv(0,2) = (A(0,1) * A(1,2) - A(0,2) * A(1,1)) * invdet;
    minv(1,0) = (A(1,2) * A(2,0) - A(1,0) * A(2,2)) * invdet;
    minv(1,1) = (A(0,0) * A(2,2) - A(0,2) * A(2,0)) * invdet;
    minv(1,2) = (A(1,0) * A(0,2) - A(0,0) * A(1,2)) * invdet;
    minv(2,0) = (A(1,0) * A(2,1) - A(2,0) * A(1,1)) * invdet;
    minv(2,1) = (A(2,0) * A(0,1) - A(0,0) * A(2,1)) * invdet;
    minv(2,2) = (A(0,0) * A(1,1) - A(1,0) * A(0,1)) * invdet;

    return minv;
}

Mat3d Mat3d::operator*(double d)
{
    Mat3d M;
    for(int i = 0; i<9; i++)
    {
        M[i] = m[i] * d;
    }
    return M;
}

Mat3d Mat3d::operator+(Mat3d d)
{
    Mat3d M;
    for(int i = 0; i<9; i++)
    {
        M[i] = m[i] + d[i];
    }
    return M;
}

Mat3d Mat3d::operator*(Mat3d B)
{
    Mat3d R;
    Mat3d &A = *this;
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

Vec3d Mat3d::operator*(Vec3d v)
{
    Vec3d x;
    x[0] = v[0]*m[0] + v[1]*m[1] + v[2]*m[2];
    x[1] = v[0]*m[3] + v[1]*m[4] + v[2]*m[5];
    x[2] = v[0]*m[6] + v[1]*m[7] + v[2]*m[8];
    return x;
}

double Vec3d::norm_sq()
{
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

double Vec3d::norm()
{
    return sqrt(norm_sq());
}

Vec3d Vec3d::normalize()
{
    Vec3d n;
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

double Vec3d::operator*(Vec3d b)
{
    return v[0]*b[0] + v[1]*b[1] + v[2]*b[2];
}

Vec3d Vec3d::operator+(Vec3d b)
{
    Vec3d sum;
    sum[0] = v[0] + b[0];
    sum[1] = v[1] + b[1];
    sum[2] = v[2] + b[2];
    return sum;
}

Vec3d Vec3d::operator-(Vec3d b)
{
    Vec3d diff;
    diff[0] = v[0] - b[0];
    diff[1] = v[1] - b[1];
    diff[2] = v[2] - b[2];
    return diff;
}

Vec3d Vec3d::operator-(double b[])
{
    Vec3d diff;
    diff[0] = v[0] - b[0];
    diff[1] = v[1] - b[1];
    diff[2] = v[2] - b[2];
    return diff;
}

Vec3d Vec3d::operator/(double b)
{
    Vec3d mult;
    mult[0] = v[0]/b;
    mult[1] = v[1]/b;
    mult[2] = v[2]/b;
    return mult;
}

Vec3d Vec3d::operator*(double b)
{
    Vec3d mult;
    mult[0] = v[0]*b;
    mult[1] = v[1]*b;
    mult[2] = v[2]*b;
    return mult;
}

void Vec3d::set(double x, double y, double z)
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

void Vec3d::set(double x[])
{
    v[0] = x[0];
    v[1] = x[1];
    v[2] = x[2];
}

void Vec3d::scale(double x)
{
    v[0] *= x;
    v[1] *= x;
    v[2] *= x;
}

Vec3d Vec3d::cross(Vec3d b)
{
    Vec3d c;
    c[0] = v[1]*b[2] - v[2]*b[1];
    c[1] = v[2]*b[0] - v[0]*b[2];
    c[2] = v[0]*b[1] - v[1]*b[0];
    return c;
}