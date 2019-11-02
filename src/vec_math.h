#pragma once
#define M_PI           3.14159265358979323846  /* pi */

class Vec3d
{
public:
    double v[3];
    double norm_sq();
    double norm();
    void scale(double s);
    double operator*(Vec3d b);
    Vec3d operator*(double b);
    Vec3d operator/(double b);
    Vec3d &operator*=(const double &b);
    Vec3d &operator/=(const double &b);
    Vec3d normalize();
    Vec3d cross(Vec3d b);
    Vec3d operator+(Vec3d b);
    Vec3d operator+(double* b);
    Vec3d operator-(Vec3d b);
    Vec3d &operator+=(const Vec3d &b);
    Vec3d &operator+=(const double* b);
    Vec3d &operator-=(const Vec3d &b);
    Vec3d operator-(double b[]);
    void set(double x, double y, double z);
    void set(double arr[]);
    double& operator[](int i)
    {
        return v[i];
    }
};

class Mat3d
{
public:
    double m[9];
    void identity();
    Mat3d operator+(Mat3d b);
    Mat3d operator*(double d);
    Mat3d operator*(Mat3d d);
    Vec3d operator*(Vec3d x);
    double& operator[](int i)
    {
        return m[i];
    }
    double& operator()(int row, int col)
    {
        return m[(row*3) + (col%3)];
    }
    Mat3d invert();
    Mat3d dot(Mat3d m);
};

/**
 * Squares the given array, assuming it is 3x1
 */
double sqr(double*V);
double sqr(const double*V);

/**
 * Computes the squared difference of the given arrays
 * assuming they are 3x1
 */
double diff_sqr(double *X, double *Y);
double diff_sqr(const double *X,const double *Y);
