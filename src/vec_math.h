#ifndef VEC_MATH_H_INCLUDED
#define VEC_MATH_H_INCLUDED

class vec3d
{
public:
    double v[3];
    double norm_sq();
    double norm();
    void scale(double s);
    double operator*(vec3d b);
    vec3d operator*(double b);
    vec3d operator/(double b);
    vec3d normalize();
    vec3d cross(vec3d b);
    vec3d operator+(vec3d b);
    vec3d operator-(vec3d b);
    vec3d operator-(double b[]);
    void set(double x, double y, double z);
    void set(double arr[]);
    double& operator[](int i)
    {
        return v[i];
    }
};

class mat3d
{
public:
    double m[9];
    void identity();
    mat3d operator+(mat3d b);
    mat3d operator*(double d);
    mat3d operator*(mat3d d);
    vec3d operator*(vec3d x);
    double& operator[](int i)
    {
        return m[i];
    }
    double& operator()(int row, int col)
    {
        return m[(row*3) + (col%3)];
    }
    mat3d invert();
    mat3d dot(mat3d m);
};


mat3d make_rot_matrix(vec3d direction, vec3d axis);

int to_hash(double x, double y, double z);

void print(mat3d R, char* header);
void print(vec3d V, char* header);

#endif // VEC_MATH_H_INCLUDED
