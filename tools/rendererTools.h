#pragma once
#include <vector>
#include "Classes/RendererClasses.h"

typedef std::vector<std::vector<float>> Matrix;

class Vector3
{
public:
    Vector3(float _x = 0, float _y = 0, float _z = 0)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    ~Vector3() {}

    float x;
    float y;
    float z;

    void print()
    {
        std::cout << "(x:" << x << ", y:" << y << ", z:" << z << ")";
    }
    float magnitude()
    {
        return sqrt(x*x + y*y + z*z);
    }
    Vector3 normalized()
    {
        float mag = sqrt(x*x + y*y + z*z);
        Vector3 ret(x / mag, y / mag, z / mag);
        return ret;
    }

private:
};

Vector3 addVec(Vector3 a, Vector3 b)
{
    Vector3 ret(a.x + b.x, a.y + b.y, a.z + b.z);
    return ret;
}

Vector3 subVec(Vector3 a, Vector3 b)
{
    Vector3 ret(a.x - b.x, a.y - b.y, a.z - b.z);
    return ret;
}

Vector3 mulVec(Vector3 a, float b)
{
    Vector3 ret(a.x * b, a.y * b, a.z * b);
    return ret;
}

Vector3 divVec(Vector3 a, float b)
{
    Vector3 ret(a.x / b, a.y / b, a.z / b);
    return ret;
}

Vector3 scalarProd(Vector3 a, Vector3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 crossProd(Vector3 a, Vector3 b)
{
    Vector3 ret(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, b.y * a.x - a.y * b.x);
    return ret;
}

float det3x3(Matrix a)
{
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[2][0]*a[1][1]*a[0][2] - a[2][1]*a[1][2]*a[0][0] - a[2][2]*a[1][0]*a[0][1];
}