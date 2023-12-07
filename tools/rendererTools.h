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

    Vector3 operator+(Vector3 const& obj)
    {
        return Vector3(x + obj.x, y + obj.y, z + obj.z);
    }
    Vector3 operator-(Vector3 const& obj)
    {
        return Vector3(x - obj.x, y - obj.y, z - obj.z);
    }
    Vector3 operator*(float const& obj)
    {
        return Vector3(x * obj, y * obj, z * obj);
    }
    Vector3 operator/(float const& obj)
    {
        return Vector3(x / obj, y / obj, z / obj);
    }

private:
};

Vector3 scalarProd(Vector3 a, Vector3 b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 crossProd(Vector3 a, Vector3 b){
    return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, b.y * a.x - a.y * b.x);
}

float det3x3(Matrix a){
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[2][0]*a[1][1]*a[0][2] - a[2][1]*a[1][2]*a[0][0] - a[2][2]*a[1][0]*a[0][1];
}