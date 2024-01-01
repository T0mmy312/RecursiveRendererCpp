#pragma once

#ifndef _MATHANDVECTORTOOLS_H_
#define _MATHANDVECTORTOOLS_H_

#include <vector>
#include <iostream>
#include <math.h>

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Math and Vector3 tools
//? ----------------------------------------------------------------------------------------------------------------------------------

#define PI 3.14159265359f
#define DEG_TO_RAD 0.0174532925199f
#define RAD_TO_DEG 57.2957795131f

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
    float magnitude() // returns the lenght of the vector
    {
        return sqrt(x*x + y*y + z*z);
    }
    Vector3 normalized() // returns the normalized vector of this vector (magnitude = 1)
    {
        float mag = sqrt(x*x + y*y + z*z);
        return Vector3(x / mag, y / mag, z / mag);
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

float sqr(float x) { // just returns x^2
    return x * x;
}

float scalarProd(Vector3 a, Vector3 b) { // returns the scalar produkt of two vectors
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vector3 crossProd(Vector3 a, Vector3 b) { // returns the normal vector of two vectors
    return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, b.y * a.x - a.y * b.x);
}

float angleInRads(Vector3 a, Vector3 b) { // returns the angle between two vectors
    return acos(scalarProd(a, b) / (a.magnitude() * b.magnitude()));
}

float det3x3(Matrix a) { // returns the determinat of a 3x3 Matrix
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[2][0]*a[1][1]*a[0][2] - a[2][1]*a[1][2]*a[0][0] - a[2][2]*a[1][0]*a[0][1];
}

float clamp(float min, float x, float max) { // clamps x between min and max
    if (x < min)
        return min;
    else if ( x > max)
        return max;
    return x;
}

float absolute(float x) { // returns the absolute of x (|x|)
    if (x < 0)
        return -x;
    return x;
}

#endif