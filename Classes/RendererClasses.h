#ifndef RENDERER_H
#define RENDERER_H

#pragma once

#include <iostream>
#include <math.h>
#include <tools/rendererTools.h>

typedef unsigned char byte;

// Color Class

class Color
{
public:
    Color(byte _r = 255, byte _g = 255, byte _b = 255)
    {
        r = _r;
        g = _g;
        b = _b;
    }
    ~Color() {}

    byte r;
    byte g;
    byte b;

    void print()
    {
        std::cout << "(r:" << (int)r << ", g:" << (int)g << ", b:" << (int)b << ")" << std::endl;
    }

private:
};

// Ray Class

class Ray
{
public:
    Ray(Vector3 _op, Vector3 _a)
    {
        op = _op;
        a = _a;
    }
    ~Ray() {}

    Vector3 op;
    Vector3 a;

    Vector3 calc(float t) {
        return Vector3(op.x + t * a.x, op.y + t * a.y, op.z + t * a.z);
    }

private:
};

// intersection Data class
class intersectionData
{
public:
    intersectionData(bool _valid = false, float _t = 0, Vector3 intersectionPoint = Vector3()) {
        valid = _valid;
        t = _t;
        point = intersectionPoint;
    }
    ~intersectionData() {}

    bool valid = false;
    float t;
    Vector3 point;
};

// Polygon Class

class Polygon
{
public:
    Polygon(Vector3 _op, Vector3 _a, Vector3 _b, Color _color)
    {
        op = _op;
        a = _a;
        b = _b;
        color = _color;
    }
    ~Polygon() {}

    Vector3 op;
    Vector3 a;
    Vector3 b;
    Color color;

    Vector3 calc(float t, float s) {
        return Vector3(op.x + t * a.x + s * b.x, op.y + t * a.y + s * b.y, op.z + t * a.z + s * b.z);
    }
    intersectionData intersect(Ray g)
    {
        // The intersection works by intersecting a infinit Plane and checking if it was inside the Boundries
        
        // We do this by setting up a system of equations from the x, y and z of the Parameter Form of the line
        // and the parameter Form of the the Plane and setting them eqal and use the determinant method to solve it
        // g: X = OP1 + t1 * a1
        // e: X = OP2 + t2 * a2 + s * b
        // X = X => OP1 + t1 * a1 = OP2 + t2 * a2 + s * b
        
        // I: OP1x + t1 * a1x = OP2x + t2 * a2x + s * bx
        // II: OP1y + t1 * a1y = OP2y + t2 * a2y + s * by
        // III: OP1z + t1 * a1z = OP2z + t2 * a2z + s * bz

        float dk = det3x3({ // Determinant of all of the coefficients
            {a.x, b.x, -g.a.x}, 
            {a.y, b.y, -g.a.y}, 
            {a.z, b.z, -g.a.z}
        });

        if (dk == 0) // checks if it intersected
            return intersectionData(); // returns a invalid intersection

        float dt2 = det3x3({ // Determinant of a 3x3 matrix where the coefficients are swapped with the answers for the wanted Variable
            {g.op.x - op.x, b.x, -g.a.x}, // t2 is the t Parameter of the Polygon
            {g.op.y - op.y, b.y, -g.a.y}, 
            {g.op.z - op.z, b.z, -g.a.z}
        });
        float dt1 = det3x3({ // t1 is the t parameter of the Ray
            {a.x, b.x, g.op.x - op.x},
            {a.y, b.y, g.op.y - op.y},
            {a.z, b.z, g.op.z - op.z}
        });
        float ds = det3x3({ // s is the s parameter of the Polygon
            {a.x, g.op.x - op.x, -g.a.x},
            {a.y, g.op.y - op.y, -g.a.y},
            {a.z, g.op.z - op.z, -g.a.z}
        });

        float t1 = dt1 / dk; // calculation of the parameter of the intersection for the Ray
        if (t1 < 0)
            return intersectionData(); // returns a invalid solution if the Ray went backwards
        
        float t2 = dt2 / dk;  // calculation of the first parameter of the intersection for the Polygon
        float s = ds / dk; // calculation of the second parameter of the intersection for the Polygon
        if (t2 + s > 1 || t1 < 0 || s < 0) // checks if the intersection Point is inside the boundaries of the Polygon
            return intersectionData(); // it can happen that it is outside, because we actually intersected it with a infinte Plane not a Polygon
        
        return intersectionData(true, t1, g.calc(t1));
    }

private:
};

// Renderer Class

class Renderer
{
public:
    Renderer();
    ~Renderer();

private:

};

#endif