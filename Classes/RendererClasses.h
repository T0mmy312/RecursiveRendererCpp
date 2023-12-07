#ifndef RENDERER_H
#define RENDERER_H

#pragma once

#include <iostream>
#include <math.h>
#include <vector>
#include "tools/rendererTools.h"

typedef unsigned char byte;

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Color Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Color
{
public:
    Color(byte _r = 255, byte _g = 255, byte _b = 255){
        r = _r;
        g = _g;
        b = _b;
    }
    ~Color() {}

    byte r; // red (max 255)
    byte g; // green (max 255)
    byte b; // blue (max 255)

    void print() // writes the color out result: (r: red, g: green, b: blue)
    {
        std::cout << "(r:" << (int)r << ", g:" << (int)g << ", b:" << (int)b << ")" << std::endl;
    }

private:
};

typedef std::vector<std::vector<Color>> picture; // a 2d vector of color Values

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Ray Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Ray
{
public:
    Ray(Vector3 _op, Vector3 _a){
        op = _op;
        a = _a;
    }
    ~Ray() {}

    Vector3 op; // origin Point of the Ray
    Vector3 a; // direction Vector of the Ray

    Vector3 calc(float t) { // returns a Point on the Ray for the parameter t
        return Vector3(op.x + t * a.x, op.y + t * a.y, op.z + t * a.z);
    }

private:
};

//? ----------------------------------------------------------------------------------------------------------------------------------
//? intersection data Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class intersectionData
{
public:
    intersectionData(bool _valid = false, float _t = 0, Vector3 intersectionPoint = Vector3()) {
        valid = _valid;
        t = _t;
        point = intersectionPoint;
    }
    ~intersectionData() {}

    bool valid = false; // says if the intersection is valid or if there is no intersection of the two shapes
    float t; // the t parameter of the Ray
    Vector3 point; // the intersection point
};

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Polygon Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Polygon
{
public:
    Polygon(Vector3 _op, Vector3 _a, Vector3 _b, Color _color){
        op = _op;
        a = _a;
        b = _b;
        color = _color;
    }
    ~Polygon() {}

    Vector3 op; // the origin Point of the Polygon
    Vector3 a; // a side of the Polygon (going out from the op Vector)
    Vector3 b; // another side of the Polygon (going out from the op Vector)
    Color color; // the Color of the Polygon
    float scatter; // 0 - 1: max 0° of random scatter (so a perfect mirror) - max 90° of random scatter (aka mattness)

    Vector3 calc(float t, float s) { // returns a point on the Plane going through the Polygon (doesn't have to be on the Polygon)
        return Vector3(op.x + t * a.x + s * b.x, op.y + t * a.y + s * b.y, op.z + t * a.z + s * b.z);
    }
    Vector3 normalVector() { // returns the normal Vector to this polygon
        return scalarProd(a, b);
    }
    intersectionData intersect(Ray g) // intersects a Polygon with a Ray
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

//! Add Light Class

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Renderer Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Renderer
{
public:
    Renderer(std::vector<Polygon> _polygons, Vector3 _c, Vector3 _f, int _xPixls, int _yPixls, float _refraction, int _maxReflections, int _splits, float _xSize) { //c is the camera Pos and f is the focal lenght and rotation of the camera
        polygons = _polygons;
        c = _c;
        f = _f;
        xPixls = _xPixls;
        yPixls = _yPixls;
        splits = _splits;
        refraction = _refraction;
        maxReflections = _maxReflections;
        xSize = _xSize;
        result = std::vector(_yPixls, std::vector(_xPixls, Color(0, 0, 0)));
    }
    ~Renderer() {}

    std::vector<Polygon> polygons; // all Polygons in the scene
    //! Add lights vector

    Vector3 c; // current position of the Camera
    Vector3 f; // focal lenght and the direction the Camera points

    int xPixls; // render x amount of pixels
    int yPixls; // render y amount of pixels

    int splits; // how many different Rays result from a reflection of a Ray
    int refraction; // refraction of the camera Lense (0 - 1: 1 is no lense and 0 is all paralell)

    int maxReflections; // maximum of bounces until it reaches a light source

    float xSize; // x size of the screen in world space

    picture result; // the resulting picture of the render

    picture render() // renders the Polygons and stores the pixl values in the result var and returns that
    {
        // The Renderer works by casting Rays and reflecting them until they reach a Light

        float pixToM = xSize / xPixls; // variable to convert pixl mesurements in global units (m)

        float addX = 1; // add to x to compensate for it being odd or even
        float addY = 1; // add to y to compensate for it being odd or even
        if (xPixls % 2 == 0)
            addX = 0.5f;
        if (yPixls % 2 == 0)
            addY = 0.5f;
        
        float hsx = xPixls / 2; // half of the x Pixls
        float hsy = yPixls / 2; // half of the y Pixls

        Vector3 fx; // the x Vector of the screens global position and rotation
        if (f.x == 0 && f.y == 1 && f.y == 0)
            fx = Vector3(1, 0, 0);
        else
            fx = scalarProd(f, Vector3(0, 0, 1)).normalized();
        Vector3 fy = scalarProd(f, fx).normalized(); // the x Vector of the screens global position and rotation
        
        for (int y = 0; y < yPixls; y++)
        {
            for (int x = 0; x < xPixls; x++)
            {
                float lx = (x - hsx + addX) * pixToM; // lenght of fx needed to get to point
                float ly = (y - hsy + addY) * pixToM; // lenght of fy needed to get to point

                Vector3 op = c + fx * lx + fy * ly; // finds the global position of the x, y Pixel
            }
        }

        return result;
    }

private:
};

#endif