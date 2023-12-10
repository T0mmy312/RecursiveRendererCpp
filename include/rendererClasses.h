#pragma once

#ifndef _RENDERERCLASSES_H_
#define _RENDERERCLASSES_H_

#include <vector>
#include <iostream>
#include <math.h>
#include <cstdlib> 
#include <ctime>
#include "mathAndVectorTools.h"

typedef unsigned char byte;

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Color Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Color
{
public:
    Color(byte _r = 255, byte _g = 255, byte _b = 255){
        r = clamp(0, _r, 255);
        g = clamp(0, _g, 255);
        b = clamp(0, _b, 255);
    }
    ~Color() {}

    byte r; // red (max 255)
    byte g; // green (max 255)
    byte b; // blue (max 255)

    void print() // writes the color out result: (r: red, g: green, b: blue)
    {
        std::cout << "(r:" << (int)r << ", g:" << (int)g << ", b:" << (int)b << ")" << std::endl;
    }

    Color operator+(Color const& obj)
    {
        return Color((byte)clamp(0, (int)(r) + (int)(obj.r), 255), (byte)clamp(0, (int)(g) + (int)(obj.g), 255), (byte)clamp(0, (int)(b) + (int)(obj.b), 255));
    }
    Color operator-(Color const& obj)
    {
        return Color((byte)clamp(0, (int)(r) - (int)(obj.r), 255), (byte)clamp(0, (int)(g) - (int)(obj.g), 255), (byte)clamp(0, (int)(b) - (int)(obj.b), 255));
    }
    Color operator/(Color const& obj)
    {
        return Color(r / obj.r, g / obj.g, b / obj.b);
    }
    Color operator*(float const& obj)
    {
        return Color((byte)clamp(0, (int)(r) * obj, 255), (byte)clamp(0, (int)(g) * obj, 255), (byte)clamp(0, (int)(b) * obj, 255));
    }

private:
};

typedef std::vector<std::vector<Color>> Picture; // a 2d vector of color Values

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Ray Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Ray
{
public:
    Ray(Vector3 _op = Vector3(), Vector3 _a = Vector3()){
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
    intersectionData(bool _valid = false, float _t = 0, Vector3 intersectionPoint = Vector3(), int _i = -1) {
        valid = _valid;
        t = _t;
        point = intersectionPoint;
        i = _i;
    }
    ~intersectionData() {}

    bool valid = false; // says if the intersection is valid or if there is no intersection of the two shapes
    float t; // the t parameter of the Ray
    Vector3 point; // the intersection point
    int i; // index used for raycasting
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

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Light Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Light
{
public:
    Light(Vector3 _position = Vector3(), Color _color = Color(), float _intensity = 1, float _radius = 0.1f, bool _valid = true) {
        position = _position;
        color = _color;
        intensity = _intensity;
        radius = _radius;
        valid = _valid;
    }
    ~Light() {}

    Vector3 position;
    Color color;
    float intensity;
    float radius;
    bool valid;

private:
};

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Renderer Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Renderer
{
public:
    intersectionData rayCast(Ray g)
    {
        intersectionData nearest = intersectionData();
        for (int i = 0; i < polygons.size(); i++)
        {
            intersectionData current = polygons[i].intersect(g);
            current.i = i;
            if (nearest.valid && current.t < nearest.t)
                nearest = current;
            else if (!nearest.valid && current.valid)
                nearest = current;
        }
        return nearest;
    }

private:
    int LightIntersect(Ray g) // returns a valid index in lights if the ray has hit a light else -1
    {
        for (int i = 0; i < lights.size(); i++)
        {
            Vector3 opp = lights[i].position - g.op;
            if (opp.magnitude() * tan(angleInRads(g.a, opp)) < lights[i].radius)
                return i;
        }
        return -1; // not intersecting with a light
    }

    Ray reflect(Ray g, intersectionData data) // calculates the perfect reflected ray of a surface
    {
        Ray ret(data.point);
        if (!data.valid)
            return ret;
    
        Vector3 nn = polygons[data.i].normalVector().normalized();
        if (angleInRads(nn, g.a * -1) > 1.57079632679)
            nn = nn * -1;
        Vector3 an = g.a.normalized();
        Vector3 nt = nn * scalarProd(an, nn);
        ret.a = nt * 2 - an;
        return ret;
    }

    Color recursvieRay(Ray g, int bounces) //! Check how to handle Colors, because it is pobably totally wrong
    {
        if (bounces > maxReflections) // checks if it has exceeded the max bounces
            return defaultBackgroundColor; // returns the default background Color
        
        int cr = 0;
        int cg = 0;
        int cb = 0;

        int light = LightIntersect(g);
        if (light != -1) // if intersected with a light returns the color
            return lights[light].color * (1/pow((lights[light].position - g.op).magnitude(), 2)) * lights[light].intensity;
        
        intersectionData mainRay = rayCast(g);
        if (!mainRay.valid)
            return defaultBackgroundColor;
        
        Ray newRay = reflect(g, mainRay);
        Color val = recursvieRay(newRay, bounces + 1)/*Add this at the end before return * (1/pow((mainRay.point - g.op).magnitude(), 2))*/;
        cr += val.r;
        cg += val.g;
        cb += val.b;

        for (int i = 0; i < splits - 1; i++)
        {
            float xRand = (((rand() % 1000) / 100) - 5) * polygons[mainRay.i].scatter;
            float yRand = (((rand() % 1000) / 100) - 5) * polygons[mainRay.i].scatter;
            float zRand = (((rand() % 1000) / 100) - 5) * polygons[mainRay.i].scatter;
            val = recursvieRay(Ray(newRay.op, Vector3(newRay.a.x + xRand, newRay.a.y + yRand, newRay.a.z + zRand)), bounces + 1) * (1/(1+xRand+yRand+zRand));
            cr += val.r;
            cg += val.g;
            cb += val.b;
        }

        return Color(cr/splits, cg/splits, cb/splits) * (1/pow((mainRay.point - g.op).magnitude(), 2)) - (Color(255, 255, 255) - polygons[mainRay.i].color);
    }

public:
    Renderer(std::vector<Polygon> _polygons, std::vector<Light> _lights, Color _defaultBackgroundColor, Vector3 _c, Vector3 _f, int _xPixls, int _yPixls, float _refraction, int _maxReflections, int _splits, float _xSize) { //c is the camera Pos and f is the focal lenght and rotation of the camera
        polygons = _polygons;
        lights = _lights;
        defaultBackgroundColor = _defaultBackgroundColor;
        c = _c;
        f = _f;
        xPixls = _xPixls;
        yPixls = _yPixls;
        splits = _splits;
        refraction = _refraction;
        maxReflections = _maxReflections;
        xSize = _xSize;
        result = Picture(_yPixls, std::vector<Color>(_xPixls, Color(0, 0, 0)));

        // set random seed
        srand(time(0));
    }
    ~Renderer() {}

    std::vector<Polygon> polygons; // all Polygons in the scene
    std::vector<Light> lights; // all Lights in the scene

    Color defaultBackgroundColor;

    Vector3 c; // current position of the Camera
    Vector3 f; // focal lenght and the direction the Camera points

    int xPixls; // render x amount of pixels
    int yPixls; // render y amount of pixels

    int splits; // how many different Rays result from a reflection of a Ray
    int refraction; // refraction of the camera Lense (0 - 1: 1 is no lense and 0 is all paralell)

    int maxReflections; // maximum of bounces until it reaches a light source

    float xSize; // x size of the screen in world space

    Picture result; // the resulting picture of the render

    Picture render() // renders the Polygons and stores the pixl values in the result var and returns that
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
            fx = crossProd(f, Vector3(0, 0, 1)).normalized();
        Vector3 fy = crossProd(f, fx).normalized(); // the x Vector of the screens global position and rotation

        Vector3 fn = f.normalized(); // normalized vector f
        
        for (int y = 0; y < yPixls; y++)
        {
            for (int x = 0; x < xPixls; x++)
            {
                float lx = (x - hsx + addX) * pixToM; // lenght of fx needed to get to point
                float ly = (y - hsy + addY) * pixToM; // lenght of fy needed to get to point

                Vector3 op = c + fx * lx + fy * ly; // finds the global position of the x, y Pixel

                Vector3 copn = (c - op).normalized(); // normalized vector going from the camera position to the Pixel Position
                Vector3 a = fn + (fn - copn * scalarProd(copn, fn)) * refraction; // direction Vector of the initial Ray according to refraction
                // sidenote: scalarProd(copn, fn) is the cos of the angle between the two vectors, because |copn| = |fn| = 1 (makes it slightly more efficient)

                Ray g(op, a);
                result[y][x] = recursvieRay(g, 0);
            }
        }

        return result;
    }

    void appendPolygons(std::vector<Polygon> polies) // simply appends a list of polygons to polygons (like from a cube funktion)
    {
        for (int i = 0; i < polies.size(); i++)
            polygons.push_back(polies[i]);
    }
};

#endif