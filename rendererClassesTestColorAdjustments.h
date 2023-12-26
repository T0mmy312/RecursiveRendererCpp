#pragma once

#ifndef _RENDERERCLASSES_H_
#define _RENDERERCLASSES_H_

#include <vector>
#include <iostream>
#include <math.h>
#include <cstdlib> 
#include <ctime>
#include "include/mathAndVectorTools.h"

#define POLYGON_ID 0
#define SPHERE_ID 1

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
    intersectionData(bool _valid = false, int _shapeID = -1, float _t = 0, Vector3 intersectionPoint = Vector3(), int _i = -1) {
        valid = _valid;
        shapeID = _shapeID;
        t = _t;
        point = intersectionPoint;
        i = _i;
    }
    ~intersectionData() {}

    bool valid = false; // says if the intersection is valid or if there is no intersection of the two shapes
    int shapeID = -1;
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

    
        float dk = det3x3({ // Determinant of a 3x3 matrix where the coefficients are swapped with the answers for the wanted Variable
            {a.x, b.x, -g.a.x}, // Determinant of all of the coefficients
            {a.y, b.y, -g.a.y}, 
            {a.z, b.z, -g.a.z}
        });

        if (dk == 0) // checks if it intersected
            return intersectionData(); // returns a invalid intersection

        float dt1 = det3x3({ // t1 is the t parameter of the Ray
            {a.x, b.x, g.op.x - op.x},
            {a.y, b.y, g.op.y - op.y},
            {a.z, b.z, g.op.z - op.z}
        });

        float t1 = dt1 / dk; // calculation of the parameter of the intersection for the Ray
        if (t1 < 0)
            return intersectionData(); // returns a invalid solution if the Ray went backwards

        float dt2 = det3x3({ // t2 is the t Parameter of the Polygon
            {g.op.x - op.x, b.x, -g.a.x},
            {g.op.y - op.y, b.y, -g.a.y}, 
            {g.op.z - op.z, b.z, -g.a.z}
        });
        float ds = det3x3({ // s is the s parameter of the Polygon
            {a.x, g.op.x - op.x, -g.a.x},
            {a.y, g.op.y - op.y, -g.a.y},
            {a.z, g.op.z - op.z, -g.a.z}
        });
        
        float t2 = dt2 / dk;  // calculation of the first parameter of the intersection for the Polygon
        float s = ds / dk; // calculation of the second parameter of the intersection for the Polygon
        if (t2 + s > 1 || t2 < 0 || s < 0) // checks if the intersection Point is inside the boundaries of the Polygon
            return intersectionData(); // it can happen that it is outside, because we actually intersected it with a infinte Plane not a Polygon
        
        return intersectionData(true, POLYGON_ID, t1, g.calc(t1));
    }

private:
};

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Sphere Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Sphere
{
public:
    Sphere(Vector3 _op = Vector3(), float _sx = 1, float _sy = 1, float _sz = 1, float _radius = 1) {
        op = _op;
        sx = _sx;
        sy = _sy;
        sz = _sz;
        radius = _radius;
    }
    ~Sphere() {}

    Vector3 op; // origin Point of the Sphere
    float sx; // scale x
    float sy; // scale y
    float sz; // scale z
    float radius; // radius of the Sphere
    // general formular for a sphere: radius^2 = (x - op.x)^2 / sx + (y - op.y)^2 / sy + (z - op.z)^2 / sz

    Vector3 normalVector(Vector3 intersectionPoint) { // calculates the normal Vector for a spezific intersection Point on the surface of the sphere
        return Vector3((intersectionPoint.x - op.x)/sx, (intersectionPoint.y - op.y)/sy, (intersectionPoint.y - op.y)/sy);
    }
    intersectionData intersect(Ray g) // intersects a ray with the sphere returns the nearest valid not negativ (on the line) intersection Point
    {
        // a, b and c are the variables for the long quadratic formular they come from:
        // g: X = op + t * a (fromular for the line)
        // s: radius^2 = (x - op.x)^2 / sx + (y - op.y)^2 / sy + (z - op.z)^2 / sz (formular for the sphere)
        // then we just need to insert the x, y and z formulars of the Line into the sphere formular and turn it into the 0 = ax^2 + bx + c form
        // -> radius^2 = (g.op.x + t * g.a.x - op.x)^2 / sx + (g.op.y + t * g.a.y - op.y)^2 / sy + (g.op.y + t * g.a.y - op.z)^2 / sz
        // -> 0 = (g.a.x^2/sx + g.a.y^2/sy + g.a.z^2/sz) * t^2 + ((2 * (g.op.x - op.x) * a.x)/sx + (2 * (g.op.y - op.y) * a.y)/sy + (2 * (g.op.z - op.z) * a.z)/sz) * t + ((g.op.x - op.x)^2/sx + (g.op.y - op.y)^2/sy + (g.op.z - op.z)^2/sz + radius^2)
        // -> a = g.a.x^2/sx + g.a.y^2/sy + g.a.z^2/sz
        //    b = (2 * (g.op.x - op.x) * a.x)/sx + (2 * (g.op.y - op.y) * a.y)/sy + (2 * (g.op.z - op.z) * a.z)/sz
        //    c = (g.op.x - op.x)^2/sx + (g.op.y - op.y)^2/sy + (g.op.z - op.z)^2/sz + radius^2
        // -> t1,t2 = (-b (+,-) sqrt(b^2 - 4*a*c))/(2*a)

        double a = (g.a.x * g.a.x) / sx + (g.a.y * g.a.y) / sy + (g.a.z * g.a.z) / sz;
        double b = (2 * (g.op.x - op.x) * g.a.x)/sx + (2 * (g.op.y - op.y) * g.a.y)/sy + (2 * (g.op.z - op.z) * g.a.z)/sz;
        double c = sqr(g.op.x - op.x)/sx + sqr(g.op.y - op.y)/sy + sqr(g.op.z - op.z)/sz + radius * radius;

        if (2*a == 0) // checks if the division is invalid if so return a invalid intersection
            return intersectionData();
        double rootContent = b*b - 4 * a * c; // calculates the content of the root
        if (rootContent < 0) // checks if it is an invalid root if so returns an invalid intersection
            return intersectionData();
        if (rootContent == 0) // checks if the content of the root is 0 if so the formular is simpler and there is only one intersection Point
        {
            float t = (-b)/(2 * a); // calcultes t with (-b (+,-) sqrt(b^2 - 4*a*c))/(2*a) but with sqrt(b^2 - 4*a*c) being 0 so t = (-b)/(2*a)
            if (t < 0) // if t is negativ it returns an invalid intersection, because the only intersection is if the ray where to travel backwards
                return intersectionData();
            return intersectionData(true, SPHERE_ID, t, g.calc(t)); // returns the valid intersection
        }
        double root = sqrt(rootContent); // calculates the root once, because it could be that it is used twice
        float t = (-b - root)/(2 * a); // calculates the nearest t to the origin
        if (t < 0) // if that t is negativ
            t = (-b + root)/(2 * a); // calcultes the other posibillity for t
            if (t < 0) // if that is also negativ it returns a invalid intersection
                return intersectionData();
            return intersectionData(true, SPHERE_ID, t, g.calc(t)); // returns the other valid intersection
        return intersectionData(true, SPHERE_ID, t, g.calc(t)); // returns the nearer valid intersection
    }

private:
};

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Light Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Light
{
public:
    Light(Vector3 _position = Vector3(), Color _color = Color(), float _intensity = 1, float _radius = 0.5f, bool _valid = true) {
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
//? Color and intensity data Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class scai // stands for Single Color And Intensity
{
public:
    scai(Color _color = Color(), double _intensity = 1) {
        color = _color;
        intensity = _intensity;
    }
    ~scai() {}

    Color color;
    double intensity;

public:
};

typedef std::vector<std::vector<scai>> iMatrix; // intensity matrix typedef

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Recursive Return Data Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class RecursiveReturnData // return data of the recursvieRay funktion
{
public:
    RecursiveReturnData(scai _lightVals, long double _distance = -1) {
        lightVals = _lightVals;
        distance = _distance;
    }
    ~RecursiveReturnData() {}

    scai lightVals;
    long double distance;

private:
};

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Renderer Class
//? ----------------------------------------------------------------------------------------------------------------------------------

class Renderer
{
public:
    intersectionData rayCast(Ray g, int except = -1, int except_shape_id = -1) // returns the nearest intersection of any object in polygons and spheres (execpt for the execption)
    {
        intersectionData nearest = intersectionData(); // sets a default invalid intersection point
        for (int i = 0; i < polygons.size(); i++) // loops through every polygon in polygons
        {
            if (i == except && except_shape_id == POLYGON_ID) // checks if this is the exception if so it goes to the next loop
                continue;
            intersectionData current = polygons[i].intersect(g); // calculates the intersection of the current Polygon
            current.i = i; // sets the i compomponent of the intersection data to the current index
            if ((!nearest.valid || current.t < nearest.t) && current.valid) // checks if the intersection is valid and nearer then the current min
                nearest = current; // replaces the nearest intersection with current intersection
        }
        for (int i = 0; i < spheres.size(); i++) // loops through all spheres
        {
            if (i == except && except_shape_id == SPHERE_ID) // if this is the execption it continues
                continue;
            intersectionData current = spheres[i].intersect(g); // calcultes the current intersection
            current.i = i; // sets the i value to the current index
            if ((!nearest.valid || current.t < nearest.t) && current.valid) // checks if the current intersection is nearer then nearest and is valid
                nearest = current; // sets nearest to current
        }
        return nearest; // returns the nearest intersection
    }

private:
    int LightIntersect(Ray g) // returns a valid index in lights if the ray has hit a light else -1
    {
        for (int i = 0; i < lights.size(); i++) // loops through all Lights
        {
            Vector3 opp = lights[i].position - g.op; // calculates the vector between the origin of the ray and the light
            if (opp.magnitude() * tan(angleInRads(g.a, opp)) < lights[i].radius) // if the lenght of the opposite side is shorter then the radius of the Light
                return i; // returns the index of the Light
        }
        return -1; // not intersecting with a light
    }

    Ray reflect(Ray g, intersectionData data) // calculates the perfect reflected ray of a surface
    {
        Ray ret(data.point); // sets a default invalid returns ray where the directional vector a of ret is (0, 0, 0)
        if (!data.valid || data.i == -1) // if the intersection data is invalid it returns the invalid ray
            return ret;

        Vector3 nn; // declares the normalized normal vector of the reflecting surface
        if (data.shapeID == POLYGON_ID) // checks if that surface is a polygon
            nn = polygons[data.i].normalVector().normalized(); // calculates the normalized normal vector of that polygon
        else if (data.shapeID == SPHERE_ID) // checks if it is a sphere
            nn = spheres[data.i].normalVector(data.point).normalized(); // calculates the normalized normal vector of that intersection point on that sphere
        float alpha = angleInRads(nn, g.a * -1); // calculates the angle between the normalized normal vector and the reversed a vector
        if (alpha > 1.57079632679) // checks if the angle between the normal vector and the reversed a vector is more then 90° (in radians) 
            nn = nn * -1; // reverses the normalized normal vector
        Vector3 an = g.a * cos(alpha) * -1; // calculates the reverse and lenght adjusted vector an
        Vector3 nt = nn - an; // calculates the vector between the normal vector and the an vector
        ret.a = nn + nt; // calculates the outgoing a vector of the reflection
        return ret; // returns ret
    }

    RecursiveReturnData recursvieRay(Ray g, int bounces, int except = -1, int except_shape_id = -1) // fires a light ray that then bounces and sents splits more rays out that do the same until it reaces a light
    {
        //! Check how to handle Colors, because it is pobably totally wrong
        if (bounces > maxReflections) // checks if it has exceeded the max bounces
            return RecursiveReturnData(scai(defaultBackgroundColor, defaultBackgroundIntensity), 0); // returns the default background Color
        
        int cr = 0; // stores the avariged red values of the next bounces
        int cg = 0; // stores the avariged green values of the next bounces
        int cb = 0; // stores the avariged blue values of the next bounces

        int light = LightIntersect(g); // checks if this ray intersects with a light
        intersectionData mainRay = rayCast(g, except, except_shape_id);
        if (!mainRay.valid) { // checks if the rayCast is valid
            if (light != -1) // if intersected with a light returns the color
                return RecursiveReturnData(scai(lights[light].color, lights[light].intensity), (lights[light].position - g.op).magnitude()); // returns the light data and the distance this ray traveled
            return scai(defaultBackgroundColor, defaultBackgroundIntensity); // returns the default background color and default background intensity
        }
        if ((lights[light].position - g.op).magnitude() < (mainRay.point - g.op).magnitude()) // checks if the light in nearer then the rayCast if so it returns the light
            return RecursiveReturnData(scai(lights[light].color, lights[light].intensity), (lights[light].position - g.op).magnitude()); // returns the light data and the distance this ray traveled

        Ray newRay = reflect(g, mainRay); // reflects the the rayCast Ray
        RecursiveReturnData newRayData = recursvieRay(newRay, bounces + 1, mainRay.i, mainRay.shapeID); // creates a new recursvieRay with the reflected Ray
        double intensity = newRayData.lightVals.intensity * 0.5f; // calculates the intensity with a weighting of 50% for this return value
        long double distance = newRayData.distance * 0.5f; // calculates the distance with a weighting of 50% for this return value
        cr += newRayData.lightVals.color.r * 0.5f; // calculates the red value with a weighting of 50% for this return value
        cg += newRayData.lightVals.color.g * 0.5f; // calculates the green value with a weighting of 50% for this return value
        cb += newRayData.lightVals.color.b * 0.5f; // calculates the blue value with a weighting of 50% for this return value

        int mulVal = 0.5f / (splits - 1); // calculates the weighting for the next rays

        for (int i = 0; i < splits - 1; i++) // generates so many other rays as splits - 1 (because one was alredy generated)
        {
            float xRand = (((rand() % 1000) / 100) - 5) * polygons[mainRay.i].scatter; // generates a random x offset for the next ray with the scatter variable of Polygon
            float yRand = (((rand() % 1000) / 100) - 5) * polygons[mainRay.i].scatter; // generates a random y offset for the next ray with the scatter variable of Polygon
            float zRand = (((rand() % 1000) / 100) - 5) * polygons[mainRay.i].scatter; // generates a random z offset for the next ray with the scatter variable of Polygon
            RecursiveReturnData rayData = recursvieRay(Ray(newRay.op, Vector3(newRay.a.x + xRand, newRay.a.y + yRand, newRay.a.z + zRand)), bounces + 1, mainRay.i, mainRay.shapeID); // generates another recursvieRay and gets the return data
            distance += rayData.distance * mulVal; // adds to the total distance with the desired weighting of this ray
            intensity += rayData.lightVals.intensity * mulVal; // adds to the total intensity with the desired weighting of this ray
            cr += rayData.lightVals.color.r * mulVal; // adds to the red value with the desired weighting of this ray
            cg += rayData.lightVals.color.g * mulVal; // adds to the green value with the desired weighting of this ray
            cb += rayData.lightVals.color.b * mulVal; // adds to the blue value with the desired weighting of this ray
        }

        return RecursiveReturnData(scai(Color(cr, cg, cb)/*- (Color(255, 255, 255) - polygons[mainRay.i].color)*/, intensity), distance + (mainRay.point - g.op).magnitude()); // returns the data collected from the other return values
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
    std::vector<Sphere> spheres; // all spheres in the scene
    std::vector<Light> lights; // all Lights in the scene

    Color defaultBackgroundColor;
    float defaultBackgroundIntensity;

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

        iMatrix preResult(yPixls, std::vector<scai>(xPixls, scai()));
        double minIntensity = 999999999;
        double maxIntensity = -999999999;
        
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
                RecursiveReturnData retData = recursvieRay(g, 0);
                preResult[y][x] = scai(retData.lightVals.color, 1 / pow(retData.distance + (1 / sqrt(retData.lightVals.intensity)), 2));
                // this "1 / pow(retData.distance + (1 / sqrt(retData.lightVals.intensity)), 2)" calculates the intensity for retData.distance where if that where 0 the intensity would be retData.lightVals.intensity
                if (preResult[y][x].intensity < minIntensity)
                    minIntensity = preResult[y][x].intensity;
                else if (preResult[y][x].intensity > maxIntensity)
                    maxIntensity = preResult[y][x].intensity;
            }
        }

        if (minIntensity == 999999999 || maxIntensity == -999999999) // returns an Error if the Intensity probably wasn't calculated correctly
        {
            std::cout << "Could not render (probably an Error with Light intensity)!" << std::endl;
            return Picture();
        }
        double dif = maxIntensity - minIntensity;
        for (int y = 0; y < yPixls; y++)
        {
            for (int x = 0; x < xPixls; x++)
            {
                result[y][x] = preResult[y][x].color * ((preResult[y][x].intensity - minIntensity) / dif);
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