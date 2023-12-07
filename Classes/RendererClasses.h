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

    Vector3 calc(float t)
    {
        Vector3 ret(op.x + t * a.x, op.y + t * a.y, op.z + t * a.z);
        return ret;
    }

private:
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

    Vector3 calc(float t, float s)
    {
        Vector3 ret(op.x + t * a.x + s * b.x, op.y + t * a.y + s * b.y, op.z + t * a.z + s * b.z);
        return ret;
    }
    Vector3 intersect(Ray g)
    {
        
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