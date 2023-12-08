#pragma once

#ifndef _RENDERERTOOLS_H_
#define _RENDERERTOOLS_H_

#include "rendererClasses.h"

//? ----------------------------------------------------------------------------------------------------------------------------------
//? Tools for easier Scene building
//? ----------------------------------------------------------------------------------------------------------------------------------

std::vector<Polygon> cube(Vector3 pos, Vector3 dir, Color color, float sx, float sy, float sz) // adds needed Polygons to get a cube a the position pos with side lenghts of sx, sy, sz and directin dir in Renderer renderer (returns from where to wher the Polygons are)
{
    std::vector<Polygon> ret;
    
    float hsx = sx / 2;
    float hsy = sy / 2;
    float hsz = sz / 2;

    Vector3 rx;
    Vector3 ry;
    Vector3 rz;
    rx = dir.normalized();
    if (rx.x == 0 && rx.y == 0 && rx.z == 1)
    {
        ry = Vector3(1, 0, 0);
        rz = Vector3(0, -1, 0);
    }
    else
    {
        ry = crossProd(rx, Vector3(0, 0, 1)).normalized();
        rz = crossProd(rx, ry).normalized();
    }

    Vector3 bfr = rx * hsx + ry * hsy - rz * hsz; // bottom front right corner of the cube
    Vector3 bbr = rx * -hsx + ry * hsy - rz * hsz; // bottom back right corner of the cube
    Vector3 bfl = rx * hsx - ry * hsy - rz * hsz; // bottom front left corner of the cube
    Vector3 bbl = rx * -hsx - ry * hsy - rz * hsz; // bottom back left corner of the cube
    Vector3 tfr = rx * hsx + ry * hsy + rz * hsz; // top front right corner of the cube
    Vector3 tbr = rx * -hsx + ry * hsy + rz * hsz; // top back right corner of the cube
    Vector3 tfl = rx * hsx - ry * hsy + rz * hsz; // top front left corner of the cube
    Vector3 tbl = rx * -hsx - ry * hsy + rz * hsz; // top back left corner of the cube

    ret.push_back(Polygon(bfl, bfr - bfl, bbl - bfl, color));
    ret.push_back(Polygon(bbr, bfr - bbr, bbl - bbr, color));
    ret.push_back(Polygon(bbl, tbl - bbl, bbr - bbl, color));
    ret.push_back(Polygon(tbr, tbl - tbr, bbr - tbr, color));
    ret.push_back(Polygon(tbr, bbr - tbr, tfr - tbr, color));
    ret.push_back(Polygon(bfr, bbr - bfr, tfr - bfr, color));
    ret.push_back(Polygon(tfr, tfl - tfr, bfr - tfr, color));
    ret.push_back(Polygon(bfl, tfl - bfl, bfr - bfl, color));
    ret.push_back(Polygon(bfl, tfl - bfl, bbl - bfl, color));
    ret.push_back(Polygon(tbl, tfl - tbl, bbl - tbl, color));
    ret.push_back(Polygon(tfl, tbl - tfl, tfr - tfl, color));
    ret.push_back(Polygon(tbr, tfr - tbr, tbl - tbr, color));

    return ret;
}

#endif