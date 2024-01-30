#include <iostream>
#include <fstream>
#include "include/renderer.h"

int main() 
{
    std::fstream file("quaternionTest.txt");
    Vector3 axis = Vector3(1, 1, 1).normalized();
    Vector3 point(1, 2, 3);
    file << "\\operatorname{vector}((0, 0, 0), ("<<axis.x*5<<","<<axis.y*5<<","<<axis.z*5<<"))" << std::endl;
    file << "p = ("<<point.x<<","<<point.y<<","<<point.z<<")" << std::endl;
    Quaternion rot;
    rot.setAxisRotationDeg(axis, 180);
    Vector3 pr = rotateActive(point, rot);
    file << "p_r = ("<<pr.x<<","<<pr.y<<","<<pr.z<<")" << std::endl;
    file.close();
}