#include <iostream>
#include <fstream>
#include <string>
#include "include/renderer.h"

int main()
{
    std::ofstream file("intersectionTest.txt");

    std::string expr = "s = [";
    for (float i = 0; i < 1; i+= 0.01)
        expr += std::to_string(i) + ",";
    expr += "1]";
    file << expr << std::endl;

    std::vector<Polygon> polies = cube(Vector3(), Vector3(1, 1, 1), Color(255, 255, 255), 1, 1, 1);
    std::vector<Light> lights = {Light(Vector3(0, 1, 3), Color(), 1000)};
    Skybox sky(Vector3(0, 0, 0), Color(0, 191, 255), 1000, 1000);
    Renderer renderer(polies, lights, sky, Color(30, 30, 30), Vector3(0, -5, 0), Vector3(0, 0.01f, 0), 200, 150, 1, 3, 2, 0.2f);
    Ray r = Ray(renderer.c, Vector3(0, 1, 0));
    intersectionData data = renderer.rayCast(r);

    // output for testing
    std::cout << "Ray: ";
    r.print();
    file << "X = (" << r.op.x << "," << r.op.y << "," << r.op.z << ") + (" << r.a.x << "," << r.a.y << "," << r.a.z << ") * t" << std::endl;
    for (int i = 0; i < polies.size(); i++)
        file << "X_{" << i << "} = ("<<polies[i].op.x<<","<<polies[i].op.y<<","<<polies[i].op.z<<") + ("<<polies[i].a.x<<","<<polies[i].a.y<<","<<polies[i].a.z<<") * t + ("<<polies[i].b.x<<","<<polies[i].b.y<<","<<polies[i].b.z<<") * s {t+s<=1}{t>=0}{s>=0}" << std::endl;
    file << sky.radius<<"^2 = (x - "<<sky.position.x<<")^2 + (y - "<<sky.position.y<<")^2 + (z - "<<sky.position.z<<")^2" << std::endl;
    if (data.valid)
        file << "S = ("<<data.point.x<<","<<data.point.y<<","<<data.point.z<<")" << std::endl;
    if (data.valid) {
        if (data.shapeID == POLYGON_ID) {
            std::cout << ", Polygon: ";
            renderer.polygons[data.i].print();
            std::cout << ", intersection Point: ";
            data.point.print();
            std::cout << std::endl;
            Ray ref = renderer.reflect(r, data);
            file << "X_r = ("<<ref.op.x<<","<<ref.op.y<<","<<ref.op.z<<") + ("<<ref.a.x<<","<<ref.a.y<<","<<ref.a.z<<") * t" << std::endl;;
        }
        else
        {
            std::cout << ", Sphere: ";
            renderer.spheres[data.i].print();
            std::cout << ", intersection Point: ";
            data.point.print();
            std::cout << std::endl;
        }
    }
    else
    {
        std::cout << " intersected with nothing" << std::endl;
    }
    file.close();

    writePNG("result.png", renderer.render());
    std::cout << "render completed!" << std::endl;
    /*
    writePNG("C:/Users/Thomas Shuttleworth/Desktop/RecursiveRendererCpp/result.png", renderer.result);
    std::cout << "finished!" << std::endl;
    */
}