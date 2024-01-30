#include <iostream>
#include <fstream>
#include <string>
#include "include/renderer.h"

int main()
{
    std::vector<Polygon> polies = cube(Vector3(), Vector3(1, 1, 1), Color(255, 255, 255), 2, 2, 2);
    std::vector<Light> lights = {Light(Vector3(0, -1, 6), Color(0, 150, 0), 5, 2)};
    Skybox sky(Vector3(0, 0, 0), Color(0, 191, 255), 100, 10);
    Renderer renderer(polies, lights, sky, Color(30, 30, 30), Vector3(0, -5, 0), Vector3(0, 0.1, 0), 400, 300, 1, 3, 2, 1);

    // output for testing
    std::ofstream file("intersectionTest.txt");

    Vector3 cf = renderer.c + renderer.f;
    file << "\\operatorname{vector}(("<<renderer.c.x<<","<<renderer.c.y<<","<<renderer.c.z<<"),("<<cf.x<<","<<cf.y<<","<<cf.z<<"))" << std::endl;
    Vector3 fx = crossProd(renderer.f, Vector3(0, 0, 1)).normalized();
    Vector3 fy = crossProd(fx, renderer.f).normalized() * -1;
    float pixToM = renderer.xSize / renderer.xPixls;
    //file << "s = [";
    std::vector<Ray> rays = {};
    for (int y = 0; y < renderer.yPixls; y++) {
        for (int x = 0; x < renderer.xPixls; x++) {
            Vector3 op = renderer.c + renderer.f + fx * (x - (renderer.xPixls/2) + 0.5) * pixToM + fy * (y - (renderer.yPixls/2) + 0.5) * pixToM; // finds the global position of the x, y Pixel
            rays.push_back(Ray(op, op - renderer.c));
            //file << "("<<op.x<<","<<op.y<<","<<op.z<<")";
            //if (y != renderer.yPixls - 1 || x != renderer.xPixls - 1)
            //    file << ",";
        }
    }
    //file << "]" << std::endl;

    for (int i = 0; i < rays.size(); i++) {
        intersectionData rayCast = renderer.rayCast(rays[i]);
        Vector3 end;
        if (!rayCast.valid) {
            end = rays[i].op + rays[i].a * 5;
            file << "\\operatorname{vector}(("<<rays[i].op.x<<","<<rays[i].op.y<<","<<rays[i].op.z<<"),("<<end.x<<","<<end.y<<","<<end.y<<"))" << std::endl;
            continue;
        }
        end = rays[i].op + rays[i].a * rayCast.t;
        file << "\\operatorname{vector}(("<<rays[i].op.x<<","<<rays[i].op.y<<","<<rays[i].op.z<<"),("<<end.x<<","<<end.y<<","<<end.y<<"))" << std::endl;
        Ray ref = renderer.reflect(rays[i], rayCast);
        end = ref.op + ref.a * 5;
        file << "\\operatorname{vector}(("<<ref.op.x<<","<<ref.op.y<<","<<ref.op.z<<"),("<<end.x<<","<<end.y<<","<<end.y<<"))" << std::endl;
    }
    file << lights[0].radius<<"^2 = (x - "<<lights[0].position.x<<")^2 + (y - "<<lights[0].position.y<<")^2 + (z - "<<lights[0].position.z<<")^2" << std::endl;
    for (int i = 0; i < polies.size(); i++) {
        Vector3 p1 = polies[i].op + polies[i].a;
        Vector3 p2 = polies[i].op + polies[i].b;
        file << "\\operatorname{triangle}(("<<polies[i].op.x<<","<<polies[i].op.y<<","<<polies[i].op.z<<"),("<<p1.x<<","<<p1.y<<","<<p1.z<<"),("<<p2.x<<","<<p2.y<<","<<p2.z<<"))" << std::endl;
    }
    file << sky.radius<<"^2 = (x - "<<sky.position.x<<")^2 + (y - "<<sky.position.y<<")^2 + (z - "<<sky.position.z<<")^2";
    file.close();

    writePNG("result.png", renderer.render());
    //writePNG("intensityGrayscale.png", renderer.intensityGrayscale);
    std::cout << "render completed!" << std::endl;
}