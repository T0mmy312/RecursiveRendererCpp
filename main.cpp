#include <iostream>
#include "include/renderer.h"

int main()
{
    std::vector<Polygon> polies = cube(Vector3(), Vector3(1, 0, 0), Color(255, 255, 255), 1, 1, 1);
    std::vector<Light> lights = {Light(Vector3(3, 0, 5), 1)};
    Renderer renderer(polies, lights, Color(30, 30, 30), Vector3(0, -5, 0), Vector3(0, 0.2f, 0), 800, 600, 0.5f, 12, 2, 0.2f);
    renderer.render();
    writePNG("C:/Users/Thomas Shuttleworth/Desktop/RecursiveRendererCpp/result.png", renderer.result);
}