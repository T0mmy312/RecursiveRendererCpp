#include <iostream>
#include "renderer2.h"

int main()
{
    std::vector<Polygon> polies = cube(Vector3(), Vector3(1, 0, 0), Color(255, 255, 255), 10, 10, 10);
    std::vector<Light> lights = {Light(Vector3(0, 1, 3), Color(), 1000)};
    Renderer renderer(polies, lights, Color(30, 30, 30), Vector3(0, -5, 0), Vector3(0, 0.2f, 0), 200, 150, 1, 3, 2, 0.2f);
    renderer.render();
    writePNG("C:/Users/Thomas Shuttleworth/Desktop/RecursiveRendererCpp/result.png", renderer.result);
    std::cout << "finished!";
}