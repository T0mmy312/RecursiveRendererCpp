#pragma once

#ifndef _PNGWRITE_H_
#define _PNGWRITE_H_

#include <iostream>
#include <vector>
#include "rendererTools.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

bool writePNG(std::string filename, Picture data) {
    const int width = data[0].size();
    const int height = data.size();

    // Create a buffer for the pixel data
    unsigned char* image = new unsigned char[width * height * 4];  // RGBA format

    // Fill the buffer with a solid red color
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            // Calculate the index of the current pixel
            int index = 4 * (y * width + x);

            image[index] = data[y][x].r;    // Red
            image[index + 1] = data[y][x].g;  // Green
            image[index + 2] = data[y][x].b;  // Blue
            image[index + 3] = 255;  // Alpha
        }
    }

    // Save the image as a PNG file
    stbi_write_png(filename.c_str(), width, height, 4, image, width * 4);

    delete[] image;  // Don't forget to free the allocated memory

    return 0;
}

#endif