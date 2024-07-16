#ifndef RENDERING_H
#define RENDERING_H

#include "object.h"
#include <GLFW/glfw3.h>
#include <vector>
#include "settings.h"
#include <glm/glm.hpp>


// Drawing functions
void drawSquare(const double x, const double y, const double r, const double g, const double b, const double alpha, const double size);
void drawCircle(const double x, const double y, const double radius, const double r, const double g, const double b, const double alpha);
void drawObject(const Object& obj);
void render_screen(const std::vector<Object>& all_objects, GLFWwindow* window, const glm::mat4& view, const glm::mat4& projection);

#endif // RENDERING_H