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
void drawSphere(const Vector3& position, double radius, const double r, const double g, const double b, const double alpha);
void drawSlimArrow(const Vector3& start, const Vector3& direction, float length, float r, float g, float b);
void drawSphereWithShadows(const Vector3& position, double radius, const double r, const double g, const double b, const double alpha);
void drawTriangle(const Vector3& position, double size, const double r, const double g, const double b, const double alpha);
void drawFlatSurface(const Vector3& position, double width, double length, const double r, const double g, const double b, const double alpha);
void drawContainer(const Vector3& position, double width, double height, double depth, const double r, const double g, const double b, const double alpha);
void drawCow(const Vector3& position, double size, const double r, const double g, const double b, const double alpha);
void drawCube(double size);

#endif // RENDERING_H