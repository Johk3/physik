#ifndef RENDERING_H
#define RENDERING_H

#include "object.h"
#include <GLFW/glfw3.h>
#include <vector>
#include "settings.h"

// Drawing functions
void drawSquare(const double x, const double y, const double r, const double g, const double b, const double alpha, const double size);
void drawCircle(const double x, const double y, const double radius, const double r, const double g, const double b, const double alpha);
void drawObject(const Object& obj);
void render_screen(const std::vector<Object>& all_objects, GLFWwindow* window);
void setup_input_callbacks(GLFWwindow* window);
void process_input(GLFWwindow* window);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void worldToScreen(double worldX, double worldY, double& screenX, double& screenY);

#endif // RENDERING_H