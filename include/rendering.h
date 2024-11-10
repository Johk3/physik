#ifndef RENDERING_H
#define RENDERING_H

#include "object.h"
#include <GLFW/glfw3.h>
#include <vector>

// Function to transform world coordinates to screen coordinates
void worldToScreen(double worldX, double worldY, double& screenX, double& screenY);

// -- Objects which can be drawn
void drawSquare(const double x, const double y, const double r, const double g, const double b, const double alpha, const double size);

void drawCircle(const double x, const double y, const double radius, const double r, const double g, const double b, const double alpha);

// Draws an object and its trail on the screen
void drawObject(const Object& obj);

// Renders the screen window based on objects on it
void render_screen(const std::vector<Object>& all_objects, GLFWwindow* window);

// Mouse button callback
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);

// Mouse movement callback
void mouse_callback(GLFWwindow* window, double xpos, double ypos);

// Scroll callback
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

// Function to handle keyboard input
void process_input(GLFWwindow* window);

// Function to set up input callbacks
void setup_input_callbacks(GLFWwindow* window);

#endif // RENDERING_H