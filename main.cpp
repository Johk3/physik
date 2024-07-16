#include <imgui.h>
#include "imgui/backends/imgui_impl_glfw.h"
#include "imgui/backends/imgui_impl_opengl3.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <cmath>

#include "include/settings.h"
#include "include/object.h"
#include "include/spatialgrid.h"
#include "include/rendering.h"
#include "include/simulationutils.h"
#include "include/controlpanel.h"
#include "include/ThreadPool.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

// Global variables for camera POV control
glm::mat4 view = glm::mat4(1.0f);
glm::mat4 projection = glm::mat4(1.0f);
float yaw = 0.0f, pitch = 0.0f;
double lastX = 500, lastY = 500;
bool firstMouse = true;
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
glm::vec3 cameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
float cameraDistance = 3.0f;
bool leftMousePressed = false;

void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (leftMousePressed) {
        if (firstMouse) {
            lastX = xpos;
            lastY = ypos;
            firstMouse = false;
        }
        float xoffset = xpos - lastX;
        float yoffset = lastY - ypos;
        lastX = xpos;
        lastY = ypos;

        float sensitivity = 0.1f;
        xoffset *= sensitivity;
        yoffset *= sensitivity;

        yaw += xoffset;
        pitch += yoffset;

        if (pitch > 89.0f) pitch = 89.0f;
        if (pitch < -89.0f) pitch = -89.0f;

        glm::vec3 front;
        front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
        front.y = sin(glm::radians(pitch));
        front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
        glm::vec3 cameraFront = glm::normalize(front);

        cameraPos = cameraTarget - cameraFront * cameraDistance;
        view = glm::lookAt(cameraPos, cameraTarget, cameraUp);
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            leftMousePressed = true;
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        } else if (action == GLFW_RELEASE) {
            leftMousePressed = false;
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            firstMouse = true;
        }
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    cameraDistance -= yoffset * 0.1f;
    cameraDistance = glm::clamp(cameraDistance, 1.0f, 10.0f);

    glm::vec3 cameraFront = glm::normalize(cameraTarget - cameraPos);
    cameraPos = cameraTarget - cameraFront * cameraDistance;
    view = glm::lookAt(cameraPos, cameraTarget, cameraUp);
}

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // Create simulation window
    GLFWwindow* simulationWindow = glfwCreateWindow(1000, 1000, "GRAVITY", nullptr, nullptr);
    if (!simulationWindow) {
        std::cerr << "Failed to create simulation window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glEnable(GL_DEPTH_TEST);

    // Create control panel window
    glfwWindowHint(GLFW_DECORATED, GLFW_FALSE);
    GLFWwindow* controlWindow = glfwCreateWindow(400, 450, "Control Panel", nullptr, nullptr);
    glfwWindowHint(GLFW_DECORATED, GLFW_TRUE);  // Reset for future windows
    if (!controlWindow) {
        std::cerr << "Failed to create control window" << std::endl;
        glfwDestroyWindow(simulationWindow);
        glfwTerminate();
        return -1;
    }

    // Set up ImGui for control panel
    glfwMakeContextCurrent(controlWindow);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(controlWindow, true);
    ImGui_ImplOpenGL3_Init("#version 110");

    // Make the simulation window's context current
    glfwMakeContextCurrent(simulationWindow);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
    // Initialize settings
    Settings::initialize();

    // Initialize opencl
    #ifdef USE_GPU
        initOpenCL();
    #endif

    // Initialize objects
    std::vector<Object> allObjects = get_objects();

    SpatialGrid grid(1.0f);
    ThreadPool pool(std::thread::hardware_concurrency());

    // Set up mouse callback
    glfwSetCursorPosCallback(simulationWindow, mouse_callback);
    glfwSetInputMode(simulationWindow, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // Set up projection matrix
    projection = glm::perspective(glm::radians(45.0f), 1000.0f / 1000.0f, 0.1f, 100.0f);

    // Set up mouse callbacks
    glfwSetCursorPosCallback(simulationWindow, mouse_callback);
    glfwSetMouseButtonCallback(simulationWindow, mouse_button_callback);
    glfwSetScrollCallback(simulationWindow, scroll_callback);

    // Initial view matrix setup
    view = glm::lookAt(cameraPos, cameraTarget, cameraUp);
    // Main loop
    while (!glfwWindowShouldClose(simulationWindow) && !glfwWindowShouldClose(controlWindow)) {
        if (Settings::g_simulate) {
            updateSimulation(allObjects, grid, pool, 1.0 / Settings::g_refreshRate);

            // Update all objects with new trail values
            for (auto& obj : allObjects) {
                obj.TRAIL_SPACING = Settings::g_trailSpacing;
                obj.MAX_TRAIL_LENGTH = Settings::g_maxTrailLength;
                obj.radius = (std::cbrt((3 * obj.mass) / (4 * 3.141592 * obj.density)) / (1000 * Settings::g_scaleFactor));
            }
        }

        // Render simulation window
        glfwMakeContextCurrent(simulationWindow);
        render_screen(allObjects, simulationWindow, view, projection);

        // Update FPS
        updateFPS();

        // Render control panel
        glfwMakeContextCurrent(controlWindow);
        renderControlPanel(controlWindow);
        glfwSwapBuffers(controlWindow);

        // Poll for and process events
        glfwPollEvents();
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(simulationWindow);
    glfwDestroyWindow(controlWindow);
    glfwTerminate();
    return 0;
}