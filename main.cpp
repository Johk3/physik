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

    // Initialize settings
    Settings::initialize();

    // Initialize objects
    std::vector<Object> allObjects = get_objects();

    SpatialGrid grid(2.0, 2.0, 0.1);  // Assuming world size is 2x2 (-1 to 1 in both dimensions)
    ThreadPool pool(std::thread::hardware_concurrency());

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
        render_screen(allObjects, simulationWindow);

        // Update FPS
        updateFPS();

        // Render control panel
        glfwMakeContextCurrent(controlWindow);
        renderControlPanel(controlWindow);
        glfwSwapBuffers(controlWindow);

        // Poll for and process events
        glfwPollEvents();

        // Set the refresh rate
        std::this_thread::sleep_for(std::chrono::milliseconds(1000/Settings::g_refreshRate));
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