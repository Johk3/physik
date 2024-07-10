// src/ControlPanel.cpp
#include "../include/ControlPanel.h"
#include "../include/Settings.h"
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <string>
#include <GLFW/glfw3.h>

// Global variables for FPS calculation
double lastTime = 0.0;
int frameCount = 0;
double fps = 0.0;

void updateFPS() {
    double currentTime = glfwGetTime();
    frameCount++;

    if (currentTime - lastTime >= 1.0) {
        fps = static_cast<double>(frameCount) / (currentTime - lastTime);
        frameCount = 0;
        lastTime = currentTime;
    }
}

std::string ui_toggle_text(bool value) {
    return value ? "ON" : "OFF";
}

void renderControlPanel(GLFWwindow* controlWindow) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    int width, height;
    glfwGetWindowSize(controlWindow, &width, &height);

    ImGui::SetNextWindowPos(ImVec2(0, 0));
    ImGui::SetNextWindowSize(ImVec2(width, height));
    ImGui::Begin("Control Panel", nullptr,
                 ImGuiWindowFlags_NoTitleBar |
                 ImGuiWindowFlags_NoResize |
                 ImGuiWindowFlags_NoMove |
                 ImGuiWindowFlags_NoCollapse |
                 ImGuiWindowFlags_NoSavedSettings |
                 ImGuiWindowFlags_NoBringToFrontOnFocus);

    ImGui::TextColored(ImVec4(0.5f, 1.5f, 0.5f, 1.0f), "Control Panel");
    ImGui::Separator();

    ImGui::Text("FPS: %.1f", fps);

    // Gravity
    ImGui::InputDouble("Gravity", &Settings::g_G, 0.001, 0.01, "%.3f");

    // MAX FORCE
    ImGui::InputDouble("Max Force", &Settings::g_MAX_FORCE, 100, 10000, "%.3f");

    // Trail Spacing
    ImGui::InputDouble("Trail Spacing", &Settings::g_trailSpacing, 0.001, 0.01, "%.3f");

    // Trail Length
    ImGui::InputDouble("Trail Length", &Settings::g_maxTrailLength, 1.0, 10.0);

    // Bounce Factor
    ImGui::InputDouble("Bounce Factor", &Settings::g_bounceFactor, 0.1, 1.0, "%.2f");

    // Scale Factor
    ImGui::InputDouble("Scale Factor", &Settings::g_scaleFactor, 0.5, 1.0, "%.2f");

    // Circle Segments
    ImGui::InputInt("Circle Segments", &Settings::g_circleSegments);

    // Refresh Rate
    ImGui::InputInt("Refresh Rate", &Settings::g_refreshRate);

    // Enable Trail
    ImGui::Checkbox("Enable Trail", &Settings::g_enableTrail);

    // Enable LOD
    ImGui::Checkbox("Enable LOD", &Settings::g_enableLOD);

    // Simulate
    ImGui::Checkbox("Physics simulation", &Settings::g_simulate);

    // Display current values
    ImGui::Text("Current Gravity: %.3f", Settings::g_G);
    ImGui::Text("Current Max Force: %.0f", Settings::g_MAX_FORCE);
    ImGui::Text("Current Trail Spacing: %.3f", Settings::g_trailSpacing);
    ImGui::Text("Current Trail Length: %.0f", Settings::g_maxTrailLength);
    ImGui::Text("Current Bounce Factor: %.2f", Settings::g_bounceFactor);
    ImGui::Text("Current Scale Factor: %.2f", Settings::g_scaleFactor);
    ImGui::Text("Current Circle Segments: %d", Settings::g_circleSegments);
    ImGui::Text("Current Refresh Rate: %d", Settings::g_refreshRate);

    const std::string sim_text = std::string("Simulation is: ") + ui_toggle_text(Settings::g_simulate);
    const std::string trail_text = std::string("Trail is: ") + ui_toggle_text(Settings::g_enableTrail);
    const std::string lod_text = std::string("LOD is: ") + ui_toggle_text(Settings::g_enableLOD);

    ImGui::Text(trail_text.c_str());
    ImGui::Text(lod_text.c_str());
    ImGui::Text(sim_text.c_str());

    ImGui::End();

    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(controlWindow, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    // After all modifications, update the settings
    Settings::update();
}