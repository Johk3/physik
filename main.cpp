#include <GLFW/glfw3.h>
#include "linmath.h"
#include <cmath>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <future>
#include <iterator>
#include <thread>
#include <future>
#include <vector>
#include <memory>
#include "settings.h"
#include "ThreadPool.h"
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <stdexcept>
#include <unordered_set>
#include <imgui.h>
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

struct Vector2 {
    double x, y;
    // Vector operations which can be performed locally
    // constexpr is for compile-time constant use
    constexpr Vector2 operator-(const Vector2& other) const {
        return {x - other.x, y - other.y};
    }

    constexpr Vector2 operator+(const Vector2& other) const {
        return {x + other.x, y + other.y};
    }

    constexpr Vector2 operator*(const float scalar) const {
        return {x * scalar, y * scalar};
    }

    constexpr Vector2 operator/(const float scalar) const {
        return {x / scalar, y / scalar};
    }

    // Length of vector, modulus
    [[nodiscard]] constexpr double length() const {
        return std::sqrt(x*x + y*y);
    }

    // Normal
    [[nodiscard]] constexpr Vector2 normal() const {
        return {x / length(), y / length()};
    }

    // Dot product
    [[nodiscard]] constexpr double dot(const Vector2& other) const {
        return (x * other.x + y * other.y);
    }

};

//Physical Object
class Object {
public:
    Object(const Vector2 pos, const Vector2 vel, const double m, const double d, const double red, const double green, const double blue)
        : position(pos), velocity(vel), mass(m), density(d), r(red), g(green), b(blue) {
        radius = (std::cbrt((3 * mass) / (4 * 3.141592 * density)) / (1000 * SCALE_FACTOR));
    }
    // Physical location and movement. Default center, no movement
    Vector2 position = {0.0f, 0.0f};
    Vector2 velocity = {0.0f, 0.0f};
    Vector2 acceleration = {0.0f, 0.0f};

    // Properties. Default 100kg, white dot
    double mass = 100.0;
    double density = 1.0;
    double r = 1.0;
    double g = 1.0;
    double b = 1.0;
    double radius = CONST_RADIUS / (1000 * SCALE_FACTOR);

    // Trail
    std::vector<Vector2> trail;
    float TRAIL_SPACING;  // Fixed distance between trail dots
    size_t MAX_TRAIL_LENGTH;  // Maximum number of trail dots
};

// Initialize static members outside the class
double g_trailSpacing = SPACING;
double g_maxTrailLength = TRAIL_LENGTH;
double g_bounceFactor = BOUNCE_FACTOR;
double g_scaleFactor = SCALE_FACTOR;
double g_circleSegments = CIRCLE_SEGMENTS;
int g_refreshRate = REFRESH_RATE;


// CONTROL PANEL --
// Global variables for FPS calculation
double lastTime = 0.0;
int frameCount = 0;
double fps = 0.0;

// Global variables for adjustable settings
float g_spacing = SPACING;
int g_trailLength = TRAIL_LENGTH;

// Function to calculate and update FPS
void updateFPS() {
    double currentTime = glfwGetTime();
    frameCount++;

    if (currentTime - lastTime >= 1.0) {
        fps = static_cast<double>(frameCount) / (currentTime - lastTime);
        frameCount = 0;
        lastTime = currentTime;
    }
}

// Function to render the control panel
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

    // TODO: The draggability doesnt work yet
    // Add a draggable area at the top of the window
    ImGui::TextColored(ImVec4(0.5f, 1.5f, 0.5f, 1.0f), "Control Panel");
    if (ImGui::IsItemActive() && ImGui::IsMouseDragging(0)) {
        int x, y;
        glfwGetWindowPos(controlWindow, &x, &y);
        ImVec2 delta = ImGui::GetIO().MouseDelta;
        glfwSetWindowPos(controlWindow, x + delta.x, y + delta.y);
    }

    ImGui::Separator();


    ImGui::Text("FPS: %.1f", fps);

    // Trail Spacing
    ImGui::InputDouble("Trail Spacing", &g_trailSpacing, 0.001, 0.01, "%.3f");
    g_trailSpacing = std::max(0.001, std::min(0.1, g_trailSpacing));

    // Trail Length
    ImGui::InputDouble("Trail Length", &g_maxTrailLength, 1.0, 10.0);
    g_maxTrailLength = std::max(10.0, std::min(1000.0, g_maxTrailLength));

    // Bounce Factor
    ImGui::InputDouble("Bounce Factor", &g_bounceFactor, 0.01, 0.1, "%.2f");
    g_bounceFactor = std::max(0.0, std::min(1.0, g_bounceFactor));

    // Scale Factor
    ImGui::InputDouble("Scale Factor", &g_scaleFactor, 0.1, 1.0, "%.2f");
    g_scaleFactor = std::max(0.1, std::min(100.0, g_scaleFactor));

    // Circle Segments
    ImGui::InputDouble("Circle Segments", &g_circleSegments, 1.0, 10.0);
    g_circleSegments = std::max(3.0, std::min(100.0, g_circleSegments));

    // Refresh Rate
    ImGui::InputInt("Refresh Rate", &g_refreshRate, 100.0, 1000.0);
    g_refreshRate = std::max(1, std::min(100000, g_refreshRate));

    ImGui::Text("Current Trail Spacing: %.3f", g_trailSpacing);
    ImGui::Text("Current Trail Length: %.0f", g_maxTrailLength);
    ImGui::Text("Current Bounce Factor: %.2f", g_bounceFactor);
    ImGui::Text("Current Scale Factor: %.2f", g_scaleFactor);
    ImGui::Text("Current Circle Segments: %.0f", g_circleSegments);
    ImGui::Text("Current Refresh Rate: %d", g_refreshRate);

    ImGui::End();

    ImGui::Render();
    int display_w, display_h;
    glfwGetFramebufferSize(controlWindow, &display_w, &display_h);
    glViewport(0, 0, display_w, display_h);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}
// CONTROL PANEL --



// -- OPTIMIZATIONS --
// Optimized SpatialGrid
class SpatialGrid {
private:
    struct Cell {
        std::vector<Object*> objects;
    };

    std::vector<Cell> grid;
    int gridWidth, gridHeight;
    double cellSize;
    double worldWidth, worldHeight;

public:
    SpatialGrid(double worldWidth, double worldHeight, double cellSize)
        : worldWidth(worldWidth), worldHeight(worldHeight), cellSize(cellSize) {
        gridWidth = static_cast<int>(std::ceil(worldWidth / cellSize));
        gridHeight = static_cast<int>(std::ceil(worldHeight / cellSize));
        grid.resize(gridWidth * gridHeight);
    }

    void clear() {
        for (auto& cell : grid) {
            cell.objects.clear();
        }
    }

    void insert(Object* obj) {
        // Calculate the range of cells this object could occupy
        int minX = static_cast<int>((obj->position.x - obj->radius + worldWidth/2) / cellSize);
        int maxX = static_cast<int>((obj->position.x + obj->radius + worldWidth/2) / cellSize);
        int minY = static_cast<int>((obj->position.y - obj->radius + worldHeight/2) / cellSize);
        int maxY = static_cast<int>((obj->position.y + obj->radius + worldHeight/2) / cellSize);

        // Clamp to grid boundaries
        minX = std::max(0, std::min(minX, gridWidth - 1));
        maxX = std::max(0, std::min(maxX, gridWidth - 1));
        minY = std::max(0, std::min(minY, gridHeight - 1));
        maxY = std::max(0, std::min(maxY, gridHeight - 1));

        // Insert the object into all cells it occupies
        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                grid[y * gridWidth + x].objects.push_back(obj);
            }
        }
    }

    std::vector<Object*> getNeighbors(const Object& obj) {
        std::vector<Object*> neighbors;
        std::unordered_set<Object*> uniqueNeighbors;  // To avoid duplicates

        // Calculate the range of cells to check
        int minX = static_cast<int>((obj.position.x - obj.radius - cellSize + worldWidth/2) / cellSize);
        int maxX = static_cast<int>((obj.position.x + obj.radius + cellSize + worldWidth/2) / cellSize);
        int minY = static_cast<int>((obj.position.y - obj.radius - cellSize + worldHeight/2) / cellSize);
        int maxY = static_cast<int>((obj.position.y + obj.radius + cellSize + worldHeight/2) / cellSize);

        // Clamp to grid boundaries
        minX = std::max(0, std::min(minX, gridWidth - 1));
        maxX = std::max(0, std::min(maxX, gridWidth - 1));
        minY = std::max(0, std::min(minY, gridHeight - 1));
        maxY = std::max(0, std::min(maxY, gridHeight - 1));

        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                for (Object* neighbor : grid[y * gridWidth + x].objects) {
                    if (neighbor != &obj && uniqueNeighbors.insert(neighbor).second) {
                        neighbors.push_back(neighbor);
                    }
                }
            }
        }

        return neighbors;
    }
};


// ----
void updateObjectTrail(Object& obj) {
    if (obj.trail.empty()) {
        obj.trail.push_back(obj.position);
        return;
    }

    const Vector2 lastPosition = obj.trail.back();
    const Vector2 diff = obj.position - lastPosition;
    double distance = diff.length();

    if (distance >= obj.TRAIL_SPACING) {
        // Calculate how many new dots we need to add
        const int numNewDots = static_cast<int>(distance / obj.TRAIL_SPACING);

        for (int i = 0; i < numNewDots; ++i) {
            const double t = (i + 1) * obj.TRAIL_SPACING / distance;
            Vector2 newDotPosition = {
                lastPosition.x + diff.x * t,
                lastPosition.y + diff.y * t
            };
            obj.trail.push_back(newDotPosition);

            // Remove oldest dot if we exceed the maximum trail length
            if (obj.trail.size() > obj.MAX_TRAIL_LENGTH) {
                obj.trail.erase(obj.trail.begin());
            }
        }
    }
}

// Check if two circular objects are overlapping
bool checkCollision(const Object& obj1, const Object& obj2) {
    const double distance = (obj1.position - obj2.position).length();
    return distance <= (obj1.radius + obj2.radius);
}

void handleCollision(Object& obj1, Object& obj2) {
    // Calculate vector from obj1 to obj2
    const Vector2 delta = obj2.position - obj1.position;
    const float distance = delta.length();

    // Normalize the delta vector
    const Vector2 normal = delta.normal();

    // Calculate relative velocity
    const Vector2 relativeVelocity = obj2.velocity - obj1.velocity;

    // Calculate relative velocity along the normal using dot product
    float velocityAlongNormal = relativeVelocity.dot(normal);

    // Do not resolve if velocities are separating
    if (velocityAlongNormal > 0)
        return;

    // Calculate impulse scalar
    float impulseScalar = - (1 + g_bounceFactor) * velocityAlongNormal;
    impulseScalar /= 1 / obj1.mass  +  1 / obj2.mass;

    // Apply impulse
    Vector2 impulse = normal * impulseScalar;

    // Update velocities
    obj1.velocity.x -= impulse.x / obj1.mass;
    obj1.velocity.y -= impulse.y / obj1.mass;
    obj2.velocity.x += impulse.x / obj2.mass;
    obj2.velocity.y += impulse.y / obj2.mass;

    // Separate the objects to prevent overlapping
    const float overlap = (obj1.radius + obj2.radius) - distance;
    if (overlap > 0) {
        Vector2 separationVector = normal * overlap;

        obj1.position.x -= separationVector.x * obj2.mass/(obj1.mass + obj2.mass);
        obj1.position.y -= separationVector.y * obj2.mass/(obj1.mass + obj2.mass);
        obj2.position.x += separationVector.x * obj1.mass/(obj1.mass + obj2.mass);
        obj2.position.y += separationVector.y * obj1.mass/(obj1.mass + obj2.mass);
    }
}

std::mutex cout_mutex;
std::atomic<bool> simulation_running(true);

// Simplified parallel gravity calculation
void calculate_gravity(Object& object1, const std::vector<Object>& objects, size_t start, size_t end) {

    // Reset acceleration
    object1.acceleration = {0.0, 0.0};

    for (size_t i = start; i < end; ++i) {
        const auto& object2 = objects[i];
        if (&object1 == &object2) continue;
        if (object2.mass < MIN_GRAVITY_MASS) continue;

        Vector2 diff = object2.position - object1.position;
        double r2 = diff.x * diff.x + diff.y * diff.y;

        if (r2 < 1e-12) continue;  // Avoid division by zero

        double r = std::sqrt(r2);

        // Avoid division by zero
        if (r < 1e-6f) continue;

        double force = 6.67430e-11 * object2.mass / (r2 * r);

        double MAX_FORCE = 10000.0f;
        force = std::min(force, MAX_FORCE);

        object1.acceleration = object1.acceleration + diff * force;
    }
}

void updateSimulation(std::vector<Object>& allObjects, SpatialGrid& grid, ThreadPool& pool, double delta_time) {

    const size_t numObjects = allObjects.size();
    // Determine the number of threads available on the system
    const size_t numThreads = std::thread::hardware_concurrency();
    // Calculate the number of objects each thread will process
    // Ensure at least 1 object per chunk to avoid division by zero
    const size_t chunkSize = std::max(size_t(1), numObjects / numThreads);

    std::vector<std::future<void>> futures;

    // Divide the objects into chunks and process each chunk in parallel
    for (size_t i = 0; i < numObjects; i += chunkSize) {
        // Calculate the end index for this chunk, ensuring we don't go past the array bounds
        size_t end = std::min(i + chunkSize, numObjects);
        // Enqueue a task to process this chunk of objects
        futures.push_back(pool.enqueue([&allObjects, i, end, delta_time]() {
            // Process each object in this chunk
            for (size_t j = i; j < end; ++j) {
                Object& obj = allObjects[j];
                calculate_gravity(obj, allObjects, 0, allObjects.size());
                obj.velocity = obj.velocity + obj.acceleration * delta_time;
                obj.position = obj.position + obj.velocity * delta_time;
                updateObjectTrail(obj);
            }
        }));
    }

    for (auto& future : futures) {
        future.get();
    }

    // Clear and rebuild the spatial grid
    grid.clear();
    for (auto& obj : allObjects) {
        grid.insert(&obj);
    }

    // Handle collisions
    for (size_t i = 0; i < allObjects.size(); i++) {
        auto neighbors = grid.getNeighbors(allObjects[i]);
        for (auto* neighbor : neighbors) {
            if (checkCollision(allObjects[i], *neighbor)) {
                handleCollision(allObjects[i], *neighbor);
            }
        }
    }
}

// -- Objects which can be drawn
void drawSquare(const double x, const double y, const double r, const double g, const double b, const double alpha, const double size) {
    mat4x4 mvp;
    mat4x4_identity(mvp);
    mat4x4_translate(mvp, x, y, 0.0f);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMultMatrixf((const GLfloat*)mvp);

    glColor4f(r, g, b, alpha);
    glPointSize(size);
    glBegin(GL_POINTS);
    glVertex2f(0.0f, 0.0f);
    glEnd();
}

void drawCircle(const double x, const double y, const double radius, const double r, const double g, const double b, const double alpha) {
    glColor4d(r, g, b, alpha);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    for (int i = 0; i <= g_circleSegments; i++) {
        const double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(g_circleSegments);
        const double dx = radius * std::cos(theta);
        const double dy = radius * std::sin(theta);
        glVertex2d(x + dx, y + dy);
    }
    glEnd();
}
// -- Objects which can be drawn

// Draws an object and its trail on the screen
void drawObject(const Object& obj) {

    // Draw trail
    if constexpr (ENABLE_TRAIL) {
        for (size_t i = 0; i < obj.trail.size(); ++i) {
            double alpha = static_cast<double>(i) / obj.trail.size();
            double trailRadius = obj.radius * 0.5 * alpha;
            drawCircle(obj.trail[i].x, obj.trail[i].y, trailRadius, obj.r, obj.g, obj.b, alpha);
        }
    };

    // Draw main object
    drawCircle(obj.position.x, obj.position.y, obj.radius, obj.r, obj.g, obj.b, 1.0);
}

// Renders the screen window based on objects on it
void render_screen(const std::vector<Object>& all_objects, GLFWwindow* window) {

    // Clear screen?
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw all objects
    for (const auto& obj : all_objects) {
        drawObject(obj);
    }

    // Swap the back buffer with the front
    glfwSwapBuffers(window);
    // Listen for any events
    glfwPollEvents();

};

std::vector<Object> get_objects() {

    std::vector<Object> allObjects;

    // Create and add objects
    for (int i=0; i < 25; i++) {
        for (int j=0; j < 25; j++) {
            const float r = float(i)/(24.0f * 0.9f) + 0.10f;
            const float g = float(j)/(24.0f * 0.9f) + 0.10f;
            const float b = 1.0f - (float(i+j)/(48.0f * 0.9f) + 0.10f);

            allObjects.emplace_back(Vector2{-0.5 + i * 0.04f, j * 0.04f}, Vector2{0.0f, 0.0f}, 1e5, 0.1, r, g, b);
        }
    }
    allObjects.emplace_back(Vector2{0.02, -0.5f}, Vector2{0.0f, 0.0f}, 5e11, 5e1, 1.0, 1.0, 1.0);

    return allObjects;

};

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        return -1;
    }

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* simulationWindow = glfwCreateWindow(1000, 1000, "GRAVITY", nullptr, nullptr);
    if (!simulationWindow) {
        glfwTerminate();
        return -1;
    }

    // Create control panel window
    glfwWindowHint(GLFW_DECORATED, GLFW_FALSE);
    GLFWwindow* controlWindow = glfwCreateWindow(350, 350, "Control Panel", nullptr, nullptr);
    glfwWindowHint(GLFW_DECORATED, GLFW_TRUE);  // Reset for future windows
    if (!controlWindow) {
        glfwDestroyWindow(simulationWindow);
        glfwTerminate();
        return -1;
    }

    // Set up ImGui for control panel
    glfwMakeContextCurrent(controlWindow);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    // Dark style
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(controlWindow, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    // Make the window's context current
    glfwMakeContextCurrent(simulationWindow);

    // Initialize objects
    std::vector<Object> allObjects = get_objects();

    SpatialGrid grid(2.0, 2.0, 0.1);  // Assuming world size is 2x2 (-1 to 1 in both dimensions)
    ThreadPool pool(std::thread::hardware_concurrency()); // Hardware conc returns the amount of supported concurrent threads


    // Loop until the user closes the window
    while (!glfwWindowShouldClose(simulationWindow) && !glfwWindowShouldClose(controlWindow)) {
        // Update Simulation
        updateSimulation(allObjects, grid, pool, 1.0 / g_refreshRate);

        // Update all objects with new trail values
        for (auto& obj : allObjects) {
            obj.TRAIL_SPACING = g_trailSpacing;
            obj.MAX_TRAIL_LENGTH = g_maxTrailLength;
            obj.radius = (std::cbrt((3 * obj.mass) / (4 * 3.141592 * obj.density)) / (1000 * g_scaleFactor));
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
        std::this_thread::sleep_for(std::chrono::milliseconds(1000/g_refreshRate));
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
/*
 * arity: 3
 * C(I(Z,Z), C(P1,P3,P2), Z)
 *
 * C(I(Z,Z), C(P1,P3,P2), Z)(x,y,z)
 * I(Z,Z)(C(P1,P3,P2), Z)(x,y,z)
 * I(Z,Z)(P1(P3,P2), Z)(x,y,z)
 * I(Z,Z)(P1(P3,P2)(x,y,z), Z(x,y,z))
 * BASE CASE: I(Z,Z)(0, Z(x,y,z)) = Z(Z(x,y,z)) = 0
 * RECURSIVE CASE: I(Z,Z)(P1(P3,P2)(x,y,z), Z(x,y,z))
 * = Z(I(Z,Z)(P1(P3,P2)(x,y,z), Z(x,y,z), P1(P3,P2)(x,y,z), Z(x,y,z)) = 0,
*/
