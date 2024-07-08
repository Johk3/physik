#include <GLFW/glfw3.h>
#include "linmath.h"
#include <cmath>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <future>
#include <iterator>
#include <thread>
#include <vector>
#include <memory>
#include <mutex>
#include <queue>
#include <bits/stl_algo.h>
#include "settings.h"
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <stdexcept>
#include <unordered_set>

struct Vector2 {
    double x, y;
    // Vector operations which can be performed locally
    Vector2 operator-(const Vector2& other) const {
        return {x - other.x, y - other.y};
    }

    Vector2 operator+(const Vector2& other) const {
        return {x + other.x, y + other.y};
    }

    Vector2 operator*(const float scalar) const {
        return {x * scalar, y * scalar};
    }

    Vector2 operator/(const float scalar) const {
        return {x / scalar, y / scalar};
    }

    // Length of vector, modulus
    [[nodiscard]] double length() const {
        return std::sqrt(x*x + y*y);
    }

    // Normal
    [[nodiscard]] Vector2 normal() const {
        return {x / length(), y / length()};
    }

    // Dot product
    [[nodiscard]] double dot(const Vector2& other) const {
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
    static constexpr double TRAIL_SPACING = SPACING;  // Fixed distance between trail dots
    static constexpr size_t MAX_TRAIL_LENGTH = TRAIL_LENGTH;  // Maximum number of trail dots
};

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

// Thread Pool implementation

class ThreadPool {
public:
    ThreadPool(size_t numThreads) : stop(false) {
        for(size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back(
                [this] {
                    while(true) {
                        std::function<void()> task;
                        {
                            std::unique_lock<std::mutex> lock(this->queue_mutex);
                            this->condition.wait(lock,
                                [this] { return this->stop || !this->tasks.empty(); });
                            if(this->stop && this->tasks.empty())
                                return;
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                        }
                        task();
                    }
                }
            );
        }
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::invoke_result<F, Args...>::type> {
        using return_type = typename std::invoke_result<F, Args...>::type;

        auto task = std::make_shared< std::packaged_task<return_type()> >(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
            );

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if(stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");
            tasks.emplace([task](){ (*task)(); });
        }
        condition.notify_one();
        return res;
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for(std::thread &worker: workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
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

    if (distance >= Object::TRAIL_SPACING) {
        // Calculate how many new dots we need to add
        const int numNewDots = static_cast<int>(distance / Object::TRAIL_SPACING);

        for (int i = 0; i < numNewDots; ++i) {
            const double t = (i + 1) * Object::TRAIL_SPACING / distance;
            Vector2 newDotPosition = {
                lastPosition.x + diff.x * t,
                lastPosition.y + diff.y * t
            };
            obj.trail.push_back(newDotPosition);

            // Remove oldest dot if we exceed the maximum trail length
            if (obj.trail.size() > Object::MAX_TRAIL_LENGTH) {
                obj.trail.erase(obj.trail.begin());
            }
        }
    }
}

// Check if two circular objects are overlapping
bool checkCollision(const Object& obj1, const Object& obj2) {
    const double distance = (obj1.position - obj2.position).length();
    return distance < (obj1.radius + obj2.radius);
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
    float impulseScalar = - (1 + BOUNCE_FACTOR) * velocityAlongNormal;
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
        Vector2 separationVector = normal * (overlap * 0.5);

        obj1.position.x -= separationVector.x * obj2.mass/(obj1.mass + obj2.mass);
        obj1.position.y -= separationVector.y * obj2.mass/(obj1.mass + obj2.mass);
        obj2.position.x += separationVector.x * obj1.mass/(obj1.mass + obj2.mass);
        obj2.position.y += separationVector.y * obj1.mass/(obj1.mass + obj2.mass);
    }
}

std::mutex cout_mutex;
std::atomic<bool> simulation_running(true);

// Simplified parallel gravity calculation
Vector2 calculateGravity(const Object& object1, const std::vector<Object>& objects, size_t start, size_t end) {
    Vector2 totalAcceleration = {0.0f, 0.0f};

    for (size_t i = start; i < end; ++i) {
        const auto& object2 = objects[i];
        if (&object1 == &object2) continue;

        Vector2 diff = object2.position - object1.position;
        double r2 = diff.x * diff.x + diff.y * diff.y;
        double r = std::sqrt(r2);

        if (r < 1e-6f) continue;

        double force = 6.67430e-11 * object2.mass / (r2 * r);
        double MAX_FORCE = 10000.0f;
        force = std::min(force, MAX_FORCE);

        totalAcceleration = totalAcceleration + diff * force;
    }

    return totalAcceleration;
}

void updateObjectsChunk(std::vector<Object>& objects, size_t start, size_t end, double delta_time) {
    for (size_t i = start; i < end; ++i) {
        Object& obj = objects[i];

        // Calculate gravity for this object
        Vector2 acceleration = calculateGravity(obj, objects, 0, objects.size());

        // Update object's state
        obj.acceleration = acceleration;
        obj.velocity = obj.velocity + obj.acceleration * delta_time;
        obj.position = obj.position + obj.velocity * delta_time;

        updateObjectTrail(obj);
    }
}

void updateSimulation(std::vector<Object>& allObjects, SpatialGrid& grid, ThreadPool& pool, double delta_time) {
    const size_t numObjects = allObjects.size();
    const size_t numThreads = std::thread::hardware_concurrency();
    const size_t chunkSize = std::max(size_t(1), numObjects / numThreads);

    std::vector<std::future<void>> futures;

    for (size_t i = 0; i < numObjects; i += chunkSize) {
        size_t end = std::min(i + chunkSize, numObjects);
        futures.push_back(pool.enqueue([&allObjects, i, end, delta_time]() {
            for (size_t j = i; j < end; ++j) {
                Object& obj = allObjects[j];
                Vector2 acceleration = calculateGravity(obj, allObjects, 0, allObjects.size());
                obj.acceleration = acceleration;
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
    for (int i = 0; i <= CIRCLE_SEGMENTS; i++) {
        const double theta = 2.0 * M_PI * static_cast<double>(i) / static_cast<double>(CIRCLE_SEGMENTS);
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

//Used to set the color of an object
void set_color(Object& obj, double r, double g, double b) {

    obj.r = r;
    obj.g = g;
    obj.b = b;

}

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        return -1;
    }

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(1000, 1000, "GRAVITY", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);
    // Initialize objects

    std::vector<Object> allObjects;

    // Create and add objects
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 25; j++) {
            allObjects.emplace_back(Vector2{-0.5 + i * 0.04f, j * 0.04f}, Vector2{1.5f, 0.0f}, 1e5, 0.1, 1.0, 1.0, 1.0);
        }
    }
    allObjects.emplace_back(Vector2{0.02, -0.5f}, Vector2{0.0f, 0.0f}, 5e11, 5e1, 0.0, 0.0, 1.0);


    double constexpr delta_time = 1.0f / REFRESH_RATE;

    SpatialGrid grid(2.0, 2.0, 0.1);  // Assuming world size is 2x2 (-1 to 1 in both dimensions)
    ThreadPool pool(std::thread::hardware_concurrency()); // Hardware conc returns the amount of supported concurrent threads


    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        updateSimulation(allObjects, grid, pool, delta_time);

        // Draw all objects
        for (const auto& obj : allObjects) {
            drawObject(obj);
        }

        // Swap the back buffer with the front
        glfwSwapBuffers(window);
        // Listen for any events
        glfwPollEvents();

        // Set the refresh rate
        std::this_thread::sleep_for(std::chrono::milliseconds(1000/REFRESH_RATE));
    }

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
