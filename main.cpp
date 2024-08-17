#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <random>
#include <thread>
#include <mutex>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <iomanip>
#include <ios>
#include <iostream>

// Constants
const int WINDOW_WIDTH = 1280;
const int WINDOW_HEIGHT = 720;
const float SPHERE_RADIUS = 0.05f;
const float WORLD_RADIUS = 1.0f;
const float GRAVITY = -9.81f;
const int MAX_PARTICLES = 100000;
const int GRID_SIZE = 10;
const float TIME_STEP = 0.016f;
const float SPAWN_INTERVAL = 0.001f; // Spawn a new particle every 0.1 seconds
const int SPAWN_COUNT = 1;
const float SPAWN_VERTICAL_OFFSET = 0.1f; // Vertical distance between spawned objects
const float PARTICLE_MASS = 0.1f; // Mass of each particle in kg
const float MIN_VELOCITY = 0.0f;
const float MAX_VELOCITY = 10.0f;
const float COEFFICIENT_OF_RESTITUTION = 0.90f;
const float CELL_SIZE = (2.0f * WORLD_RADIUS) / GRID_SIZE;

// Camera constants
const float CAMERA_INITIAL_DISTANCE = 5.0f;
const float CAMERA_MIN_DISTANCE = 1.5f;
const float CAMERA_MAX_DISTANCE = 10.0f;
const float CAMERA_ROTATION_SPEED = 0.005f;
const float CAMERA_ZOOM_SPEED = 0.1f;

// Camera struct
struct Camera {
    float distance;
    float theta;
    float phi;
};


// Particle struct
struct Particle {
    glm::vec3 position;
    glm::vec3 oldPosition;
    glm::vec3 acceleration;
    glm::vec3 velocity;
    float radius;
    bool active;
};

// Spatial grid cell
struct GridCell {
    std::vector<int> particleIndices;
    std::mutex mutex;
};

// Global variables
std::vector<Particle> particles(MAX_PARTICLES);
int activeParticles = 0;
float timeSinceLastSpawn = 0.0f;
float lastFrameTime = 0.0f;
bool gravityEnabled = true;
float totalEnergy = 0.0f;
float initialTotalEnergy = 0.0f;
float kineticEnergy = 0.0f;
float potentialEnergy = 0.0f;
float maxEnergy = 100.0f; // Initial max energy, will be updated dynamically
float energyDifference = 0.0f;
glm::vec2 forceFieldPosition(0.0f, 0.0f);
float forceFieldRadius = 200.0f;
float forceFieldStrength = 30.0f;
bool isForceFieldActive = false;

float fps = 0.0f;
float fpsUpdateInterval = 0.5f;  // Update FPS every 0.5 seconds
float fpsLastUpdate = 0.0f;
int frameCount = 0;

GridCell grid[GRID_SIZE][GRID_SIZE][GRID_SIZE];
GLuint VBO, VAO;
GLuint shaderProgram;

Camera camera;
bool leftMousePressed = false;
bool mouseOverButton = false;
ImVec2 buttonPos;
ImVec2 buttonSize;
double lastMouseX, lastMouseY;

// Function declarations
void initOpenGL(); 
void initShaders();
void initParticles();
void renderParticles();
void applyConstraints(Particle& p);
void handleCollisions();
void updateGrid();
void spawnParticles();
glm::ivec3 getGridCoords(const glm::vec3& position);
void addToGrid(int particleIndex, const glm::ivec3& gridCoords);
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos);
void scrollCallback(GLFWwindow* window, double xoffset, double yoffset);
glm::mat4 getViewMatrix();
void updateActiveParticles();
void calculateEnergy();
std::string formatEnergy(float energy);
float calculateParticleEnergy(const Particle& p);

// Shader source code
const char* vertexShaderSource = R"(
    #version 330 core
    layout (location = 0) in vec3 aPos;
    layout (location = 1) in float aVelocity;
    uniform mat4 model;
    uniform mat4 view;
    uniform mat4 projection;
    uniform float pointSize;
    uniform float minVelocity;
    uniform float maxVelocity;
    out vec3 particleColor;
    void main() {
        vec4 viewPos = view * model * vec4(aPos, 1.0);
        gl_Position = projection * viewPos;
        gl_PointSize = pointSize / -viewPos.z;

        // Calculate color based on velocity
        float t = clamp((aVelocity - minVelocity) / (maxVelocity - minVelocity), 0.0, 1.0);
        particleColor = mix(vec3(0.1, 0.0, 1.0), vec3(1.0, 0.0, 0.1), t);
    }
)";

const char* fragmentShaderSource = R"(
    #version 330 core
    in vec3 particleColor;
    out vec4 FragColor;
    void main() {
        vec2 circCoord = 2.0 * gl_PointCoord - 1.0;
        if (dot(circCoord, circCoord) > 1.0) {
            discard;
        }
        FragColor = vec4(particleColor, 1.0);
    }
)";

void initImGui(GLFWwindow* window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;   // Enable Gamepad Controls
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");
}


void cleanupImGui() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void updateFPS(float deltaTime) {
    frameCount++;
    fpsLastUpdate += deltaTime;

    if (fpsLastUpdate >= fpsUpdateInterval) {
        fps = frameCount / fpsLastUpdate;
        frameCount = 0;
        fpsLastUpdate = 0;
    }
}

void renderImGui() {
    ImGui::SetNextWindowPos(ImVec2(10, 10));
    ImGui::SetNextWindowBgAlpha(0.35f);
    ImGui::Begin("Stats and Controls", nullptr, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav);

    ImGui::Text("FPS: %.1f", fps);
    ImGui::Text("Objects: %d", activeParticles);

    // Energy bars
    float barWidth = 200.0f;
    float barHeight = 20.0f;
    ImVec2 cursor = ImGui::GetCursorScreenPos();

    // Total Energy
    ImGui::Text("Total Energy:");
    float totalEnergyNormalized = std::min(totalEnergy / maxEnergy, 1.0f);
    ImGui::ProgressBar(totalEnergyNormalized, ImVec2(barWidth, barHeight));
    ImGui::SameLine();
    ImGui::Text("%s", formatEnergy(totalEnergy).c_str());

    // Kinetic Energy
    ImGui::Text("Kinetic Energy:");
    float kineticEnergyNormalized = std::min(kineticEnergy / maxEnergy, 1.0f);
    ImGui::ProgressBar(kineticEnergyNormalized, ImVec2(barWidth, barHeight));
    ImGui::SameLine();
    ImGui::Text("%s", formatEnergy(kineticEnergy).c_str());

    // Potential Energy
    ImGui::Text("Potential Energy:");
    float potentialEnergyNormalized = std::min(potentialEnergy / maxEnergy, 1.0f);
    ImGui::ProgressBar(potentialEnergyNormalized, ImVec2(barWidth, barHeight));
    ImGui::SameLine();
    ImGui::Text("%s", formatEnergy(potentialEnergy).c_str());

    // Initial Total Energy
    ImGui::Text("Initial Total Energy: %s", formatEnergy(initialTotalEnergy).c_str());
    ImGui::Text("Max. allocated energy for system: %s", formatEnergy(maxEnergy).c_str());
    // Energy Conservation
    ImGui::Text("Energy Difference: %s", formatEnergy(energyDifference).c_str());
    if (std::abs(energyDifference) > 0.01f * initialTotalEnergy) {
        ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "Warning: Energy not conserved!");
    } else {
        ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), "Energy conserved within tolerance.");
    }

    // Store button position and size for hit testing
    buttonPos = ImGui::GetCursorScreenPos();
    if (ImGui::Button(gravityEnabled ? "Gravity: ON" : "Gravity: OFF")) {
        gravityEnabled = !gravityEnabled;
    }
    buttonSize = ImGui::GetItemRectSize();

    ImGui::Text("Left-click & drag to rotate camera");
    ImGui::Text("Scroll to zoom");

    ImGui::End();
}

void renderForceField() {
    if (!isForceFieldActive) return;

    glUseProgram(shaderProgram);

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT, 0.1f, 100.0f);
    glm::mat4 view = getViewMatrix();
    glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(forceFieldPosition.x, forceFieldPosition.y, 0.0f));
    model = glm::scale(model, glm::vec3(forceFieldRadius));

    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));

    // Draw force field circle
    const int segments = 64;
    std::vector<glm::vec3> circleVertices;
    for (int i = 0; i <= segments; ++i) {
        float theta = 2.0f * glm::pi<float>() * float(i) / float(segments);
        circleVertices.push_back(glm::vec3(cos(theta), sin(theta), 0.0f));
    }

    GLuint vbo, vao;
    glGenBuffers(1, &vbo);
    glGenVertexArrays(1, &vao);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, circleVertices.size() * sizeof(glm::vec3), circleVertices.data(), GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);

    glDrawArrays(GL_LINE_LOOP, 0, segments);

    glDeleteBuffers(1, &vbo);
    glDeleteVertexArrays(1, &vao);
}


void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    renderParticles();
    renderForceField();
    renderImGui();

    // Render ImGui
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    glfwSwapBuffers(glfwGetCurrentContext());
}


int main() {
    initOpenGL();
    initShaders();

    GLFWwindow* window = glfwGetCurrentContext();
    initImGui(window);  // Initialize ImGui

    // Initialize camera
    camera.distance = CAMERA_INITIAL_DISTANCE;
    camera.theta = glm::radians(45.0f);
    camera.phi = glm::radians(45.0f);

    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPositionCallback);
    glfwSetScrollCallback(window, scrollCallback);

    lastFrameTime = glfwGetTime();

    while (!glfwWindowShouldClose(window)) {
        float currentTime = glfwGetTime();
        float deltaTime = currentTime - lastFrameTime;
        lastFrameTime = currentTime;

        glfwPollEvents();

        // Start the ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        updateFPS(deltaTime);

        timeSinceLastSpawn += deltaTime;
        if (timeSinceLastSpawn >= SPAWN_INTERVAL && activeParticles < MAX_PARTICLES) {
            spawnParticles();
            timeSinceLastSpawn = 0.0f;
        }

        updateActiveParticles();
        updateGrid();
        handleCollisions();
        calculateEnergy();

        renderScene();
    }

    cleanupImGui();  // Cleanup ImGui
    glfwTerminate();
    return 0;
}

void initOpenGL() {
    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Physics Engine", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK) {
        exit(EXIT_FAILURE);
    }

    glEnable(GL_DEPTH_TEST);
    glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

    // Set the clear color to grey
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
}

void initShaders() {
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
    glCompileShader(vertexShader);

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fragmentShader);

    shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, MAX_PARTICLES * (sizeof(glm::vec3) + sizeof(float)), nullptr, GL_DYNAMIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3) + sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(glm::vec3) + sizeof(float), (void*)(sizeof(glm::vec3)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void initParticles() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-0.1f, 0.1f);

    for (int i = 0; i < MAX_PARTICLES; ++i) {
        glm::vec3 position(WORLD_RADIUS - SPHERE_RADIUS, WORLD_RADIUS - SPHERE_RADIUS, WORLD_RADIUS - SPHERE_RADIUS);
        position += glm::vec3(dis(gen), dis(gen), dis(gen));

        Particle p;
        p.position = position;
        p.oldPosition = position;
        p.acceleration = glm::vec3(0.0f, GRAVITY, 0.0f);
        p.radius = SPHERE_RADIUS;

        particles.push_back(p);
    }

    updateGrid();
}

void spawnParticles() {
    if (activeParticles + SPAWN_COUNT <= MAX_PARTICLES) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<float> dis(-0.1f, 0.1f);

        for (int i = 0; i < SPAWN_COUNT; ++i) {
            glm::vec3 spawnPosition(
                -WORLD_RADIUS + SPHERE_RADIUS,
                (i - (SPAWN_COUNT - 1) / 2.0f) * SPAWN_VERTICAL_OFFSET,
                dis(gen) * WORLD_RADIUS
            );
            glm::vec3 initialVelocity(0.1f + dis(gen), dis(gen), dis(gen)); // Shoot towards the right

            particles[activeParticles].position = spawnPosition;
            particles[activeParticles].oldPosition = spawnPosition - initialVelocity * TIME_STEP;
            particles[activeParticles].acceleration = gravityEnabled ? glm::vec3(0.0f, GRAVITY, 0.0f) : glm::vec3(0.0f);
            particles[activeParticles].radius = SPHERE_RADIUS;
            particles[activeParticles].active = true;

            // Calculate and add the energy of the new particle
            float newParticleEnergy = calculateParticleEnergy(particles[activeParticles]);
            initialTotalEnergy += newParticleEnergy;

            // Update maxEnergy if necessary
            if (initialTotalEnergy > maxEnergy) {
                maxEnergy = initialTotalEnergy * 2.0f; // Add 200% buffer
            }

            activeParticles++;
        }
    }
}

float calculateParticleEnergy(const Particle& p) {
    glm::vec3 velocity = (p.position - p.oldPosition) / TIME_STEP;
    float kinetic = 0.5f * PARTICLE_MASS * glm::dot(velocity, velocity);
    float potential = PARTICLE_MASS * std::abs(GRAVITY) * (p.position.y + WORLD_RADIUS);
    return kinetic + potential;
}



void updateActiveParticles() {
    std::vector<std::thread> threads;
    int particlesPerThread = activeParticles / std::thread::hardware_concurrency();

    for (int i = 0; i < std::thread::hardware_concurrency(); ++i) {
        int start = i * particlesPerThread;
        int end = (i == std::thread::hardware_concurrency() - 1) ? activeParticles : (i + 1) * particlesPerThread;

        threads.emplace_back([start, end]() {
            for (int j = start; j < end; ++j) {
                Particle& p = particles[j];
                if (p.active) {
                    // Update position
                    glm::vec3 newPosition = p.position + p.velocity * TIME_STEP + 0.5f * p.acceleration * TIME_STEP * TIME_STEP;

                    // Apply force field effect
                    if (isForceFieldActive) {
                        glm::vec2 particlePos2D(p.position.x, p.position.y);
                        glm::vec2 forceDir = particlePos2D - forceFieldPosition;
                        float distance = glm::length(forceDir);

                        if (distance < forceFieldRadius) {
                            float forceMagnitude = forceFieldStrength * (1.0f - distance / forceFieldRadius);
                            glm::vec2 force = glm::normalize(forceDir) * forceMagnitude;
                            p.velocity += glm::vec3(force.x, force.y, 0.0f) * TIME_STEP;
                        }
                    }

                    // Update velocity
                    p.velocity += p.acceleration * TIME_STEP;

                    // Store old position
                    p.oldPosition = p.position;

                    // Update position
                    p.position = newPosition;

                    // Apply gravity only if it's enabled
                    if (gravityEnabled) {
                        p.acceleration = glm::vec3(0.0f, GRAVITY, 0.0f);
                    } else {
                        p.acceleration = glm::vec3(0.0f, 0.0f, 0.0f);
                    }

                    applyConstraints(p);
                }
            }
        });
    }

    for (auto& thread : threads) {
        thread.join();
    }
}


void applyConstraints(Particle& p) {
    glm::vec3 toCenter = p.position;
    float dist = glm::length(toCenter);

    if (dist > WORLD_RADIUS - p.radius) {
        glm::vec3 n = glm::normalize(toCenter);
        p.position = n * (WORLD_RADIUS - p.radius);

        // Reflect velocity off the sphere surface
        glm::vec3 reflectedVelocity = glm::reflect(p.velocity, n);

        // Apply coefficient of restitution to the normal component of velocity
        glm::vec3 normalVelocity = glm::dot(reflectedVelocity, n) * n;
        glm::vec3 tangentialVelocity = reflectedVelocity - normalVelocity;
        p.velocity = tangentialVelocity + COEFFICIENT_OF_RESTITUTION * normalVelocity;
    }
}



void renderParticles() {
    glUseProgram(shaderProgram);

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)WINDOW_WIDTH / (float)WINDOW_HEIGHT, 0.1f, 100.0f);
    glm::mat4 view = getViewMatrix();

    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniform1f(glGetUniformLocation(shaderProgram, "minVelocity"), MIN_VELOCITY);
    glUniform1f(glGetUniformLocation(shaderProgram, "maxVelocity"), MAX_VELOCITY);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    std::vector<float> particleData;
    for (int i = 0; i < activeParticles; ++i) {
        if (particles[i].active) {
            particleData.push_back(particles[i].position.x);
            particleData.push_back(particles[i].position.y);
            particleData.push_back(particles[i].position.z);
            particleData.push_back(glm::length(particles[i].velocity));
        }
    }

    glBufferSubData(GL_ARRAY_BUFFER, 0, particleData.size() * sizeof(float), particleData.data());

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SPRITE);

    for (int i = 0; i < activeParticles; ++i) {
        if (particles[i].active) {
            glm::mat4 model = glm::translate(glm::mat4(1.0f), particles[i].position);
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));

            float pointSize = particles[i].radius * 1000.0f;
            glUniform1f(glGetUniformLocation(shaderProgram, "pointSize"), pointSize);

            glDrawArrays(GL_POINTS, i, 1);
        }
    }

    glDisable(GL_POINT_SPRITE);
    glDisable(GL_PROGRAM_POINT_SIZE);

    glBindVertexArray(0);
}

void handleCollisions() {
    for (int x = 0; x < GRID_SIZE; ++x) {
        for (int y = 0; y < GRID_SIZE; ++y) {
            for (int z = 0; z < GRID_SIZE; ++z) {
                auto& cell = grid[x][y][z];
                for (size_t i = 0; i < cell.particleIndices.size(); ++i) {
                    for (size_t j = i + 1; j < cell.particleIndices.size(); ++j) {
                        int p1Index = cell.particleIndices[i];
                        int p2Index = cell.particleIndices[j];
                        Particle& p1 = particles[p1Index];
                        Particle& p2 = particles[p2Index];

                        if (p1.active && p2.active) {
                            glm::vec3 collisionAxis = p1.position - p2.position;
                            float dist = glm::length(collisionAxis);

                            if (dist < p1.radius/1.2 + p2.radius/1.2) {
                                glm::vec3 collisionNormal = glm::normalize(collisionAxis);

                                // Calculate relative velocity
                                glm::vec3 relativeVelocity = p1.velocity - p2.velocity;

                                // Calculate relative velocity along the normal
                                float velAlongNormal = glm::dot(relativeVelocity, collisionNormal);

                                // Do not resolve if velocities are separating
                                if (velAlongNormal > 0)
                                    continue;

                                // Calculate restitution
                                float e = COEFFICIENT_OF_RESTITUTION;

                                // Calculate impulse scalar
                                float j = -(1 + e) * velAlongNormal;
                                j /= 1 / PARTICLE_MASS + 1 / PARTICLE_MASS;

                                // Apply impulse
                                glm::vec3 impulse = j * collisionNormal;
                                p1.velocity += impulse / PARTICLE_MASS;
                                p2.velocity -= impulse / PARTICLE_MASS;

                                // Positional correction
                                float percent = 0.2f; // usually 20% to 80%
                                float slop = 0.01f; // usually 0.01 to 0.1
                                glm::vec3 correction = glm::max(dist - (p1.radius + p2.radius) - slop, 0.0f) / (1/PARTICLE_MASS + 1/PARTICLE_MASS) * percent * collisionNormal;
                                p1.position += correction / PARTICLE_MASS;
                                p2.position -= correction / PARTICLE_MASS;
                            }
                        }
                    }
                }
            }
        }
    }
}

void calculateEnergy() {
    kineticEnergy = 0.0f;
    potentialEnergy = 0.0f;

    for (int i = 0; i < activeParticles; ++i) {
        if (particles[i].active) {
            // Kinetic energy: 1/2 * m * v^2
            kineticEnergy += 0.5f * PARTICLE_MASS * glm::dot(particles[i].velocity, particles[i].velocity);

            // Potential energy: m * g * h
            potentialEnergy += PARTICLE_MASS * std::abs(GRAVITY) * (particles[i].position.y + WORLD_RADIUS);
        }
    }

    totalEnergy = kineticEnergy + potentialEnergy;

    // Calculate energy difference
    energyDifference = totalEnergy - initialTotalEnergy;
}


std::string formatEnergy(float energy) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(3);

    if (energy < 1e3) oss << energy << " J";
    else if (energy < 1e6) oss << energy / 1e3 << " kJ";
    else if (energy < 1e9) oss << energy / 1e6 << " MJ";
    else if (energy < 1e12) oss << energy / 1e9 << " GJ";
    else oss << energy / 1e12 << " TJ";

    return oss.str();
}

void updateForceFieldPosition(double xpos, double ypos) {
    int width, height;
    glfwGetWindowSize(glfwGetCurrentContext(), &width, &height);

    // Convert screen coordinates to world coordinates
    glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)width / (float)height, 0.1f, 100.0f);
    glm::mat4 view = getViewMatrix();
    glm::mat4 invVP = glm::inverse(projection * view);

    glm::vec4 screenPos = glm::vec4(
        (2.0f * xpos) / width - 1.0f,
        1.0f - (2.0f * ypos) / height,
        1.0f,
        1.0f
    );

    glm::vec4 worldPos = invVP * screenPos;
    worldPos /= worldPos.w;

    forceFieldPosition = glm::vec2(worldPos.x, worldPos.y);
}


void updateGrid() {
    for (auto& row : grid) {
        for (auto& col : row) {
            for (auto& cell : col) {
                cell.particleIndices.clear();
            }
        }
    }

    for (int i = 0; i < activeParticles; ++i) {
        if (particles[i].active) {
            glm::ivec3 gridCoords = getGridCoords(particles[i].position);
            addToGrid(i, gridCoords);
        }
    }
}

glm::ivec3 getGridCoords(const glm::vec3& position) {
    glm::vec3 normalized = (position + glm::vec3(WORLD_RADIUS)) / (2.0f * WORLD_RADIUS);
    return glm::ivec3(
        glm::clamp(int(normalized.x * GRID_SIZE), 0, GRID_SIZE - 1),
        glm::clamp(int(normalized.y * GRID_SIZE), 0, GRID_SIZE - 1),
        glm::clamp(int(normalized.z * GRID_SIZE), 0, GRID_SIZE - 1)
    );
}

void addToGrid(int particleIndex, const glm::ivec3& gridCoords) {
    grid[gridCoords.x][gridCoords.y][gridCoords.z].particleIndices.push_back(particleIndex);
}

bool isMouseOverButton(double xpos, double ypos) {
    return (xpos >= buttonPos.x && xpos <= buttonPos.x + buttonSize.x &&
            ypos >= buttonPos.y && ypos <= buttonPos.y + buttonSize.y);
}

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            if (!mouseOverButton) {
                leftMousePressed = true;
                lastMouseX = xpos;
                lastMouseY = ypos;
            } else if (mouseOverButton) {
                gravityEnabled = !gravityEnabled;  // Toggle gravity when clicking the button
            }
        }else if (action == GLFW_RELEASE) {
                leftMousePressed = false;
            }
        } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        if (action == GLFW_PRESS) {
            isForceFieldActive = true;
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            updateForceFieldPosition(xpos, ypos);
        } else if (action == GLFW_RELEASE) {
            isForceFieldActive = false;
        }
    }
}

void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
    mouseOverButton = isMouseOverButton(xpos, ypos);

    if (leftMousePressed && !mouseOverButton) {
        double deltaX = xpos - lastMouseX;
        double deltaY = ypos - lastMouseY;

        camera.theta -= CAMERA_ROTATION_SPEED * deltaX;
        camera.phi -= CAMERA_ROTATION_SPEED * deltaY;

        camera.phi = glm::clamp(camera.phi, 0.1f, glm::pi<float>() - 0.1f);

        lastMouseX = xpos;
        lastMouseY = ypos;
    }
    if (isForceFieldActive) {
        updateForceFieldPosition(xpos, ypos);
    }
}

void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) {
        return;  // Let ImGui handle this event
    }

    camera.distance -= CAMERA_ZOOM_SPEED * yoffset;
    camera.distance = glm::clamp(camera.distance, CAMERA_MIN_DISTANCE, CAMERA_MAX_DISTANCE);
}

glm::mat4 getViewMatrix() {
    float x = camera.distance * sin(camera.phi) * sin(camera.theta);
    float y = camera.distance * cos(camera.phi);
    float z = camera.distance * sin(camera.phi) * cos(camera.theta);

    return glm::lookAt(
        glm::vec3(x, y, z),
        glm::vec3(0.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 1.0f, 0.0f)
    );
}

