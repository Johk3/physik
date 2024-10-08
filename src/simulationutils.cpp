#include "../include/simulationutils.h"
#include <iostream>
#include "../include/settings.h"
#ifdef USE_AVX2
#include <immintrin.h>
#endif
#include <fstream>
#ifdef USE_GPU
#include <CL/cl2.hpp>
#include <stdexcept>
// OpenCL global variables
cl::Context context;
cl::CommandQueue queue;
cl::Kernel gravityKernel;
cl::Kernel updatePositionsKernel;
cl::Kernel collisionKernel;
cl::Program program;
#endif

#ifdef USE_AVX2
// Helper function to sum up a __m256d vector
double _mm256_reduce_add_pd(__m256d v) {
    __m128d vlow  = _mm256_castpd256_pd128(v);
    __m128d vhigh = _mm256_extractf128_pd(v, 1);
    vlow  = _mm_add_pd(vlow, vhigh);
    __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
    return _mm_cvtsd_f64(_mm_add_sd(vlow, high64));
}
// SingleInstructionMultipleData Gravity using AVX architecture
void calculate_gravity_simd(Object& object1, const std::vector<Object>& objects, size_t start, size_t end) {
    const __m256d obj1_x = _mm256_set1_pd(object1.position.x);
    const __m256d obj1_y = _mm256_set1_pd(object1.position.y);
    const __m256d epsilon = _mm256_set1_pd(Settings::EPSILON);
    const __m256d max_force = _mm256_set1_pd(Settings::g_MAX_FORCE);
    const __m256d g_const = _mm256_set1_pd(Settings::g_G);

    __m256d acc_x = _mm256_setzero_pd();
    __m256d acc_y = _mm256_setzero_pd();

    for (size_t i = start; i < end; i += 4) {
        __m256d obj2_x, obj2_y, mass;

        if (i + 3 < end) {
            obj2_x = _mm256_setr_pd(objects[i].position.x, objects[i+1].position.x, objects[i+2].position.x, objects[i+3].position.x);
            obj2_y = _mm256_setr_pd(objects[i].position.y, objects[i+1].position.y, objects[i+2].position.y, objects[i+3].position.y);
            mass = _mm256_setr_pd(objects[i].mass, objects[i+1].mass, objects[i+2].mass, objects[i+3].mass);
        } else {
            // Handle the last incomplete batch
            double x_array[4] = {0}, y_array[4] = {0}, m_array[4] = {0};
            for (size_t j = 0; j < end - i; ++j) {
                x_array[j] = objects[i+j].position.x;
                y_array[j] = objects[i+j].position.y;
                m_array[j] = objects[i+j].mass;
            }
            obj2_x = _mm256_loadu_pd(x_array);
            obj2_y = _mm256_loadu_pd(y_array);
            mass = _mm256_loadu_pd(m_array);
        }

        __m256d dx = _mm256_sub_pd(obj2_x, obj1_x);
        __m256d dy = _mm256_sub_pd(obj2_y, obj1_y);
        __m256d r2 = _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy));
        __m256d r = _mm256_sqrt_pd(r2);
        __m256d r3 = _mm256_mul_pd(r, _mm256_mul_pd(r, r));

        __m256d force = _mm256_div_pd(_mm256_mul_pd(g_const, mass), r3);
        force = _mm256_min_pd(force, max_force);

        __m256d mask = _mm256_cmp_pd(r2, epsilon, _CMP_GT_OQ);
        force = _mm256_and_pd(force, mask);

        acc_x = _mm256_add_pd(acc_x, _mm256_mul_pd(force, dx));
        acc_y = _mm256_add_pd(acc_y, _mm256_mul_pd(force, dy));
    }

    // Sum up the vector components
    double acc_x_sum = _mm256_reduce_add_pd(acc_x);
    double acc_y_sum = _mm256_reduce_add_pd(acc_y);

    // Set the final acceleration
    object1.acceleration = {acc_x_sum, acc_y_sum};
}
#endif
// Non-AVX version of gravity calculation
void calculate_gravity_normal(Object& object1, const std::vector<Object>& objects, size_t start, size_t end) {
    Vector2 acc = {0.0, 0.0};
    for (size_t i = start; i < end; ++i) {
        const auto& object2 = objects[i];
        if (&object1 == &object2) continue;
        if (object2.mass < Settings::EPSILON) continue;

        Vector2 diff = object2.position - object1.position;
        double r2 = diff.x * diff.x + diff.y * diff.y;

        if (r2 < 1e-12) continue;

        double r = std::sqrt(r2);
        if (r < 1e-6f) continue;

        double force = 6.67430e-11 * object2.mass / (r2 * r);
        double MAX_FORCE = 10000.0f;
        force = std::min(force, MAX_FORCE);

        acc = acc + diff * force;
    }
    object1.acceleration = acc;
}


// Function to choose between AVX and normal calculation
void calculate_gravity(Object& object1, const std::vector<Object>& objects, size_t start, size_t end) {
#ifdef USE_AVX2
    calculate_gravity_simd(object1, objects, start, end);
#else
    calculate_gravity_normal(object1, objects, start, end);
#endif
}

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
    float impulseScalar = - (1 + Settings::g_bounceFactor) * velocityAlongNormal;
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

// Helper function to read OpenCL kernel source
std::string readKernelSource(const char* filename) {
    std::ifstream file(filename);
    return std::string(std::istreambuf_iterator<char>(file),
                       std::istreambuf_iterator<char>());
}

#ifdef USE_GPU
// Helper function to check OpenCL errors
void checkError(cl_int error, const std::string& message) {
    if (error != CL_SUCCESS) {
        throw std::runtime_error(message + ": " + std::to_string(error));
    }
}

void initializeOpenCL() {
    try {
        // Get available platforms
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if (platforms.empty()) {
            throw std::runtime_error("No OpenCL platforms found");
        }

        // Select the first platform
        cl::Platform platform = platforms[0];

        // Get available devices
        std::vector<cl::Device> devices;
        platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
        if (devices.empty()) {
            throw std::runtime_error("No OpenCL devices found");
        }

        // Select the first device
        cl::Device device = devices[0];

        // Create OpenCL context
        context = cl::Context({device});

        // Create command queue
        queue = cl::CommandQueue(context, device);

        // Read kernel source
        const char* kernelSource = R"(
typedef struct {
    float2 position;
    float2 velocity;
    float2 acceleration;
    float mass;  // Changed from double to float
    float radius;  // Changed from double to float
} Object;

__kernel void calculateGravity(__global Object* objects,
                               const unsigned int numObjects,
                               const float G,
                               const float MAX_FORCE) {
    int i = get_global_id(0);
    if (i >= numObjects) return;

    float2 acc = (float2)(0.0f, 0.0f);

    for (int j = 0; j < numObjects; j++) {
        if (i != j) {
            float2 r = objects[j].position - objects[i].position;
            float dist2 = dot(r, r);
            if (dist2 > 1e-10f) {
                float dist = sqrt(dist2);
                float minDist = objects[i].radius + objects[j].radius;
                dist = max(dist, minDist);  // Prevent objects from getting too close
                float invDist3 = 1.0f / (dist * dist2);
                float force = G * objects[i].mass * objects[j].mass * invDist3;
                force = min(force, MAX_FORCE);
                acc += (force / objects[i].mass) * r;
            }
        }
    }
    objects[i].acceleration = acc;
}

__kernel void updatePositionsAndDetectCollisions(
    __global Object* objects,
    const unsigned int numObjects,
    const float dt,
    const float worldSize,
    const float cellSize,
    __global int* gridCounters,
    __global int* grid,
    const int gridWidth
) {
    int i = get_global_id(0);
    if (i >= numObjects) return;

    // Update velocity and position
    objects[i].velocity += objects[i].acceleration * dt;
    objects[i].position += objects[i].velocity * dt;

    // Wrap around world boundaries
    objects[i].position = (float2)(
        fmod(objects[i].position.x + 1.0f, 2.0f) - 1.0f,
        fmod(objects[i].position.y + 1.0f, 2.0f) - 1.0f
    );

    // Calculate grid cell
    int cellX = (int)((objects[i].position.x + 1.0f) / cellSize);
    int cellY = (int)((objects[i].position.y + 1.0f) / cellSize);
    int cellIndex = cellY * gridWidth + cellX;

    // Add object to grid
    int objectIndex = atomic_inc(&gridCounters[cellIndex]);
    if (objectIndex < 8) { // Limit objects per cell to prevent buffer overflow
        grid[cellIndex * 8 + objectIndex] = i;
    }
}
__kernel void resolveCollisions(
    __global Object* objects,
    const unsigned int numObjects,
    const float worldSize,
    const float cellSize,
    __global int* gridCounters,
    __global int* grid,
    const int gridWidth,
    const float bounceFactor
) {
    int i = get_global_id(0);
    if (i >= numObjects) return;

    int cellX = (int)((objects[i].position.x + 1.0f) / cellSize);
    int cellY = (int)((objects[i].position.y + 1.0f) / cellSize);

    // Check neighboring cells
    for (int offsetY = -1; offsetY <= 1; offsetY++) {
        for (int offsetX = -1; offsetX <= 1; offsetX++) {
            int neighborX = (cellX + offsetX + gridWidth) % gridWidth;
            int neighborY = (cellY + offsetY + gridWidth) % gridWidth;
            int neighborCellIndex = neighborY * gridWidth + neighborX;

            int neighborCount = gridCounters[neighborCellIndex];
            for (int k = 0; k < min(neighborCount, 8); k++) {
                int j = grid[neighborCellIndex * 8 + k];
                if (i == j) continue;

                float2 r = objects[j].position - objects[i].position;
                r = r - round(r / worldSize) * (float2)(worldSize, worldSize); // Adjust for wraparound
                float dist = length(r) + 1e-6f;
                float minDist = objects[i].radius + objects[j].radius;

                if (dist < minDist) {
                    float2 normal = r / (dist + 1e-6f);
                    float2 relVel = objects[j].velocity - objects[i].velocity;
                    float velAlongNormal = dot(relVel, normal);

                    if (velAlongNormal > 0) continue;

                    float jj = -(1 + bounceFactor) * velAlongNormal;
                    jj /= 1/objects[i].mass + 1/objects[(int)j].mass;

                    float2 impulse = jj * normal;
                    objects[i].velocity -= impulse / objects[i].mass;
                    objects[(int)j].velocity += impulse / objects[(int)j].mass;
                    objects[i].velocity *= 0.99f;  // Small energy dissipation
                    objects[(int)j].velocity *= 0.99f;


                    // Separate overlapping objects
                    float overlap = minDist - dist;
                    float2 separation = overlap * normal * 0.5f;
                    objects[i].position -= separation;
                    objects[(int)j].position += separation;
                }
            }
        }
    }
}
)";
        // Create program from source
        program = cl::Program(context, kernelSource);

        // Build program
        cl_int buildErr = program.build({device});
        if (buildErr != CL_SUCCESS) {
            std::string buildLog = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
            throw std::runtime_error("Error building program: " + buildLog);
        }

        // Create kernels
        cl_int kernelErr;

        // Create updatePositionsAndDetectCollisions kernel
        gravityKernel = cl::Kernel(program, "calculateGravity", &kernelErr);
        if (kernelErr != CL_SUCCESS) {
            throw std::runtime_error("Error creating resolveCollisions kernel: " + std::to_string(kernelErr));
        }

        updatePositionsKernel = cl::Kernel(program, "updatePositionsAndDetectCollisions", &kernelErr);
        if (kernelErr != CL_SUCCESS) {
            throw std::runtime_error("Error creating updatePositionsAndDetectCollisions kernel: " + std::to_string(kernelErr));
        }

        collisionKernel = cl::Kernel(program, "resolveCollisions", &kernelErr);
        if (kernelErr != CL_SUCCESS) {
            throw std::runtime_error("Error creating resolveCollisions kernel: " + std::to_string(kernelErr));
        }

    } catch (const std::exception& e) {
        std::cerr << "Error in initializeOpenCL: " << e.what() << std::endl;
        exit(1);
    }
}
#endif



void updateSimulation(std::vector<Object>& allObjects, SpatialGrid& grid, ThreadPool& pool, double delta_time) {
    #ifdef USE_GPU
    try {
        size_t numObjects = allObjects.size();
        const float worldSize = 2.0f; // Total size is still 2, but ranges from -1 to 1
        const float cellSize = worldSize / 20.0f; // Adjust cell size as needed
        const int gridWidth = static_cast<int>(worldSize / cellSize);
        const int gridSize = gridWidth * gridWidth;

        // Prepare GPU data
        std::vector<Object::GPUData> gpuData(numObjects);
        for (size_t i = 0; i < numObjects; ++i) {
            gpuData[i] = allObjects[i].getGPUData();
        }

        // Create OpenCL buffers
        cl::Buffer objectBuffer(context, CL_MEM_READ_WRITE, sizeof(Object::GPUData) * numObjects);
        cl::Buffer gridCountersBuffer(context, CL_MEM_READ_WRITE, sizeof(int) * gridSize);
        cl::Buffer gridBuffer(context, CL_MEM_READ_WRITE, sizeof(int) * gridSize * 8); // Assuming max 8 objects per cell

        // Copy data to GPU
        queue.enqueueWriteBuffer(objectBuffer, CL_TRUE, 0, sizeof(Object::GPUData) * numObjects, gpuData.data());

        // Clear grid counters
        std::vector<int> zeroCounters(gridSize, 0);
        queue.enqueueWriteBuffer(gridCountersBuffer, CL_TRUE, 0, sizeof(int) * gridSize, zeroCounters.data());

        cl::Kernel gravityKernel(program, "calculateGravity");
        gravityKernel.setArg(0, objectBuffer);
        gravityKernel.setArg(1, static_cast<cl_uint>(numObjects));
        gravityKernel.setArg(2, static_cast<cl_float>(Settings::g_G));
        gravityKernel.setArg(3, static_cast<cl_float>(Settings::g_MAX_FORCE));
        queue.enqueueNDRangeKernel(gravityKernel, cl::NullRange, cl::NDRange(numObjects), cl::NullRange);

        // Update positions and detect collisions
        cl::Kernel updateKernel(program, "updatePositionsAndDetectCollisions");
        updateKernel.setArg(0, objectBuffer);
        updateKernel.setArg(1, static_cast<cl_uint>(numObjects));
        updateKernel.setArg(2, static_cast<cl_float>(delta_time));
        updateKernel.setArg(3, static_cast<cl_float>(worldSize));
        updateKernel.setArg(4, static_cast<cl_float>(cellSize));
        updateKernel.setArg(5, gridCountersBuffer);
        updateKernel.setArg(6, gridBuffer);
        updateKernel.setArg(7, static_cast<cl_int>(gridWidth));

        queue.enqueueNDRangeKernel(updateKernel, cl::NullRange, cl::NDRange(numObjects), cl::NullRange);

        // Resolve collisions
        collisionKernel.setArg(0, objectBuffer);
        collisionKernel.setArg(1, static_cast<cl_uint>(numObjects));
        collisionKernel.setArg(2, static_cast<cl_float>(worldSize));
        collisionKernel.setArg(3, static_cast<cl_float>(cellSize));
        collisionKernel.setArg(4, gridCountersBuffer);
        collisionKernel.setArg(5, gridBuffer);
        collisionKernel.setArg(6, static_cast<cl_int>(gridWidth));
        collisionKernel.setArg(7, static_cast<cl_float>(Settings::g_bounceFactor));

        queue.enqueueNDRangeKernel(collisionKernel, cl::NullRange, cl::NDRange(numObjects), cl::NullRange);

        // Read results back
        queue.enqueueReadBuffer(objectBuffer, CL_TRUE, 0, sizeof(Object::GPUData) * numObjects, gpuData.data());

        // Update CPU-side objects and spatial grid
        grid.clear();
        for (size_t i = 0; i < numObjects; ++i) {
            allObjects[i].updateFromGPUData(gpuData[i]);
            updateObjectTrail(allObjects[i]);
            grid.insert(&allObjects[i]);
        }


    } catch (const std::exception& e) {
        std::cerr << "Error in updateSimulation: " << e.what() << std::endl;
    }
#else
    // CPU-based simulation update
    const size_t numObjects = allObjects.size();
    const size_t numThreads = pool.getThreadCount();
    const size_t batchSize = std::max(size_t(1), numObjects / (numThreads * 4));

    // Step 1: Gravity Calculation
    std::vector<std::future<void>> gravityfutures;
    for (size_t i = 0; i < numObjects; i += batchSize) {
        size_t end = std::min(i + batchSize, numObjects);
        gravityfutures.push_back(pool.enqueue([&allObjects, i, end]() {
            for (size_t j = i; j < end; ++j) {
                calculate_gravity(allObjects[j], allObjects, 0, allObjects.size());
            }
        }));
    }

    // Wait for gravity calculations to complete
    for (auto& future : gravityfutures) {
        future.get();
    }

    // Step 2: Update positions and trails
    std::vector<std::future<void>> updateFutures;
    for (size_t i = 0; i < numObjects; i += batchSize) {
        size_t end = std::min(i + batchSize, numObjects);
        updateFutures.push_back(pool.enqueue([&allObjects, i, end, delta_time]() {
            for (size_t j = i; j < end; ++j) {
                Object& obj = allObjects[j];
                obj.velocity = obj.velocity + obj.acceleration * delta_time;
                obj.position = obj.position + obj.velocity * delta_time;
                updateObjectTrail(obj);
            }
        }));
    }

    // Wait for position updates to complete
    for (auto& future : updateFutures) {
        future.get();
    }

    // Step 3: Rebuild spatial grid
    grid.clear();
    for (auto& obj : allObjects) {
        grid.insert(&obj);
    }

    // Step 4: Handle collisions
    std::vector<std::future<void>> collisionFutures;
    for (size_t i = 0; i < numObjects; i += batchSize) {
        size_t end = std::min(i + batchSize, numObjects);
        collisionFutures.push_back(pool.enqueue([&allObjects, &grid, i, end]() {
            for (size_t j = i; j < end; ++j) {
                auto neighbors = grid.getNeighbors(allObjects[j]);
                for (auto* neighbor : neighbors) {
                    if (checkCollision(allObjects[j], *neighbor)) {
                        handleCollision(allObjects[j], *neighbor);
                    }
                }
            }
        }));
    }

    // Wait for collision handling to complete
    for (auto& future : collisionFutures) {
        future.get();
    }
#endif

}



std::vector<Object> get_objects() {
    std::vector<Object> allObjects;

    // Create and add objects
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 25; j++) {
            const float r = float(i) / (24.0f * 0.9f) + 0.10f;
            const float g = float(j) / (24.0f * 0.9f) + 0.10f;
            const float b = 1.0f - (float(i + j) / (48.0f * 0.9f) + 0.10f);

            allObjects.emplace_back(
                Vector2{-0.5f + i * 0.04f, j * 0.04f},  // position
                Vector2{0.2f, 0.0f},                    // velocity
                1e4,                                    // mass
                0.1,                                    // density
                r, g, b                                 // color
            );
        }
    }

    // Add a larger central object
    allObjects.emplace_back(
        Vector2{0.02f, 0.0001f},  // position
        Vector2{0.2f, 0.0f},    // velocity
        5e14,                   // mass
        5e6,                    // density
        1.0f, 1.0f, 1.0f        // color (white)
    );

    return allObjects;
}