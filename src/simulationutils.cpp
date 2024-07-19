#include "../include/simulationutils.h"
#include <iostream>
#include "../include/settings.h"
#ifndef USE_GPU
    #include <cpuid.h>
#endif
#include <fstream>
#ifdef USE_AVX2
#include <immintrin.h>
#endif
#ifdef USE_GPU
#include <CL/opencl.hpp>
// Global OpenCL variables
cl::Context context;
cl::CommandQueue queue;
cl::Kernel kernel;
cl::Program program;
cl::Buffer d_positions;
cl::Buffer d_masses;
cl::Buffer d_accelerations;

// Function to initialize OpenCL
void initOpenCL() {
    // Get available platforms
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);

    if (platforms.empty()) {
        std::cerr << "No OpenCL platforms found" << std::endl;
        exit(1);
    }

    // Select the first platform
    cl::Platform platform = platforms[0];

    // Get available devices
    std::vector<cl::Device> devices;
    platform.getDevices(CL_DEVICE_TYPE_GPU, &devices);

    if (devices.empty()) {
        std::cerr << "No OpenCL devices found" << std::endl;
        exit(1);
    }

    // Select the first device
    cl::Device device = devices[0];

    // Create a context
    context = cl::Context({device});

    // Create a command queue
    queue = cl::CommandQueue(context, device);

    // Read the kernel source
    const char* kernelSource = R"(
__kernel void calculate_gravity(
    __global const float4* positions,
    __global const float* masses,
    __global float4* accelerations,
    const float G,
    const float epsilon,
    const float max_force,
    const int num_objects
) {
    int gid = get_global_id(0);
    if (gid >= num_objects) return;

    float4 pos1 = positions[gid];
    float4 acc = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

    for (int i = 0; i < num_objects; i++) {
        if (i == gid) continue;

        float4 pos2 = positions[i];
        float4 r = pos2 - pos1;
        float dist2 = r.x*r.x + r.y*r.y + r.z*r.z;

        if (dist2 > epsilon) {
            float inv_dist = rsqrt(dist2);
            float inv_dist3 = inv_dist * inv_dist * inv_dist;
            float force = min(G * masses[i] * inv_dist3, max_force);
            acc += force * r;
        }
    }

    accelerations[gid] = acc;
}
)";
    // Create the program
    program = cl::Program(context, kernelSource);

    // Build the program
    if (program.build({device}) != CL_SUCCESS) {
        std::cerr << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device) << std::endl;
        exit(1);
    }

    // Create the kernel
    kernel = cl::Kernel(program, "calculate_gravity");
}

void calculate_gravity_gpu(std::vector<Object>& objects) {
    int num_objects = objects.size();

    // Prepare data for GPU
    std::vector<cl_float4> positions(num_objects);
    std::vector<float> masses(num_objects);
    std::vector<cl_float4> accelerations(num_objects);

    for (int i = 0; i < num_objects; i++) {
        positions[i] = {(float)objects[i].position.x, (float)objects[i].position.y, (float)objects[i].position.z, 0.0f};
        masses[i] = (float)objects[i].mass;
    }

    // Create buffers
    d_positions = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_float4) * num_objects, positions.data());
    d_masses = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * num_objects, masses.data());
    d_accelerations = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_float4) * num_objects);

    // Set kernel arguments
    kernel.setArg(0, d_positions);
    kernel.setArg(1, d_masses);
    kernel.setArg(2, d_accelerations);
    kernel.setArg(3, (float)Settings::g_G);
    kernel.setArg(4, (float)Settings::EPSILON);
    kernel.setArg(5, (float)Settings::g_MAX_FORCE);
    kernel.setArg(6, num_objects);

    // Execute the kernel
    queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(num_objects), cl::NullRange);

    // Read the results
    queue.enqueueReadBuffer(d_accelerations, CL_TRUE, 0, sizeof(cl_float4) * num_objects, accelerations.data());

    // Update object accelerations
    for (int i = 0; i < num_objects; i++) {
        objects[i].acceleration = Vector3{accelerations[i].s[0], accelerations[i].s[1], accelerations[i].s[2]};
    }
}
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

// SingleInstructionMultipleData Gravity calculation using AVX architecture
void calculate_gravity_simd(Object& object1, const std::vector<Object>& objects, size_t start, size_t end) {
    // Load object1's position into AVX registers
    const __m256d obj1_x = _mm256_set1_pd(object1.position.x);
    const __m256d obj1_y = _mm256_set1_pd(object1.position.y);
    const __m256d obj1_z = _mm256_set1_pd(object1.position.z);  // New Z component

    // Load constants into AVX registers
    const __m256d epsilon = _mm256_set1_pd(Settings::EPSILON);
    const __m256d max_force = _mm256_set1_pd(Settings::g_MAX_FORCE);
    const __m256d g_const = _mm256_set1_pd(Settings::g_G);

    // Initialize acceleration accumulators
    __m256d acc_x = _mm256_setzero_pd();
    __m256d acc_y = _mm256_setzero_pd();
    __m256d acc_z = _mm256_setzero_pd();  // New Z component

    // Process objects in batches of 4
    for (size_t i = start; i < end; i += 4) {
        __m256d obj2_x, obj2_y, obj2_z, mass;

        if (i + 3 < end) {
            // Load 4 objects' data into AVX registers
            obj2_x = _mm256_setr_pd(objects[i].position.x, objects[i+1].position.x, objects[i+2].position.x, objects[i+3].position.x);
            obj2_y = _mm256_setr_pd(objects[i].position.y, objects[i+1].position.y, objects[i+2].position.y, objects[i+3].position.y);
            obj2_z = _mm256_setr_pd(objects[i].position.z, objects[i+1].position.z, objects[i+2].position.z, objects[i+3].position.z);
            mass = _mm256_setr_pd(objects[i].mass, objects[i+1].mass, objects[i+2].mass, objects[i+3].mass);
        } else {
            // Handle the last incomplete batch
            double x_array[4] = {0}, y_array[4] = {0}, z_array[4] = {0}, m_array[4] = {0};
            for (size_t j = 0; j < end - i; ++j) {
                x_array[j] = objects[i+j].position.x;
                y_array[j] = objects[i+j].position.y;
                z_array[j] = objects[i+j].position.z;
                m_array[j] = objects[i+j].mass;
            }
            obj2_x = _mm256_loadu_pd(x_array);
            obj2_y = _mm256_loadu_pd(y_array);
            obj2_z = _mm256_loadu_pd(z_array);
            mass = _mm256_loadu_pd(m_array);
        }

        // Calculate distance vectors
        __m256d dx = _mm256_sub_pd(obj2_x, obj1_x);
        __m256d dy = _mm256_sub_pd(obj2_y, obj1_y);
        __m256d dz = _mm256_sub_pd(obj2_z, obj1_z);

        // Calculate squared distance
        __m256d r2 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)), _mm256_mul_pd(dz, dz));

        // Calculate distance and cube of distance
        __m256d r = _mm256_sqrt_pd(r2);
        __m256d r3 = _mm256_mul_pd(r, _mm256_mul_pd(r, r));

        // Calculate gravitational force
        __m256d force = _mm256_div_pd(_mm256_mul_pd(g_const, mass), r3);
        force = _mm256_min_pd(force, max_force);

        // Apply epsilon check to avoid division by zero
        __m256d mask = _mm256_cmp_pd(r2, epsilon, _CMP_GT_OQ);
        force = _mm256_and_pd(force, mask);

        // Accumulate acceleration components
        acc_x = _mm256_add_pd(acc_x, _mm256_mul_pd(force, dx));
        acc_y = _mm256_add_pd(acc_y, _mm256_mul_pd(force, dy));
        acc_z = _mm256_add_pd(acc_z, _mm256_mul_pd(force, dz));
    }

    // Sum up the vector components
    double acc_x_sum = _mm256_reduce_add_pd(acc_x);
    double acc_y_sum = _mm256_reduce_add_pd(acc_y);
    double acc_z_sum = _mm256_reduce_add_pd(acc_z);

    // Set the final acceleration
    object1.acceleration = Vector3{acc_x_sum, acc_y_sum, acc_z_sum};
}
#endif
// Non-AVX version of gravity calculation

void calculate_gravity_normal(Object& object1, const std::vector<Object>& objects, size_t start, size_t end) {
    Vector3 acc = {0.0, 0.0, 0.0};
    for (size_t i = start; i < end; ++i) {
        const auto& object2 = objects[i];
        if (&object1 == &object2) continue;
        if (object2.mass < Settings::EPSILON) continue;

        Vector3 diff = object2.position - object1.position;
        double r2 = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;

        if (r2 < 1e-12) continue;

        double r = std::sqrt(r2);
        if (r < 1e-6f) continue;

        double force = Settings::g_G * object2.mass / (r2 * r);
        force = std::min(force, Settings::g_MAX_FORCE);

        acc = acc + diff * force;
    }
    object1.acceleration = acc;
}

void updateObjectTrail(Object& obj) {
    if (obj.trail.empty()) {
        obj.trail.push_back(obj.position);
        return;
    }

    const Vector3 lastPosition = obj.trail.back();
    const Vector3 diff = obj.position - lastPosition;
    double distance = diff.length();

    if (distance >= obj.TRAIL_SPACING) {
        // Calculate how many new dots we need to add
        const int numNewDots = static_cast<int>(distance / obj.TRAIL_SPACING);

        for (int i = 0; i < numNewDots; ++i) {
            const double t = (i + 1) * obj.TRAIL_SPACING / distance;
            Vector3 newDotPosition = {
                lastPosition.x + diff.x * t,
                lastPosition.y + diff.y * t,
                lastPosition.z + diff.z * t
            };
            obj.trail.push_back(newDotPosition);

            // Remove oldest dot if we exceed the maximum trail length
            if (obj.trail.size() > obj.MAX_TRAIL_LENGTH) {
                obj.trail.erase(obj.trail.begin());
            }
        }
    }
}

// Helper function to get the bounding sphere of an object
std::pair<Vector3, double> getBoundingSphere(const Object& obj) {
    Vector3 center = obj.position;
    double radius = obj.radius;

    switch(obj.shape) {
        case ObjectShape::SPHERE:
                break;
        case ObjectShape::TRIANGLE: {
            double height = obj.radius * std::sqrt(3) / 2;
            double depth = obj.radius / 4;
            center.y += height / 3;
            radius = std::sqrt(std::pow(2*height/3, 2) + std::pow(obj.radius/2, 2) + std::pow(depth/2, 2));
            break;
        }
        case ObjectShape::FLAT_SURFACE:
                    center.y += obj.radius / 4;
        // Radius is half the diagonal of the prism
        radius = std::sqrt(std::pow(obj.radius * 2, 2) + std::pow(obj.radius * 2, 2) + std::pow(obj.radius, 2)) / 2;
        break;
        case ObjectShape::CONTAINER:
                center.y += obj.radius;
        radius = std::sqrt(3) * obj.radius;
        break;
        case ObjectShape::COW:
            // Increase the radius to encompass the cow shape
                radius *= 2;
        break;
    }

    return {center, radius};
}

// Check if two objects are overlapping using their bounding spheres
bool checkCollision(const Object& obj1, const Object& obj2) {
    auto [center1, radius1] = getBoundingSphere(obj1);
    auto [center2, radius2] = getBoundingSphere(obj2);

    const double distance = (center1 - center2).length();
    return distance <= (radius1 + radius2);
}

void handleCollision(Object& obj1, Object& obj2) {
    auto [center1, radius1] = getBoundingSphere(obj1);
    auto [center2, radius2] = getBoundingSphere(obj2);

    // Calculate vector from obj1 to obj2
    Vector3 delta = center2 - center1;
    const float distance = delta.length();

    // Calculate overlap
    const float overlap = (radius1 + radius2) - distance;

    if (overlap > 0) {
        // Normalize the delta vector
        const Vector3 normal = delta.normal();

        // Separate the objects
        float mass_ratio1 = obj2.mass / (obj1.mass + obj2.mass);
        float mass_ratio2 = obj1.mass / (obj1.mass + obj2.mass);

        obj1.position = obj1.position - normal * (overlap * mass_ratio1);
        obj2.position = obj2.position + normal * (overlap * mass_ratio2);

        // Calculate relative velocity
        const Vector3 relativeVelocity = obj2.velocity - obj1.velocity;

        // Calculate relative velocity along the normal
        float velocityAlongNormal = relativeVelocity.dot(normal);

        // Do not resolve if velocities are separating
        if (velocityAlongNormal > 0)
            return;

        // Calculate impulse scalar
        float impulseScalar = -(1 + Settings::g_bounceFactor) * velocityAlongNormal;
        impulseScalar /= 1 / obj1.mass + 1 / obj2.mass;

        // Apply impulse
        Vector3 impulse = normal * impulseScalar;

        // Update linear velocities
        obj1.velocity = obj1.velocity - impulse * (1 / obj1.mass);
        obj2.velocity = obj2.velocity + impulse * (1 / obj2.mass);

        // Calculate point of impact
        Vector3 impactPoint = center1 + normal * radius1;

        // Only apply angular impulses if rotation is enabled
        if (Settings::g_enableRotation) {
            // Calculate angular impulse
            Vector3 r1 = impactPoint - obj1.position;
            Vector3 r2 = impactPoint - obj2.position;

            Vector3 angularImpulse1 = r1.cross(impulse);
            Vector3 angularImpulse2 = r2.cross(impulse);

            // Calculate impact speed
            float impactSpeed = std::abs(velocityAlongNormal);

            // Scale angular impulse based on impact speed
            float angularImpulseScale = std::min(impactSpeed / 10.0f, 1.0f);  // Cap at 1.0
            angularImpulse1 = angularImpulse1 * angularImpulseScale;
            angularImpulse2 = angularImpulse2 * angularImpulseScale;

            // Update angular velocities
            obj1.angularVelocity = obj1.angularVelocity + angularImpulse1 * (1 / obj1.momentOfInertia);
            obj2.angularVelocity = obj2.angularVelocity - angularImpulse2 * (1 / obj2.momentOfInertia);
        }
    }
}


void updateSimulation(std::vector<Object>& allObjects, SpatialGrid& grid, ThreadPool& pool, double delta_time) {
    const size_t numObjects = allObjects.size();
    const size_t numThreads = pool.getThreadCount();
    const size_t batchSize = std::max(size_t(1), numObjects / (numThreads * 4));

    // Step 1: Gravity Calculation
    #ifdef USE_GPU
        calculate_gravity_gpu(allObjects);
    #else
        std::vector<std::future<void>> gravityfutures;
        for (size_t i = 0; i < numObjects; i += batchSize) {
            size_t end = std::min(i + batchSize, numObjects);
            gravityfutures.push_back(pool.enqueue([&allObjects, i, end]() {
                for (size_t j = i; j < end; ++j) {
                    #ifdef USE_AVX2
                        calculate_gravity_simd(allObjects[j], allObjects, 0, allObjects.size());
                    #else
                        calculate_gravity_normal(allObjects[j], allObjects, 0, allObjects.size());
                    #endif
                }
            }));
        }

        for (auto& future : gravityfutures) {
            future.get();
        }
    #endif

    // Step 2: Update positions and trails
    std::vector<std::future<void>> updateFutures;
    for (size_t i = 0; i < numObjects; i += batchSize) {
        size_t end = std::min(i + batchSize, numObjects);
        updateFutures.push_back(pool.enqueue([&allObjects, i, end, delta_time]() {
            for (size_t j = i; j < end; ++j) {
                Object& obj = allObjects[j];
                obj.velocity = obj.velocity + obj.acceleration * delta_time;
                obj.position = obj.position + obj.velocity * delta_time;
                if (Settings::g_enableRotation) {// Only update rotation if enabled
                    obj.updateRotation(delta_time);
                }
                updateObjectTrail(obj);
            }
        }));
    }

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
                    if (&allObjects[j] != neighbor && checkCollision(allObjects[j], *neighbor)) {
                        handleCollision(allObjects[j], *neighbor);
                    }
                }
            }
        }));
    }

    for (auto& future : collisionFutures) {
        future.get();
    }
}

std::vector<Object> get_objects() {
    std::vector<Object> allObjects;
    allObjects.emplace_back(Vector3{0.0, 0.0, 0.0}, Vector3{0.0f, 0.0f, 0.0f}, 5e11, 5e1, 1.0, 1.0, 1.0, ObjectShape::SPHERE);

    // Add examples of new shapes
    allObjects.emplace_back(Vector3{0.5f, 0.0f, 0.0f}, Vector3{0.0f, 0.0f, 0.0f}, 1e3, 0.00001, 1.0, 0.0, 0.0, ObjectShape::TRIANGLE);
    allObjects.emplace_back(Vector3{-0.5f, -0.0f, 0.0f}, Vector3{0.0f, 0.0f, 0.0f}, 1e3, 0.00001, 0.0, 1.0, 0.0, ObjectShape::FLAT_SURFACE);
    allObjects.emplace_back(Vector3{0.0f, -0.5f, -0.5f}, Vector3{0.0f, 0.0f, 0.0f}, 1e3, 0.00001, 0.0, 0.0, 1.0, ObjectShape::CONTAINER);
    allObjects.emplace_back(Vector3{-0.0f, -0.5f, -0.0f}, Vector3{0.0f, 0.0f, 0.0f}, 1e3, 0.00001, 1.0, 1.0, 0.0, ObjectShape::COW);

    return allObjects;
}