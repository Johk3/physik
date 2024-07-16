cmake_minimum_required(VERSION 3.28)
project(physik CXX)

set(CMAKE_CXX_STANDARD 17)

# Find required packages
find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(glad REQUIRED)
find_package(OpenCL REQUIRED)
find_package(OpenGL REQUIRED)


# Function to check for CPU features on macOS
function(check_cpu_feature feature result_var)
    execute_process(
            COMMAND sysctl -n machdep.cpu.features
            OUTPUT_VARIABLE CPU_FEATURES
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    string(FIND "${CPU_FEATURES}" "${feature}" feature_found)
    if(feature_found GREATER -1)
        set(${result_var} TRUE PARENT_SCOPE)
    else()
        set(${result_var} FALSE PARENT_SCOPE)
    endif()
endfunction()

# Check for AVX and AVX2
include(CheckCXXSourceRuns)

if(APPLE)
    target_compile_definitions(physik PRIVATE __APPLE__)
    check_cpu_feature("AVX1.0" CPU_HAS_AVX)
    check_cpu_feature("AVX2" CPU_HAS_AVX2)
    message(STATUS "User is using an APPLE device")
else()
    set(AVX_CODE "
        #include <immintrin.h>
        int main() {
            __m256 a = _mm256_set1_ps(1.0f);
            __m256 b = _mm256_set1_ps(2.0f);
            __m256 c = _mm256_add_ps(a, b);
            float result = _mm256_cvtss_f32(c);
            return result == 3.0f ? 0 : 1;
        }
    ")

    set(AVX2_CODE "
        #include <immintrin.h>
        int main() {
            __m256i a = _mm256_set1_epi32(1);
            __m256i b = _mm256_set1_epi32(2);
            __m256i c = _mm256_add_epi32(a, b);
            return _mm256_extract_epi32(c, 0) == 3 ? 0 : 1;
        }
    ")

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_REQUIRED_FLAGS "-mavx")
    elseif(MSVC)
        set(CMAKE_REQUIRED_FLAGS "/arch:AVX")
    endif()

    check_cxx_source_runs("${AVX_CODE}" CPU_HAS_AVX)

    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_REQUIRED_FLAGS "-mavx2")
    elseif(MSVC)
        set(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
    endif()

    check_cxx_source_runs("${AVX2_CODE}" CPU_HAS_AVX2)
endif()

# ImGui setup
set(IMGUI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/imgui") # Adjust this path as needed

# Add ImGui source files
set(IMGUI_SOURCES
        ${IMGUI_DIR}/imgui.cpp
        ${IMGUI_DIR}/imgui_draw.cpp
        ${IMGUI_DIR}/imgui_tables.cpp
        ${IMGUI_DIR}/imgui_widgets.cpp
        ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
        ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)

# Define header and source files
set(HEADERS
        include/object.h
        include/vector2.h
        include/spatialgrid.h
        include/rendering.h
        include/controlpanel.h
        include/simulationutils.h
        include/settings.h
        include/ThreadPool.h
        include/constants.h
        include/linmath.h
)

set(SOURCES
        main.cpp
        src/object.cpp
        src/rendering.cpp
        src/simulationutils.cpp
        src/controlpanel.cpp
        src/spatialgrid.cpp
        src/settings.cpp
)

# Add executable
add_executable(physik
        ${IMGUI_SOURCES}
        ${HEADERS}
        ${SOURCES}
)

# Include directories
target_include_directories(physik PRIVATE
        ${IMGUI_DIR}
        ${IMGUI_DIR}/backends
)

# Link libraries
target_link_libraries(physik
        PRIVATE
        glfw
        GLEW::GLEW
        glad::glad
        OpenCL::OpenCL
        OpenGL::GL
)

# Set compiler flags and definitions based on CPU support
if(CPU_HAS_AVX2)
    if(MSVC)
        target_compile_options(physik PRIVATE /O2 /arch:AVX2)
    else()
        target_compile_options(physik PRIVATE -O3 -march=native -mavx2 -mfma)
    endif()
    target_compile_definitions(physik PRIVATE USE_AVX2)
    message(STATUS "AVX2 support detected and enabled")
elseif(CPU_HAS_AVX)
    if(MSVC)
        target_compile_options(physik PRIVATE /O2 /arch:AVX)
    else()
        target_compile_options(physik PRIVATE -O3 -march=native -mavx)
    endif()
    target_compile_definitions(physik PRIVATE USE_AVX)
    message(STATUS "AVX support detected and enabled")
else()
    message(STATUS "Neither AVX2 nor AVX support detected, using standard optimizations")
    if(MSVC)
        target_compile_options(physik PRIVATE /O2)
    else()
        target_compile_options(physik PRIVATE -O3)
    endif()
endif()

# ImGui compile definitions
target_compile_definitions(physik PRIVATE IMGUI_IMPL_OPENGL_LOADER_GLAD)
