cmake_minimum_required(VERSION 3.28)
project(physik CXX)

set(CMAKE_CXX_STANDARD 17)


# Find required packages
find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(glad REQUIRED)


# Check AVX2 Support
include(CheckCXXSourceRuns)
set(AVX2_CODE "
    #include <immintrin.h>
    int main() {
        __m256i a = _mm256_set1_epi32(1);
        __m256i b = _mm256_set1_epi32(2);
        __m256i c = _mm256_add_epi32(a, b);
        return _mm256_extract_epi32(c, 0) == 3 ? 0 : 1;
    }
")
include(CheckCXXSourceRuns)
include(CheckCXXCompilerFlag)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    set(CMAKE_REQUIRED_FLAGS "-mavx2")
elseif(MSVC)
    set(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
endif()

check_cxx_source_runs("${AVX2_CODE}" COMPILER_SUPPORTS_AVX2)

if(MSVC AND COMPILER_SUPPORTS_AVX2)
    message(STATUS "AVX support detected, enabling AVX MSVC")
    add_compile_options(/O2 /arch:AVX2)
endif ()
if(NOT MSVC AND COMPILER_SUPPORTS_AVX2)
    message(STATUS "AVX support detected, enabling AVX MINGW")

    add_compile_options(-O3 -march=native -mavx2 -mfma)
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
)

# AVX2 support detection
if(COMPILER_SUPPORTS_AVX2)
    target_compile_definitions(physik PRIVATE USE_AVX2)
    message(STATUS "AVX support detected, enabling AVX optimizations")
else()
    message(STATUS "AVX support not detected, using standard optimizations")
    # These optimization flags can be enabled
#    if(MSVC)
#        add_compile_options(/O2)
#    else()
#        add_compile_options(-O3)
#    endif()
endif()

# ImGui compile definitions (if needed)
target_compile_definitions(physik PRIVATE IMGUI_IMPL_OPENGL_LOADER_GLAD)