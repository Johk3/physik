cmake_minimum_required(VERSION 3.28)
project(physik)

set(CMAKE_CXX_STANDARD 17)

# Find required packages
find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(glad REQUIRED)

# SIMD
if(MSVC)
    add_compile_options(/O2 /arch:AVX2)
else()
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

# Add executable
add_executable(physik
        main.cpp
        linmath.h
        settings.h
        ThreadPool.h
        ${IMGUI_SOURCES}
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

# ImGui compile definitions (if needed)
target_compile_definitions(physik PRIVATE IMGUI_IMPL_OPENGL_LOADER_GLAD)