cmake_minimum_required(VERSION 3.28)
project(physik)

set(CMAKE_CXX_STANDARD 17)
find_package(glfw3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(glad REQUIRED)



add_executable(physik main.cpp
        linmath.h)

target_link_libraries(physik glfw)
target_link_libraries(physik GLEW::GLEW)
target_link_libraries(physik glad::glad)


