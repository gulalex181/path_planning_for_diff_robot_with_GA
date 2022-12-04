cmake_minimum_required(VERSION 3.14)

# set the project name and version
project(object_control_with_genetic_algorithm)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# include_directories(
#     include
# )

add_subdirectory(matplotplusplus)

# add the executable
add_executable(object_control
    object_control.cpp
    genetic_algorithm.cpp
    gradient_descent_method.cpp
    graphics.h
)
target_link_libraries(object_control PUBLIC
    matplot
)