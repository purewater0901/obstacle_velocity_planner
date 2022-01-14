cmake_minimum_required(VERSION 3.21)
project(obstacle_velocity_planner)

#set(CMAKE_CXX_STANDARD 14)
add_compile_options(-std=c++14 -Wall)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})
find_package(Python3 COMPONENTS Development NumPy)

find_package(osqp REQUIRED)

include_directories(include)
include_directories(${Boost_INCLUDE_DIR})

add_executable(obstacle_velocity_planner main.cpp
                                         src/interpolation/linear_interpolation.cpp
                                         src/interpolation/spline_interpolation.cpp
                                         src/csc_matrix_conv.cpp
                                         src/osqp_interface.cpp
                                         src/velocity_optimizer.cpp
                                         src/velocity_smoother.cpp)
target_include_directories(obstacle_velocity_planner PRIVATE ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(obstacle_velocity_planner Python3::Python Python3::NumPy osqp::osqp)
