cmake_minimum_required(VERSION 3.16)
project(Matrix_dynamics)

set(CMAKE_CXX_STANDARD 20)

find_package(PythonLibs 2.7)

include_directories(${PYTHON_INCLUDE_DIRS})

add_executable(Matrix_dynamics main.cpp dynamics.h matplotlibcpp.h traj_generate.h traj_generate.cpp)

target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})