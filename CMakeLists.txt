cmake_minimum_required(VERSION 3.10)
project(DGFEM1D VERSION 0.1.0
        DESCRIPTION "One dimensional Discontinuous Galerkin Finite Element Method Solver"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/.)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

find_package (Eigen3 3.3 REQUIRED)

add_subdirectory(src)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

