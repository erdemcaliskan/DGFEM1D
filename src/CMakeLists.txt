cmake_minimum_required(VERSION 3.10) 

project(solver VERSION 1.0)

include(CTest) 
include(CPack) 
enable_testing()

aux_source_directory(util DIR_UTIL_SRCS)

add_executable(solver main.cpp ${DIR_UTIL_SRCS})

target_link_libraries(solver PUBLIC Eigen3::Eigen)