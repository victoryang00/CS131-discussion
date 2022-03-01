cmake_minimum_required(VERSION 3.16)

project(lambda_test)

set(CMAKE_CXX_STANDARD 20)

#find_package(OpenCV REQUIRED CONFIG)
#find_package(benchmark REQUIRED)
add_executable(lambda_test main.cpp)
#include_directories(${CMAKE_SOURCE_DIR}/include
#        ${OpenCV_INCLUDE_DIRS})
add_subdirectory(benchmark)
target_link_libraries(lambda_test benchmark::benchmark)