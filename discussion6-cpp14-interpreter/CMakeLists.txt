cmake_minimum_required(VERSION 3.17)
project(SFINAE)

set(CMAKE_CXX_STANDARD 20)
FIND_PACKAGE(Boost)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
add_executable(SFINAE main.cpp)

