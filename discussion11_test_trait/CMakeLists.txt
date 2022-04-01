cmake_minimum_required(VERSION 3.21)
project(test_trait)

set(CMAKE_CXX_STANDARD 20)

add_executable(test_trait main.cpp)
include(FindRust.cmake)
include(CMakeCargo.cmake)
include(FetchContent)

macro(find_fetch package repo)
    find_package(${package} QUIET)
    if(NOT ${${package}_FOUND})
        message(STATUS "fetching ${package} from ${repo}")
        FetchContent_Declare(${package}
                GIT_REPOSITORY ${repo}
                GIT_SHALLOW true
                )
        FetchContent_MakeAvailable(${package})
        FetchContent_GetProperties(${package})
        if(NOT ${${package}_POPULATED})
            FetchContent_Populate(${package})
            add_subdirectory(${${package}_SOURCE_DIR} ${${package}_BINARY_DIR})
        endif(NOT ${${package}_POPULATED})
        list(APPEND CMAKE_PREFIX_PATH ${${package}_BINARY_DIR})
        find_package(${package} REQUIRED)
    endif(NOT ${${package}_FOUND})
endmacro()

find_fetch(xtl https://github.com/xtensor-stack/xtl)
find_fetch(xtensor https://github.com/xtensor-stack/xtensor)

if (NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/develop/conan.cmake"
            "${CMAKE_BINARY_DIR}/conan.cmake"
            TLS_VERIFY ON)
endif ()

include(${CMAKE_BINARY_DIR}/conan.cmake)

conan_cmake_configure(REQUIRES fmt/8.1.1
        GENERATORS cmake_find_package)

conan_cmake_configure(REQUIRES xtensor/0.23.10
        GENERATORS cmake_find_package)

conan_cmake_autodetect(settings)

conan_cmake_install(PATH_OR_REFERENCE .
        BUILD missing
        REMOTE conancenter
        SETTINGS ${settings})
find_package(fmt)
find_package(xtensor)
target_link_libraries(test_trait xtensor fmt::fmt)

cargo_build(NAME test_trait_rs)
