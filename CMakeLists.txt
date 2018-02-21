cmake_minimum_required(VERSION 2.8.12)
project(htt)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -O2")
include_directories("include")


set(DELPHES_LIBRARIES "$ENV{DELPHES_ROOT}/libDelphes.so")
include_directories("$ENV{DELPHES_ROOT}" "$ENV{DELPHES_ROOT}/external")

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

find_package(Boost 1.63 COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})


add_subdirectory(src)
add_subdirectory(prog)
