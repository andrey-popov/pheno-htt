cmake_minimum_required(VERSION 2.8.12)
project(htt)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


include_directories("include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wextra -O2")


# Sources for the main library
set(LIBRARY_SOURCES "")
foreach(file AnalysisPlugin.cpp DelphesReader.cpp LJetsSelection.cpp NuReco.cpp Plugin.cpp
  Processor.cpp)
    list(APPEND LIBRARY_SOURCES "${PROJECT_SOURCE_DIR}/src/${file}")
endforeach(file)


# External dependencies
set(DELPHES_LIBRARIES "$ENV{DELPHES_ROOT}/libDelphes.so")
include_directories("$ENV{DELPHES_ROOT}" "$ENV{DELPHES_ROOT}/external")

find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIRS})

find_package(Boost 1.63.0 COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})


# Main library
add_library(htt SHARED ${LIBRARY_SOURCES})
target_link_libraries(htt ${DELPHES_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY})

# Main program
add_executable(htt-tuples "prog/htt-tuples.cpp")
target_link_libraries(htt-tuples htt)
