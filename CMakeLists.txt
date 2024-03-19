cmake_minimum_required(VERSION 3.10)

project(space)

# Export compile commands for VSCode
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# GraphBLAS
set(GRAPHBLAS_DIR ${CMAKE_SOURCE_DIR}/external/GraphBLAS)
set(GRAPHBLAS_INCLUDE_DIR ${GRAPHBLAS_DIR}/Include)
set(GRAPHBLAS_LIBRARY_DIR ${GRAPHBLAS_DIR}/build)
find_library(GRAPHBLAS_LIB graphblas PATHS ${GRAPHBLAS_LIBRARY_DIR})
include_directories(${GRAPHBLAS_INCLUDE_DIR})

# Print the GraphBLAS Paths
message(STATUS "GRAPHBLAS_DIR: ${GRAPHBLAS_DIR}")
message(STATUS "GRAPHBLAS_INCLUDE_DIR: ${GRAPHBLAS_INCLUDE_DIR}")
message(STATUS "GRAPHBLAS_LIBRARY_DIR: ${GRAPHBLAS_LIBRARY_DIR}")
message(STATUS "GRAPHBLAS_LIB: ${GRAPHBLAS_LIB}")

# don't omit frame pointers for better perfomance analysis
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")

# Include
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
# get all the headers
file(GLOB_RECURSE HEADERS ${INCLUDE_DIR}/*.h)

include_directories(include)

add_subdirectory(tests)
