cmake_minimum_required(VERSION 3.16)

project("TREX - Scherman-Morrison random generator" VERSION "0.1.0")

set(INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/headers)
set(GLOBAL_BIN_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(GLOBAL_LIB_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Handle versioning
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/versioning.h.in" "${PROJECT_BINARY_DIR}/versioning.h")
include_directories(${PROJECT_BINARY_DIR})

find_package(HDF5 REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(OpenMP REQUIRED)

add_subdirectory(src)