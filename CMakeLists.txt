cmake_minimum_required(VERSION 3.18.1)

# Include the cmake scripts directory
LIST( APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
)

# create the new project
project(NP_gen)

# Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -DMY_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math")

# Requires OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# install the libraries as well as the executables in lib & bin folder in build/
# This has to be after project is defined
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

add_executable(NP_gen main.cpp basic_struct_operations.cpp load_residue.cpp place_sulfur.cpp)

add_subdirectory(tools)

target_include_directories(NP_gen PUBLIC ${CMAKE_SOURCE_DIR})
