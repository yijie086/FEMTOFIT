cmake_minimum_required(VERSION 3.10)
project(CFF_Fitter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- Include Directories ---
include_directories(include)

# --- Source Files ---
file(GLOB SOURCES src/*.cpp)

# --- Add pybind11 as subdirectory instead of using find_package ---
add_subdirectory(external/pybind11)

# --- Define Executable ---
add_executable(CFF_Fitter main.cpp ${SOURCES})

# --- Link pybind11 embedding support ---
target_link_libraries(CFF_Fitter PRIVATE pybind11::embed)
