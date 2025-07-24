cmake_minimum_required(VERSION 3.10)
project(CFF_Fitter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Ensure CLAS12ROOT is set
if(NOT DEFINED ENV{CLAS12ROOT})
    message(FATAL_ERROR "CLAS12ROOT environment variable not set!")
endif()

# --- ROOT Configuration ---

# Get ROOT include directory
execute_process(
    COMMAND root-config --incdir
    OUTPUT_VARIABLE ROOT_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get ROOT libraries (includes -lCore, -lMinuit, -ldl, -pthread, etc.)
# ⚠️ NOTE: Do NOT add "Minuit" — it's not a valid argument on this system
execute_process(
    COMMAND root-config --glibs
    OUTPUT_VARIABLE ROOT_LIB_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Convert to list of flags
separate_arguments(ROOT_LIBS UNIX_COMMAND "${ROOT_LIB_FLAGS}")
list(APPEND ROOT_LIBS -lMinuit)
# --- Include Directories ---
include_directories(${ROOT_INCLUDE_DIR})
include_directories($ENV{CLAS12ROOT}/Clas12Root)
include_directories(include)

# --- Library Directories ---
link_directories($ENV{CLAS12ROOT}/lib)

# --- Source Files ---
file(GLOB SOURCES src/*.cpp)

# --- Add pybind11 ---
add_subdirectory(external/pybind11)

# --- Define Executable ---
add_executable(CFF_Fitter main.cpp ${SOURCES})

# --- Link Libraries ---
target_link_libraries(CFF_Fitter PRIVATE
    pybind11::embed
    ${ROOT_LIBS}
    pthread
    Clas12Root
)

# --- Optional Debug Info ---
message(STATUS "CLAS12ROOT: $ENV{CLAS12ROOT}")
message(STATUS "ROOT INCLUDE DIR: ${ROOT_INCLUDE_DIR}")
message(STATUS "ROOT LIBS: ${ROOT_LIBS}")
