# ==============================================================================
# Project Configuration
# ==============================================================================
cmake_minimum_required (VERSION 3.8)

# Set project name
project ("MaxPR")

# Set C++17 standard (as specified in c_cpp_properties)
set (CMAKE_CXX_STANDARD 17)

# ==============================================================================
# Include Directories
# ==============================================================================
# Declare directory containing .h files
include_directories ("${CMAKE_SOURCE_DIR}/include")

# ==============================================================================
# Source Files
# ==============================================================================
# Find all .cpp files
file(GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/*.cpp")

# ==============================================================================
# Build Target
# ==============================================================================
# Create executable
add_executable (main ${SOURCE_FILES})

# ==============================================================================
# Set working directory for Visual Studio (Debug & Release)
# ==============================================================================
set_target_properties(main PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)